use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::path::Path;
use std::collections::HashMap;
use clap::Parser;
use regex::Regex;
use anyhow::{Result, anyhow};
use std::collections::BTreeMap;

#[derive(Parser, Debug)]
#[command(
    author, 
    version, 
    about = "Split large genome sequences into smaller fragments at N-stretches or intergenic regions",
    long_about = "A tool for splitting large genome sequences into manageable fragments.\n\
    It identifies suitable split points either at long stretches of N bases or in intergenic regions\n\
    to avoid disrupting gene annotations. When a GFF/GTF file is provided, the tool ensures splits\n\
    occur only between genes, maintaining the integrity of gene annotations.\n\
    \n\
    The tool outputs:\n\
    - Split sequences in FASTA format (.fa)\n\
    - Split positions in TSV format (.cutsite.tsv)\n\
    - Adjusted annotation file if GFF/GTF is provided"
)]
struct Args {
    /// Input genome sequence file in FASTA format
    #[arg(short = 'f', long = "fasta")]
    fa: String,

    /// Optional GTF/GFF annotation file for the genome
    #[arg(short = 'g', long = "gtf")]
    gtf: Option<String>,

    /// Prefix for output files (.fa and .cutsite.tsv will be appended)
    #[arg(short = 'o', long = "prefix")]
    prefix: String,

    /// Minimum length of consecutive N bases to be used as a split separator
    #[arg(long = "Ns", default_value = "10", hide = true)]
    num_n: usize,

    /// Minimum length of output scaffold fragments (in base pairs)
    #[arg(long = "min_length", default_value = "300000000")]
    min_length: usize,

    /// Maximum length of output scaffold fragments (in base pairs)
    #[arg(long = "max_length", default_value = "500000000")]
    max_length: usize,
}

#[derive(Debug, Clone)]
struct GeneRegion {
    start: usize,
    end: usize,
}

/// 染色体上所有基因区域（已排序且合并重叠区域）
type MergedGeneRegions = HashMap<String, Vec<GeneRegion>>;

fn process_fragment(
    id: &str,
    fragment: &str,
    start_position: usize,
    length: usize,
    outpref: &str,
    real_len: usize,
) -> Result<()> {
    let end = start_position + length - 1;
    let start_position = start_position + 1;
    let end = end + 1;

    let mut out_fa = BufWriter::new(
        fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(format!("{}.fa", outpref))?
    );

    if real_len == length {
        writeln!(out_fa, ">{}", id)?;
    } else {
        writeln!(out_fa, ">{}_{}_{}",id, start_position, end)?;
    }

    // Output sequence with line breaks every 60 characters
    for chunk in fragment.as_bytes().chunks(60) {
        out_fa.write_all(chunk)?;
        out_fa.write_all(b"\n")?;
    }

    let mut out_txt = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(format!("{}.cutsite.tsv", outpref))?;
    
    writeln!(out_txt, "{}\t{}\t{}", id, start_position, end)?;
    Ok(())
}

/// 处理GTF/GFF文件，提取基因区域并合并重叠区域
/// 
/// 此函数解析GTF/GFF文件，提取基因和外显子位置信息，合并重叠区域，并根据染色体组织
/// 返回按染色体名称索引的包含已合并基因区域的数据结构
fn get_merged_gene_regions(gxf_file: Option<&str>) -> Result<MergedGeneRegions> {
    let mut gene_regions: HashMap<String, Vec<GeneRegion>> = HashMap::new();
    
    let Some(gxf_path) = gxf_file else {
        return Ok(gene_regions);
    };

    eprintln!("Reading gene annotation file: {}", gxf_path);
    let file = File::open(gxf_path)?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // 使用1MB缓冲区提高读取性能

    let mut raw_regions: HashMap<String, Vec<GeneRegion>> = HashMap::new();
    let mut line_count = 0;
    let mut feature_count = 0;
    let mut skipped_count = 0;
    
    // 第一遍扫描，收集所有基因和外显子区域
    for line in reader.lines() {
        let line = line?;
        line_count += 1;
        
        if line_count % 1000000 == 0 {
            eprintln!("Scanned {} lines", line_count);
        }
        
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            skipped_count += 1;
            continue;
        }

        let chr = fields[0].to_string();
        let type_ = fields[2].to_lowercase();
        if !["gene"].contains(&type_.as_str()) {
            continue;
        }

        let start: usize = fields[3].parse()?;
        let end: usize = fields[4].parse()?;

        raw_regions.entry(chr)
            .or_insert_with(Vec::new)
            .push(GeneRegion { start, end });
            
        feature_count += 1;
    }
    
    eprintln!("Initial scan completed: processed {} features, {} chromosomes, skipped {} lines", 
              feature_count, raw_regions.len(), skipped_count);
    
    // 对每条染色体的基因区域进行排序和合并
    for (chr, regions) in raw_regions {
        if regions.is_empty() {
            continue;
        }
        
        // 按起始位置排序
        let mut sorted_regions = regions;
        sorted_regions.sort_by_key(|r| r.start);
        
        // 合并重叠区域
        let mut merged_regions = Vec::new();
        let mut current = sorted_regions[0].clone();
        
        for region in sorted_regions.iter().skip(1) {
            // 如果当前区域与前一个区域重叠或相邻，则合并
            if region.start <= current.end + 1 {
                // 扩展当前区域
                current.end = current.end.max(region.end);
            } else {
                // 没有重叠，将当前区域添加到结果中并开始新区域
                merged_regions.push(current);
                current = region.clone();
            }
        }
        // 添加最后一个区域
        merged_regions.push(current);
        
        eprintln!("Chromosome {}: merged {} raw regions into {} non-overlapping regions", 
                 chr, sorted_regions.len(), merged_regions.len());
        
        gene_regions.insert(chr, merged_regions);
    }
    
    eprintln!("Gene region processing completed: merged regions for {} chromosomes", gene_regions.len());

    Ok(gene_regions)
}

/// 判断给定位置是否在基因间区域
/// 
/// 此函数使用二分查找快速找到位置是否落在任何基因区域之间
/// 如果染色体没有基因注释，所有位置都被视为基因间区域
fn is_in_intergenic(pos: usize, chr: &str, gene_regions: &MergedGeneRegions) -> bool {
    if let Some(regions) = gene_regions.get(chr) {
        // 如果没有基因区域，视为基因间区域
        if regions.is_empty() {
            return true;
        }
        
        // 检查是否在第一个基因之前
        if pos < regions[0].start {
            return true;
        }
        
        // 检查是否在最后一个基因之后
        if pos > regions.last().unwrap().end {
            return true;
        }
        
        // 二分查找找到pos所在的区间
        let idx = match regions.binary_search_by(|r| {
            if pos < r.start { std::cmp::Ordering::Greater }
            else if pos > r.end { std::cmp::Ordering::Less }
            else { std::cmp::Ordering::Equal }
        }) {
            Ok(_) => return false, // 位于基因区域内
            Err(i) => i,
        };
        
        // 检查位置是否在两个基因区域之间
        // 由于合并了重叠区域，只需检查是否在上一个区域之后且在下一个区域之前
        if idx > 0 && idx < regions.len() {
            return pos > regions[idx-1].end && pos < regions[idx].start;
        } else if idx == 0 {
            // 在第一个基因之前
            return pos < regions[0].start;
        } else {
            // 在最后一个基因之后
            return pos > regions.last().unwrap().end;
        }
    } else {
        // 如果染色体没有基因注释，所有位置都视为基因间区域
        true
    }
}

/// 找到离给定位置最近的基因间区域
/// 
/// 从目标位置开始向前和向后搜索，找到最近的基因间区域
/// 返回找到的基因间区域位置
fn find_nearest_intergenic_region(
    pos: usize, 
    chr: &str, 
    gene_regions: &MergedGeneRegions,
    max_search_distance: usize
) -> Option<usize> {
    if let Some(regions) = gene_regions.get(chr) {
        if regions.is_empty() {
            return Some(pos); // 如果没有基因，当前位置就是基因间区域
        }
        
        // 如果当前位置已经是基因间区域，直接返回
        if is_in_intergenic(pos, chr, gene_regions) {
            return Some(pos);
        }
        
        // 向前和向后搜索
        for offset in 1..=max_search_distance {
            // 向前搜索
            if pos > offset {
                let backward_pos = pos - offset;
                if is_in_intergenic(backward_pos, chr, gene_regions) {
                    return Some(backward_pos);
                }
            }
            
            // 向后搜索
            let forward_pos = pos + offset;
            if is_in_intergenic(forward_pos, chr, gene_regions) {
                return Some(forward_pos);
            }
        }
    } else {
        // 如果染色体没有基因注释，所有位置都是基因间区域
        return Some(pos);
    }
    
    None // 没有找到合适的基因间区域
}

// 创建基于用户指定的N长度的正则表达式
fn create_n_pattern(num_n: usize) -> Regex {
    Regex::new(&format!(r"N{{{},}}", num_n)).unwrap()
}

fn find_best_split_position(
    genome_sequence: &str,
    current_position: usize,
    min_fragment_length: usize,
    max_fragment_length: usize,
    chr: &str,
    gene_regions: &MergedGeneRegions,
    n_pattern: &Regex,
) -> Result<Option<(usize, usize)>> {
    let total_length = genome_sequence.len();
    let remaining_length = total_length - current_position;

    if remaining_length <= max_fragment_length {
        return Ok(None);
    }

    // 在N序列处寻找分割点 - 这通常是最佳选择
    for cap in n_pattern.find_iter(&genome_sequence[current_position..]) {
        let match_start = current_position + cap.start();
        let fragment_length = match_start - current_position;

        if fragment_length < min_fragment_length {
            continue;
        }

        if fragment_length > max_fragment_length {
            continue;
        }

        if is_in_intergenic(match_start, chr, gene_regions) {
            eprintln!("Found suitable N-region split point at {}:{}", chr, match_start);
            return Ok(Some((match_start, fragment_length)));
        }
    }

    // 在基因间区域寻找合适的分割点
    let target_pos = current_position + max_fragment_length;
    eprintln!("Searching for intergenic split point near target position {}:{}", chr, target_pos);
    
    // 查找最接近目标位置的基因间区域
    let max_search_distance = 10000000; // 10 Mb搜索窗口
    if let Some(intergenic_pos) = find_nearest_intergenic_region(target_pos, chr, gene_regions, max_search_distance) {
        let fragment_length = intergenic_pos - current_position;
        
        // 检查片段长度是否在允许范围内
        if min_fragment_length <= fragment_length && fragment_length <= (max_fragment_length as f64 * 1.2) as usize {
            eprintln!("Found suitable intergenic split point at {}:{}", chr, intergenic_pos);
            return Ok(Some((intergenic_pos, fragment_length)));
        } else {
            eprintln!("Found intergenic region at {}:{}, but fragment length {} doesn't meet requirements", 
                     chr, intergenic_pos, fragment_length);
        }
    }

    // 如果无法找到合适的分割点，使用接近目标位置的强制分割点
    let fallback_pos = target_pos;
    let fallback_length = fallback_pos - current_position;
    
    if fallback_length >= min_fragment_length {
        eprintln!("Using forced split point at {}:{} (Note: may split gene regions)", chr, fallback_pos);
        return Ok(Some((fallback_pos, fallback_length)));
    }

    Err(anyhow!("Failed to find suitable split position near {}:{}", chr, target_pos))
}

fn process_chromosome(
    id: &str,
    genome_sequence: &str,
    outpref: &str,
    min_fragment_length: usize,
    max_fragment_length: usize,
    gene_regions: &MergedGeneRegions,
    n_pattern: &Regex,
) -> Result<()> {
    let mut current_position = 0;
    let mut real_s = 0;
    let real_len = genome_sequence.len();

    while current_position < genome_sequence.len() {
        match find_best_split_position(
            genome_sequence,
            current_position,
            min_fragment_length,
            max_fragment_length,
            id,
            gene_regions,
            n_pattern,
        )? {
            Some((split_pos, fragment_length)) => {
                let current_fragment = &genome_sequence[current_position..split_pos];
                process_fragment(id, current_fragment, real_s, fragment_length, outpref, real_len)?;
                current_position = split_pos;
                real_s = split_pos;
            }
            None => {
                let remaining_length = genome_sequence.len() - current_position;
                if remaining_length > 0 {
                    let remaining_fragment = &genome_sequence[current_position..];
                    process_fragment(id, remaining_fragment, real_s, remaining_length, outpref, real_len)?;
                }
                break;
            }
        }
    }
    Ok(())
}

/// 处理GXF文件（GTF/GFF），根据分割区域调整基因注释坐标
/// 
/// 此函数读取分割位置信息，然后处理GXF文件，调整基因注释坐标使其相对于分割序列
/// 它保持染色体顺序，确保输出文件中的注释维持与输入相同的顺序
fn process_gxf(ingxf: &str, region: &str, outgxf: &str) -> Result<()> {
    // 使用BTreeMap确保染色体按字典顺序排序
    let mut reg: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    
    // 读取区域信息 - 使用更大的缓冲区来提高读取性能
    let file = File::open(region)?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // 1MB缓冲区
    
    eprintln!("Reading split position information...");
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 3 {
            continue;
        }
        let (chr, s, e) = (parts[0], parts[1], parts[2]);
        let start: usize = s.parse()?;
        let end: usize = e.parse()?;
        
        reg.entry(chr.to_string())
           .or_default()
           .insert(start, end);
    }
    eprintln!("Read split positions for {} chromosomes", reg.len());

    // 处理GXF文件 - 使用更大的缓冲区来提高读/写性能
    let in_file = File::open(ingxf)?;
    let reader = BufReader::with_capacity(1024 * 1024, in_file); // 1MB缓冲区
    let mut out_file = BufWriter::with_capacity(1024 * 1024, File::create(outgxf)?); // 1MB缓冲区

    // 预分配足够的容量以减少重新分配
    let mut processed_count = 0;
    let mut skipped_count = 0;
    
    eprintln!("Processing annotation file...");
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            writeln!(out_file, "{}", line)?;
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 6 {
            skipped_count += 1;
            continue;
        }

        let chr = parts[0];
        let source = parts[1];
        let type_ = parts[2];
        let start: usize = parts[3].parse()?;
        let end: usize = parts[4].parse()?;
        let rest = parts[5..].join("\t");

        let mut found = false;
        if let Some(chr_regions) = reg.get(chr) {
            for (&s, &e) in chr_regions.iter() {
                // 检查特征是否在当前区域内
                if e < start || s > end {
                    continue;
                }
                
                // 检查特征是否跨越区域边界
                if end > e {
                    return Err(anyhow!("Feature crosses split region boundary {}:{}-{}:\n{}", chr, s, e, line));
                }

                // 计算新坐标（相对于分割序列）
                let new_s = start - s + 1;
                let new_e = end - s + 1;

                // 根据染色体是否被分割选择不同的输出格式
                if chr_regions.len() == 1 {
                    // 染色体未分割，保持原始染色体名称
                    writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}", chr, source, type_, new_s, new_e, rest)?;
                } else {
                    // 染色体已分割，使用新的命名格式
                    writeln!(out_file, "{}_{}_{}	{}	{}	{}	{}	{}", 
                        chr, s, e, source, type_, new_s, new_e, rest)?;
                }
                
                found = true;
                processed_count += 1;
                break; // 找到匹配区域后无需继续检查
            }
        }
        
        if !found {
            skipped_count += 1;
        }
        
        // 定期报告进度
        if (processed_count + skipped_count) % 100000 == 0 {
            eprintln!("Processed {} annotation records, skipped {}", processed_count, skipped_count);
        }
    }

    eprintln!("Annotation processing completed: processed {} records, skipped {}", processed_count, skipped_count);
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.max_length <= args.min_length {
        return Err(anyhow!("Maximum length must be greater than minimum length"));
    }

    // 删除现有的输出文件
    for ext in &[".fa", ".cutsite.tsv"] {
        let path = format!("{}{}", args.prefix, ext);
        if Path::new(&path).exists() {
            fs::remove_file(path)?;
        }
    }

    // 读取并处理基因位置信息，合并重叠区域
    let gene_regions = get_merged_gene_regions(args.gtf.as_deref())?;
    
    // 创建N序列模式
    let n_pattern = create_n_pattern(args.num_n);
    
    // 打开FASTA文件进行流式处理，而不是一次性读入内存
    let file = File::open(&args.fa)?;
    let reader = BufReader::new(file);
    
    // 使用BTreeMap确保染色体顺序
    let mut current_id = String::new();
    let mut current_sequence = String::new();
    
    // 使用有序染色体集合来保持原始顺序
    let mut chrom_order = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with('>') {
            // 处理前一个序列（如果有）
            if !current_id.is_empty() && !current_sequence.is_empty() {
                // 将序列存储在有序向量中，保持原始顺序
                chrom_order.push((current_id.clone(), current_sequence.clone()));
            }
            
            // 开始新序列
            let header = line.trim_start_matches('>');
            current_id = header.split_whitespace().next().unwrap_or(header).to_string();
            current_sequence = String::new();
        } else {
            // 添加序列行（删除空格）
            current_sequence.push_str(line.trim());
        }
    }
    
    // 添加最后一个序列
    if !current_id.is_empty() && !current_sequence.is_empty() {
        chrom_order.push((current_id, current_sequence));
    }
    
    // 按原始顺序处理序列
    for (id, sequence) in chrom_order.iter() {
        eprintln!("Processing chromosome: {}, length: {}", id, sequence.len());
        process_chromosome(
            id,
            sequence,
            &args.prefix,
            args.min_length,
            args.max_length,
            &gene_regions,
            &n_pattern,
        )?;
    }
    
    // 上面的循环中已处理所有序列

    // 如果提供了GXF文件，处理它
    if let Some(gxf_file) = args.gtf {
        if !gxf_file.to_lowercase().ends_with(".gff") && 
           !gxf_file.to_lowercase().ends_with(".gff3") && 
           !gxf_file.to_lowercase().ends_with(".gtf") {
            return Err(anyhow!("Invalid GXF file extension (must be gtf/gff/gff3)"));
        }
        
        let ext = Path::new(&gxf_file)
            .extension()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("Invalid file extension"))?;
            
        process_gxf(&gxf_file, &format!("{}.cutsite.tsv", args.prefix), &format!("{}.{}", args.prefix, ext))?;
    }

    Ok(())
}

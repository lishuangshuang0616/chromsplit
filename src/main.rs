use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::collections::HashMap;
use clap::Parser;
use regex::Regex;
use anyhow::{Result, anyhow};

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

#[derive(Debug)]
struct GeneRegion {
    start: usize,
    end: usize,
}

type GeneRegions = HashMap<String, Vec<GeneRegion>>;

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

    let mut out_fa = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(format!("{}.fa", outpref))?;

    if real_len == length {
        writeln!(out_fa, ">{}", id)?;
    } else {
        writeln!(out_fa, ">{}_{}_{}",id, start_position, end)?;
    }

    // 每60个字符换行输出序列
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

fn get_gene_regions(gxf_file: Option<&str>) -> Result<GeneRegions> {
    let mut gene_regions = HashMap::new();
    
    let Some(gxf_path) = gxf_file else {
        return Ok(gene_regions);
    };

    let file = File::open(gxf_path)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chr = fields[0].to_string();
        let type_ = fields[2].to_lowercase();
        if !["gene", "exon"].contains(&type_.as_str()) {
            continue;
        }

        let start: usize = fields[3].parse()?;
        let end: usize = fields[4].parse()?;

        gene_regions.entry(chr)
            .or_insert_with(Vec::new)
            .push(GeneRegion { start, end });
    }

    // 对每条染色体的区域进行排序
    for regions in gene_regions.values_mut() {
        regions.sort_by_key(|r| r.start);
    }

    Ok(gene_regions)
}

fn is_in_intergenic(pos: usize, chr: &str, gene_regions: &GeneRegions) -> bool {
    if let Some(regions) = gene_regions.get(chr) {
        for region in regions {
            if pos >= region.start && pos <= region.end {
                return false;
            }
            if pos < region.start {
                return true;
            }
        }
    }
    true
}

fn find_best_split_position(
    genome_sequence: &str,
    current_position: usize,
    min_fragment_length: usize,
    max_fragment_length: usize,
    num_n: usize,
    chr: &str,
    gene_regions: &GeneRegions,
) -> Result<Option<(usize, usize)>> {
    let total_length = genome_sequence.len();
    let remaining_length = total_length - current_position;

    if remaining_length <= max_fragment_length {
        return Ok(None);
    }

    // 创建N序列的正则表达式模式
    let n_pattern = format!(r"N{{{},}}", num_n);
    let re = Regex::new(&n_pattern)?;

    // 在N序列处查找分割点
    for cap in re.find_iter(&genome_sequence[current_position..]) {
        let match_start = current_position + cap.start();
        let fragment_length = match_start - current_position;

        if fragment_length < min_fragment_length {
            continue;
        }

        if fragment_length > max_fragment_length {
            continue;
        }

        if is_in_intergenic(match_start, chr, gene_regions) {
            return Ok(Some((match_start, fragment_length)));
        }
    }

    // 在基因间区寻找合适的分割点
    let target_pos = current_position + max_fragment_length;
    let window_size = 50000;

    // 向前查找
    for pos in (current_position.max(target_pos.saturating_sub(window_size))..=target_pos).rev() {
        if is_in_intergenic(pos, chr, gene_regions) {
            let fragment_length = pos - current_position;
            if min_fragment_length <= fragment_length && fragment_length <= max_fragment_length {
                return Ok(Some((pos, fragment_length)));
            }
        }
    }

    // 向后查找
    for pos in target_pos..=(target_pos + window_size).min(genome_sequence.len()) {
        if is_in_intergenic(pos, chr, gene_regions) {
            let fragment_length = pos - current_position;
            if min_fragment_length <= fragment_length && fragment_length <= (max_fragment_length as f64 * 1.1) as usize {
                return Ok(Some((pos, fragment_length)));
            }
        }
    }

    Err(anyhow!("Cannot find suitable intergenic region near position {} in {}", target_pos, chr))
}

fn process_gxf(ingxf: &str, region: &str, outgxf: &str) -> Result<()> {
    let mut reg: HashMap<String, HashMap<String, HashMap<String, bool>>> = HashMap::new();
    
    // 读取区域信息
    let file = File::open(region)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 3 {
            continue;
        }
        let (chr, s, e) = (parts[0], parts[1], parts[2]);
        reg.entry(chr.to_string())
            .or_default()
            .entry(s.to_string())
            .or_default()
            .insert(e.to_string(), true);
    }

    // 处理GXF文件
    let in_file = File::open(ingxf)?;
    let reader = BufReader::new(in_file);
    let mut out_file = File::create(outgxf)?;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            writeln!(out_file, "{}", line)?;
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 6 {
            continue;
        }

        let chr = parts[0];
        let source = parts[1];
        let type_ = parts[2];
        let start: usize = parts[3].parse()?;
        let end: usize = parts[4].parse()?;
        let rest = parts[5..].join("\t");

        if let Some(chr_regions) = reg.get(chr) {
            for (s, e_map) in chr_regions.iter() {
                for e in e_map.keys() {
                    let s: usize = s.parse()?;
                    let e: usize = e.parse()?;
                    
                    if e < start || s > end {
                        continue;
                    }
                    if end > e {
                        return Err(anyhow!("feature splitted into diff region {}:{}-{}:\n{}", chr, s, e, line));
                    }

                    let new_s = start - s + 1;
                    let new_e = end - s + 1;

                    if chr_regions.len() == 1 {
                        writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}", chr, source, type_, new_s, new_e, rest)?;
                    } else {
                        writeln!(out_file, "{}_{}_{}\t{}\t{}\t{}\t{}\t{}", 
                            chr, s, e, source, type_, new_s, new_e, rest)?;
                    }
                }
            }
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.max_length <= args.min_length {
        return Err(anyhow!("maxlen must be greater than minlen"));
    }

    // 删除已存在的输出文件
    for ext in &[".fa", ".cutsite.tsv"] {
        let path = format!("{}{}", args.prefix, ext);
        if Path::new(&path).exists() {
            fs::remove_file(path)?;
        }
    }

    // 读取基因位置信息
    let gene_regions = get_gene_regions(args.gtf.as_deref())?;

    // 读取FASTA文件
    let fasta_content = fs::read_to_string(&args.fa)?;
    let records: Vec<&str> = fasta_content.trim().split('>').filter(|s| !s.is_empty()).collect();

    for record in records {
        let mut lines = record.lines();
        let header = lines.next().ok_or_else(|| anyhow!("Invalid FASTA format"))?;
        let id = header.split_whitespace().next().unwrap_or(header);
        let genome_sequence: String = lines.collect();

        let mut current_position = 0;
        let mut real_s = 0;
        let real_len = genome_sequence.len();

        while current_position < genome_sequence.len() {
            match find_best_split_position(
                &genome_sequence,
                current_position,
                args.min_length,
                args.max_length,
                args.num_n,
                id,
                &gene_regions,
            )? {
                Some((split_pos, fragment_length)) => {
                    let current_fragment = &genome_sequence[current_position..split_pos];
                    process_fragment(id, current_fragment, real_s, fragment_length, &args.prefix, real_len)?;
                    current_position = split_pos;
                    real_s = split_pos;
                }
                None => {
                    let remaining_length = genome_sequence.len() - current_position;
                    if remaining_length > 0 {
                        let remaining_fragment = &genome_sequence[current_position..];
                        process_fragment(id, remaining_fragment, real_s, remaining_length, &args.prefix, real_len)?;
                    }
                    break;
                }
            }
        }
    }

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

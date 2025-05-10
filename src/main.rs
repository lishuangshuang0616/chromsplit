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

/// Read gene region information from GXF file (GTF/GFF)
/// 
/// This function parses the GXF file, extracts gene and exon position information, and organizes by chromosome
/// Returns a data structure indexed by chromosome name containing gene regions
fn get_gene_regions(gxf_file: Option<&str>) -> Result<GeneRegions> {
    let mut gene_regions = HashMap::new();
    
    let Some(gxf_path) = gxf_file else {
        return Ok(gene_regions);
    };

    eprintln!("Reading gene annotation file: {}", gxf_path);
    let file = File::open(gxf_path)?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // Using 1MB buffer to improve reading performance

    let mut line_count = 0;
    let mut feature_count = 0;
    let mut skipped_count = 0;
    
    // Pre-allocate a reasonable size vector to reduce reallocation
    let mut chr_capacity: HashMap<String, usize> = HashMap::new();

    // First scan, count the number of features on each chromosome
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
        if !["gene", "exon"].contains(&type_.as_str()) {
            continue;
        }
        
        *chr_capacity.entry(chr).or_insert(0) += 1;
        feature_count += 1;
    }
    
    eprintln!("First scan completed: total {} lines, {} chromosomes, {} features, skipped {} lines", 
              line_count, chr_capacity.len(), feature_count, skipped_count);
    
    // Pre-allocate vector capacity for each chromosome
    for (chr, count) in &chr_capacity {
        gene_regions.insert(chr.clone(), Vec::with_capacity(*count));
    }
    
    // Second scan, actually read features
    let file = File::open(gxf_path)?;
    let reader = BufReader::with_capacity(1024 * 1024, file);
    line_count = 0;
    feature_count = 0;
    
    for line in reader.lines() {
        let line = line?;
        line_count += 1;
        
        if line_count % 1000000 == 0 {
            eprintln!("Processed {} lines", line_count);
        }
        
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
            
        feature_count += 1;
    }

    // Sort regions for each chromosome
    for regions in gene_regions.values_mut() {
        regions.sort_by_key(|r| r.start);
    }
    
    eprintln!("Gene region reading completed: processed {} features, {} chromosomes", feature_count, gene_regions.len());

    Ok(gene_regions)
}

/// Determine if a given position is in an intergenic region
/// 
/// This function uses binary search to quickly locate the position and checks if it is within or near any gene region
/// If the chromosome has no gene annotations, all positions are considered intergenic regions
fn is_in_intergenic(pos: usize, chr: &str, gene_regions: &GeneRegions) -> bool {
    if let Some(regions) = gene_regions.get(chr) {
        // If there are no gene regions, consider it intergenic
        if regions.is_empty() {
            return true;
        }
        
        // 使用二分查找快速定位位置
        match regions.binary_search_by_key(&pos, |r| r.start) {
            Ok(_) => false, // At gene start position
            Err(i) => {
                // Check previous gene - if position is within range of previous gene, it's not intergenic
                if i > 0 && pos <= regions[i-1].end {
                    return false;
                }
                
                // Check current gene (if exists) - if position is within range of current gene, it's not intergenic
                if i < regions.len() && pos >= regions[i].start && pos <= regions[i].end {
                    return false;
                }
                
                // Check if within safe distance - increase safe distance to avoid cutting important regulatory regions near genes
                let safe_distance = 5000; // Increased to 5kb safe distance
                
                // Check distance to the next gene
                if i < regions.len() && pos >= regions[i].start.saturating_sub(safe_distance) {
                    return false;
                }
                
                // Check distance to the previous gene
                if i > 0 && pos <= regions[i-1].end + safe_distance {
                    return false;
                }
                
                // Passed all checks, confirmed as intergenic region
                true
            }
        }
    } else {
        // If chromosome has no gene annotations, consider all positions as intergenic regions
        true
    }
}

// Use dynamically created regular expression based on user-specified N length
fn create_n_pattern(num_n: usize) -> Regex {
    Regex::new(&format!(r"N{{{},}}", num_n)).unwrap()
}

fn find_best_split_position(
    genome_sequence: &str,
    current_position: usize,
    min_fragment_length: usize,
    max_fragment_length: usize,
    chr: &str,
    gene_regions: &GeneRegions,
    n_pattern: &Regex,
) -> Result<Option<(usize, usize)>> {
    let total_length = genome_sequence.len();
    let remaining_length = total_length - current_position;

    if remaining_length <= max_fragment_length {
        return Ok(None);
    }

    // Find split points at N sequences - this is usually the best choice
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
            eprintln!("Found suitable N-stretch split point at {}:{}", chr, match_start);
            return Ok(Some((match_start, fragment_length)));
        }
    }

    // Find suitable split points in intergenic regions
    let target_pos = current_position + max_fragment_length;
    // 增加查找窗口大小以提高找到合适分割点的概率
    let forward_window_size = 10000000; // 增加到10Mb
    let backward_window_size = 5000000; // 增加到5Mb

    eprintln!("Searching for intergenic split point near target position {}:{}", chr, target_pos);
    
    // Search forward - prioritize finding before target position
    let forward_start = current_position.max(target_pos.saturating_sub(forward_window_size));
    eprintln!("Searching forward from {} to {}", forward_start, target_pos);
    
    for pos in (forward_start..=target_pos).rev() {
        let is_intergenic = is_in_intergenic(pos, chr, gene_regions);
        if pos % 1000000 == 0 {
            eprintln!("Checking position {}: intergenic = {}", pos, is_intergenic);
        }
        if is_intergenic {
            let fragment_length = pos - current_position;
            if min_fragment_length <= fragment_length && fragment_length <= max_fragment_length {
                eprintln!("Found suitable forward intergenic split point at {}:{}", chr, pos);
                return Ok(Some((pos, fragment_length)));
            }
        }
    }

    // Search backward
    let backward_end = (target_pos + backward_window_size).min(genome_sequence.len());
    eprintln!("Searching backward from {} to {}", target_pos, backward_end);
    
    for pos in target_pos..=backward_end {
        if pos % 1000000 == 0 {
            eprintln!("Checking position {}", pos);
        }
        if is_in_intergenic(pos, chr, gene_regions) {
            let fragment_length = pos - current_position;
            // Allow slightly longer fragments when searching backward
            if min_fragment_length <= fragment_length && fragment_length <= (max_fragment_length as f64 * 1.2) as usize {
                eprintln!("Found suitable backward intergenic split point at {}:{}", chr, pos);
                return Ok(Some((pos, fragment_length)));
            }
        }
    }

    // 如果找不到合适的基因间区，尝试找到最接近目标位置且距离基因最远的位置
    eprintln!("No suitable intergenic region found, searching for optimal fallback position");
    let mut best_fallback_pos = target_pos;
    let mut max_distance_to_gene = 0;
    let fallback_search_range = 2000000; // 增加到2Mb
    
    // 在目标位置附近搜索最佳切割点
    for offset in 0..fallback_search_range {
        // 向前搜索
        if target_pos > offset {
            let pos = target_pos - offset;
            let distance = find_distance_to_nearest_gene(pos, chr, gene_regions);
            if distance > max_distance_to_gene {
                max_distance_to_gene = distance;
                best_fallback_pos = pos;
                if pos % 100000 == 0 {
                    eprintln!("New best fallback position: {}:{} (distance: {}bp)", chr, pos, distance);
                }
            }
        }
        
        // 向后搜索
        let pos = target_pos + offset;
        if pos < genome_sequence.len() {
            let distance = find_distance_to_nearest_gene(pos, chr, gene_regions);
            if distance > max_distance_to_gene {
                max_distance_to_gene = distance;
                best_fallback_pos = pos;
                if pos % 100000 == 0 {
                    eprintln!("New best fallback position: {}:{} (distance: {}bp)", chr, pos, distance);
                }
            }
        }
    }
    
    let fallback_pos = best_fallback_pos;
    let fallback_length = fallback_pos - current_position;
    
    if fallback_length >= min_fragment_length {
        eprintln!("Using fallback split position at {}:{} (distance to nearest gene: {}bp)", 
                 chr, fallback_pos, max_distance_to_gene);
        return Ok(Some((fallback_pos, fallback_length)));
    }

    Err(anyhow!("Unable to find suitable split position near {}:{}", chr, target_pos))
}

fn process_chromosome(
    id: &str,
    genome_sequence: &str,
    outpref: &str,
    min_fragment_length: usize,
    max_fragment_length: usize,
    gene_regions: &GeneRegions,
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

/// Process GXF file (GTF/GFF), adjust gene annotation coordinates based on split regions
/// 
/// This function reads split position information, then processes the GXF file, adjusting gene annotation coordinates relative to the split sequences
/// It maintains chromosome order, ensuring annotations in the output file maintain the same order as the input
fn process_gxf(ingxf: &str, region: &str, outgxf: &str) -> Result<()> {
    // Use BTreeMap to ensure chromosomes are in dictionary order
    let mut reg: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    
    // Read region information - use larger buffer to improve reading performance
    let file = File::open(region)?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // 1MB buffer
    
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

    // Process GXF file - use larger buffer to improve read/write performance
    let in_file = File::open(ingxf)?;
    let reader = BufReader::with_capacity(1024 * 1024, in_file); // 1MB buffer
    let mut out_file = BufWriter::with_capacity(1024 * 1024, File::create(outgxf)?); // 1MB buffer

    // Pre-allocate sufficient capacity to reduce reallocation
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
                // Check if feature is within current region
                if e < start || s > end {
                    continue;
                }
                
                // Check if feature crosses region boundary
                if end > e {
                    return Err(anyhow!("Feature crosses split region boundary {}:{}-{}:\n{}", chr, s, e, line));
                }

                // Calculate new coordinates (relative to split sequence)
                let new_s = start - s + 1;
                let new_e = end - s + 1;

                // Choose different output format based on whether chromosome is split
                if reg[chr].len() == 1 {
                    // Chromosome not split, keep original chromosome name
                    writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}", chr, source, type_, new_s, new_e, rest)?;
                } else {
                    // Chromosome split, use new naming format
                    writeln!(out_file, "{}_{}_{}	{}	{}	{}	{}	{}", 
                        chr, s, e, source, type_, new_s, new_e, rest)?;
                }
                
                found = true;
                processed_count += 1;
                break; // No need to continue checking after finding matching region
            }
        }
        
        if !found {
            skipped_count += 1;
        }
        
        // Periodically report progress
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
        return Err(anyhow!("maxlen must be greater than minlen"));
    }

    // Delete existing output files
    for ext in &[".fa", ".cutsite.tsv"] {
        let path = format!("{}{}", args.prefix, ext);
        if Path::new(&path).exists() {
            fs::remove_file(path)?;
        }
    }

    // Read gene position information
    let gene_regions = get_gene_regions(args.gtf.as_deref())?;
    
    // Create N sequence pattern
    let n_pattern = create_n_pattern(args.num_n);
    
    // Open FASTA file for streaming processing, rather than reading all at once into memory
    let file = File::open(&args.fa)?;
    let reader = BufReader::new(file);
    
    // Use BTreeMap to ensure chromosome order
    let mut current_id = String::new();
    let mut current_sequence = String::new();
    
    // Use ordered chromosome collection to preserve original order
    let mut chrom_order = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with('>') {
            // Process previous sequence (if any)
            if !current_id.is_empty() && !current_sequence.is_empty() {
                // Store sequence in ordered vector, preserving original order
                chrom_order.push((current_id.clone(), current_sequence.clone()));
            }
            
            // Start new sequence
            let header = line.trim_start_matches('>');
            current_id = header.split_whitespace().next().unwrap_or(header).to_string();
            current_sequence = String::new();
        } else {
            // Add sequence line (remove whitespace)
            current_sequence.push_str(line.trim());
        }
    }
    
    // Add last sequence
    if !current_id.is_empty() && !current_sequence.is_empty() {
        chrom_order.push((current_id, current_sequence));
    }
    
    // Process sequences in original order
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
    
    // All sequences have been processed in the loop above

    // If GXF file is provided, process it
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

/// Calculate the distance from a position to the nearest gene
/// 
/// This function finds the minimum distance from the given position to any gene on the chromosome
/// Returns the distance in base pairs, or 0 if the position is within a gene
/// Calculate the distance from a position to the nearest gene
/// 
/// This function finds the minimum distance from the given position to any gene boundary
/// Returns usize::MAX if there are no genes on the chromosome
/// Returns 0 if the position is within a gene
fn find_distance_to_nearest_gene(pos: usize, chr: &str, gene_regions: &GeneRegions) -> usize {
    if let Some(regions) = gene_regions.get(chr) {
        if regions.is_empty() {
            return usize::MAX; // No genes on this chromosome
        }
        
        // Use binary search to find the nearest genes
        match regions.binary_search_by_key(&pos, |r| r.start) {
            Ok(_) => 0, // Position is exactly at a gene start
            Err(i) => {
                let mut min_distance = usize::MAX;
                
                // Check distance to previous gene (if exists)
                if i > 0 {
                    let prev_gene = &regions[i-1];
                    if pos <= prev_gene.end {
                        return 0; // Position is within the previous gene
                    }
                    min_distance = min_distance.min(pos - prev_gene.end);
                }
                
                // Check distance to next gene (if exists)
                if i < regions.len() {
                    let next_gene = &regions[i];
                    if pos >= next_gene.start && pos <= next_gene.end {
                        return 0; // Position is within the next gene
                    }
                    min_distance = min_distance.min(next_gene.start.saturating_sub(pos));
                }
                
                // Check additional nearby genes to ensure we find the truly closest one
                // This handles cases where gene order might not perfectly correlate with position
                let search_range = 5; // Check 5 genes in each direction
                
                // Check additional previous genes
                for j in 2..=search_range {
                    if i >= j {
                        let prev_gene = &regions[i-j];
                        if pos <= prev_gene.end {
                            return 0; // Position is within this gene
                        }
                        min_distance = min_distance.min(pos - prev_gene.end);
                    } else {
                        break;
                    }
                }
                
                // Check additional next genes
                for j in 1..search_range {
                    if i+j < regions.len() {
                        let next_gene = &regions[i+j];
                        if pos >= next_gene.start && pos <= next_gene.end {
                            return 0; // Position is within this gene
                        }
                        min_distance = min_distance.min(next_gene.start.saturating_sub(pos));
                    } else {
                        break;
                    }
                }
                
                min_distance
            }
        }
    } else {
        usize::MAX // No genes on this chromosome
    }
}

use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Result, anyhow};
use regex::Regex;

use crate::args::Args;
use crate::gene::{MergedGeneRegions, get_merged_gene_regions, is_in_intergenic, find_nearest_intergenic_region};
use crate::utils::{create_n_pattern, process_fragment};

/// Process FASTA file and split it based on criteria
pub fn process_fasta(args: &Args) -> Result<()> {
    // Read and process gene position information, merge overlapping regions
    let gene_regions = get_merged_gene_regions(args.gtf.as_deref())?;
    
    // Create N sequence pattern
    let n_pattern = create_n_pattern(args.num_n);
    
    // Open FASTA file for streaming processing
    let file = File::open(&args.fa)?;
    let reader = BufReader::new(file);
    
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
    
    Ok(())
}

/// Process a single chromosome, dividing it into suitable fragments
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

/// Find the best position to split a sequence
/// 
/// This function tries to find suitable split points in the following order:
/// 1. At N-stretches in intergenic regions
/// 2. In intergenic regions near the target position
/// 3. At a fallback position if no suitable intergenic region is found
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
            eprintln!("Found suitable N-region split point at {}:{}", chr, match_start);
            return Ok(Some((match_start, fragment_length)));
        }
    }

    // Find suitable split points in intergenic regions
    let target_pos = current_position + max_fragment_length;
    eprintln!("Searching for intergenic split point near target position {}:{}", chr, target_pos);
    
    // Find the nearest intergenic region
    let max_search_distance = 10000000; // 10 Mb search window
    if let Some(intergenic_pos) = find_nearest_intergenic_region(target_pos, chr, gene_regions, max_search_distance) {
        let fragment_length = intergenic_pos - current_position;
        
        // Check if fragment length is within acceptable range
        if min_fragment_length <= fragment_length && fragment_length <= (max_fragment_length as f64 * 1.2) as usize {
            eprintln!("Found suitable intergenic split point at {}:{}", chr, intergenic_pos);
            return Ok(Some((intergenic_pos, fragment_length)));
        } else {
            eprintln!("Found intergenic region at {}:{}, but fragment length {} doesn't meet requirements", 
                     chr, intergenic_pos, fragment_length);
        }
    }

    // If no suitable split point is found, use a fallback position near the target
    let fallback_pos = target_pos;
    let fallback_length = fallback_pos - current_position;
    
    if fallback_length >= min_fragment_length {
        eprintln!("Using forced split point at {}:{} (Note: may split gene regions)", chr, fallback_pos);
        return Ok(Some((fallback_pos, fallback_length)));
    }

    Err(anyhow!("Failed to find suitable split position near {}:{}", chr, target_pos))
} 
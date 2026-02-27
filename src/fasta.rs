use anyhow::{anyhow, Result};
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};

use crate::args::Args;
use crate::gene::{
    find_nearest_intergenic_region_leftward, find_nearest_intergenic_region_rightward,
    get_merged_gene_regions, is_in_intergenic, MergedGeneRegions,
};
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

    // Open output files
    let mut out_fa = BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(format!("{}.fa", args.prefix))?,
    );
    let mut out_txt = BufWriter::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(format!("{}.cutsite.tsv", args.prefix))?,
    );

    let mut current_id = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;

        if line.starts_with('>') {
            // Process previous sequence (if any)
            if !current_id.is_empty() && !current_sequence.is_empty() {
                eprintln!(
                    "Processing chromosome: {}, length: {}",
                    current_id,
                    current_sequence.len()
                );
                process_chromosome(
                    &current_id,
                    &current_sequence,
                    &mut out_fa,
                    &mut out_txt,
                    args.min_length,
                    args.max_length,
                    &gene_regions,
                    &n_pattern,
                )?;
            }

            // Start new sequence
            let header = line.trim_start_matches('>');
            current_id = header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string();
            current_sequence.clear();
        } else {
            // Add sequence line (remove whitespace)
            current_sequence.push_str(line.trim());
        }
    }

    // Add last sequence
    if !current_id.is_empty() && !current_sequence.is_empty() {
        eprintln!(
            "Processing chromosome: {}, length: {}",
            current_id,
            current_sequence.len()
        );
        process_chromosome(
            &current_id,
            &current_sequence,
            &mut out_fa,
            &mut out_txt,
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
    out_fa: &mut impl std::io::Write,
    out_txt: &mut impl std::io::Write,
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
                process_fragment(
                    id,
                    current_fragment,
                    real_s,
                    fragment_length,
                    out_fa,
                    out_txt,
                    real_len,
                )?;
                current_position = split_pos;
                real_s = split_pos;
            }
            None => {
                let remaining_length = genome_sequence.len() - current_position;
                if remaining_length > 0 {
                    let remaining_fragment = &genome_sequence[current_position..];
                    process_fragment(
                        id,
                        remaining_fragment,
                        real_s,
                        remaining_length,
                        out_fa,
                        out_txt,
                        real_len,
                    )?;
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
            eprintln!(
                "Found suitable N-region split point at {}:{}",
                chr, match_start
            );
            return Ok(Some((match_start, fragment_length)));
        }
    }

    // Find suitable split points in intergenic regions
    let target_pos = current_position + max_fragment_length;
    eprintln!(
        "Searching for intergenic split point near target position {}:{}",
        chr, target_pos
    );

    // Find the nearest intergenic region (LEFTWARD ONLY)
    let max_search_distance_left = 10000000; // 10 Mb search window
    if let Some(intergenic_pos_left) = find_nearest_intergenic_region_leftward(
        target_pos,
        chr,
        gene_regions,
        max_search_distance_left,
    ) {
        // Since we only search leftward, intergenic_pos is guaranteed to be <= target_pos.
        // Thus fragment_length is unconditionally <= max_fragment_length
        let fragment_length = intergenic_pos_left - current_position;

        // Check if fragment length is within acceptable range (must be >= min_fragment_length)
        if min_fragment_length <= fragment_length && fragment_length <= max_fragment_length {
            eprintln!(
                "Found suitable intergenic split point at {}:{} (Leftward Search)",
                chr, intergenic_pos_left
            );
            return Ok(Some((intergenic_pos_left, fragment_length)));
        } else {
            eprintln!("Leftward intergenic region at {}:{} was too short (length: {}, min req: {}). Falling back to rightward...", 
                     chr, intergenic_pos_left, fragment_length, min_fragment_length);
        }
    }

    // RIGHTWARD FALLBACK (WITH 10% TOLERANCE)
    // If leftward search found a segment that was too short, or found nothing,
    // we search rightward but with a strict max tolerance distance of up to 10% of max_length.
    let max_search_distance_right = (max_fragment_length as f64 * 0.1) as usize;

    if let Some(intergenic_pos_right) = find_nearest_intergenic_region_rightward(
        target_pos,
        chr,
        gene_regions,
        max_search_distance_right,
    ) {
        let fragment_length = intergenic_pos_right - current_position;
        if fragment_length <= (max_fragment_length as f64 * 1.1) as usize {
            eprintln!(
                "Found suitable intergenic split point at {}:{} (Rightward Fallback, using 1.1x tolerance)",
                chr, intergenic_pos_right
            );
            return Ok(Some((intergenic_pos_right, fragment_length)));
        } else {
            eprintln!(
                "Rightward intergenic region at {}:{} exceeded 1.1x tolerance limit (length: {})",
                chr, intergenic_pos_right, fragment_length
            );
        }
    }

    // If no suitable split point is found, use a fallback position near the target
    let fallback_pos = target_pos;
    let fallback_length = fallback_pos - current_position;

    if fallback_length >= min_fragment_length {
        eprintln!(
            "Using forced split point at {}:{} (Note: may split gene regions)",
            chr, fallback_pos
        );
        return Ok(Some((fallback_pos, fallback_length)));
    }

    Err(anyhow!(
        "Failed to find suitable split position near {}:{}",
        chr,
        target_pos
    ))
}

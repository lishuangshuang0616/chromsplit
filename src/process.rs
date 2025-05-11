use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::BTreeMap;
use anyhow::{Result, anyhow};

/// Process GXF file (GTF/GFF), adjust gene annotation coordinates based on split regions
/// 
/// This function reads split position information, then processes the GXF file, adjusting gene annotation coordinates
/// relative to the split sequences. It maintains chromosome order, ensuring annotations maintain the same order.
pub fn process_gxf(ingxf: &str, region: &str, outgxf: &str) -> Result<()> {
    // Use BTreeMap to ensure chromosomes are sorted in dictionary order
    let mut reg: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    
    // Read region information with larger buffer for better performance
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

    // Process GXF file with larger buffers for better performance
    let in_file = File::open(ingxf)?;
    let reader = BufReader::with_capacity(1024 * 1024, in_file); // 1MB buffer
    let mut out_file = BufWriter::with_capacity(1024 * 1024, File::create(outgxf)?); // 1MB buffer

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

                // Choose output format based on whether chromosome is split
                if chr_regions.len() == 1 {
                    // Chromosome not split, keep original name
                    writeln!(out_file, "{}\t{}\t{}\t{}\t{}\t{}", chr, source, type_, new_s, new_e, rest)?;
                } else {
                    // Chromosome split, use new naming format
                    writeln!(out_file, "{}_{}_{}	{}	{}	{}	{}	{}", 
                        chr, s, e, source, type_, new_s, new_e, rest)?;
                }
                
                found = true;
                processed_count += 1;
                break; // No need to check other regions after finding a match
            }
        }
        
        if !found {
            skipped_count += 1;
        }
        
        // Report progress periodically
        if (processed_count + skipped_count) % 100000 == 0 {
            eprintln!("Processed {} annotation records, skipped {}", processed_count, skipped_count);
        }
    }

    eprintln!("Annotation processing completed: processed {} records, skipped {}", processed_count, skipped_count);
    Ok(())
}

/// Process genome using predefined cut sites from a file
/// 
/// This function reads cut site information from a file and splits the genome accordingly
/// It outputs split sequences in FASTA format
pub fn process_with_cut_sites(
    fa_file: &str,
    cut_site_file: &str,
    outpref: &str,
) -> Result<()> {
    // Read cut site information
    let mut cut_sites: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    
    eprintln!("Reading cut site information from {}", cut_site_file);
    let file = File::open(cut_site_file)?;
    let reader = BufReader::with_capacity(1024 * 1024, file);
    
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 3 {
            continue;
        }
        
        let chr = parts[0].to_string();
        let start: usize = parts[1].parse()?;
        let end: usize = parts[2].parse()?;
        
        cut_sites.entry(chr)
            .or_insert_with(Vec::new)
            .push((start, end));
    }
    
    // Read and process FASTA file
    eprintln!("Processing FASTA file with predefined cut sites");
    let file = File::open(fa_file)?;
    let reader = BufReader::new(file);
    
    let mut current_id = String::new();
    let mut current_sequence = String::new();
    
    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with('>') {
            // Process previous sequence if exists
            if !current_id.is_empty() && !current_sequence.is_empty() {
                process_sequence_with_cut_sites(&current_id, &current_sequence, &cut_sites, outpref)?;
            }
            
            // Start new sequence
            let header = line.trim_start_matches('>');
            current_id = header.split_whitespace().next().unwrap_or(header).to_string();
            current_sequence = String::new();
        } else {
            // Add sequence line (remove spaces)
            current_sequence.push_str(line.trim());
        }
    }
    
    // Process the last sequence
    if !current_id.is_empty() && !current_sequence.is_empty() {
        process_sequence_with_cut_sites(&current_id, &current_sequence, &cut_sites, outpref)?;
    }
    
    eprintln!("Completed processing with predefined cut sites");
    Ok(())
}

/// Process a single sequence using predefined cut sites
fn process_sequence_with_cut_sites(
    id: &str,
    sequence: &str,
    cut_sites: &BTreeMap<String, Vec<(usize, usize)>>,
    outpref: &str,
) -> Result<()> {
    if let Some(regions) = cut_sites.get(id) {
        eprintln!("Processing chromosome {} with {} predefined regions", id, regions.len());
        
        for &(start, end) in regions {
            // Adjust for 1-based coordinates in cut site file
            let start_idx = start.saturating_sub(1);
            let end_idx = end.min(sequence.len());
            
            if start_idx >= sequence.len() || end_idx <= start_idx {
                eprintln!("Warning: Invalid region {}:{}-{}, skipping", id, start, end);
                continue;
            }
            
            let fragment = &sequence[start_idx..end_idx];
            let _fragment_length = end_idx - start_idx;
            
            // Output the fragment
            let mut out_fa = BufWriter::new(
                fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(format!("{}.fa", outpref))?
            );
            
            if regions.len() == 1 {
                writeln!(out_fa, ">{}", id)?;
            } else {
                writeln!(out_fa, ">{}_{}_{}",id, start, end)?;
            }
            
            // Output sequence with line breaks every 60 characters
            for chunk in fragment.as_bytes().chunks(60) {
                out_fa.write_all(chunk)?;
                out_fa.write_all(b"\n")?;
            }
            
            // Also write to cut site file for consistency
            let mut out_txt = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(format!("{}.cutsite.tsv", outpref))?;
            
            writeln!(out_txt, "{}\t{}\t{}", id, start, end)?;
        }
    } else {
        eprintln!("No cut sites defined for chromosome {}, skipping", id);
    }
    
    Ok(())
} 
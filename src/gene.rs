use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use anyhow::Result;

#[derive(Debug, Clone)]
pub struct GeneRegion {
    pub start: usize,
    pub end: usize,
}

/// Chromosome gene regions (sorted and merged overlapping regions)
pub type MergedGeneRegions = HashMap<String, Vec<GeneRegion>>;

/// Process GTF/GFF file, extract gene regions and merge overlapping regions
/// 
/// This function parses the GTF/GFF file, extracts gene position information, merges overlapping regions, and organizes by chromosome
/// Returns a data structure indexed by chromosome name containing merged gene regions
pub fn get_merged_gene_regions(gxf_file: Option<&str>) -> Result<MergedGeneRegions> {
    let mut gene_regions: HashMap<String, Vec<GeneRegion>> = HashMap::new();
    
    let Some(gxf_path) = gxf_file else {
        return Ok(gene_regions);
    };

    eprintln!("Reading gene annotation file: {}", gxf_path);
    let file = File::open(gxf_path)?;
    let reader = BufReader::with_capacity(1024 * 1024, file); // Using 1MB buffer to improve reading performance

    let mut raw_regions: HashMap<String, Vec<GeneRegion>> = HashMap::new();
    let mut line_count = 0;
    let mut feature_count = 0;
    let mut skipped_count = 0;
    
    // First scan, collect all gene regions
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
    
    // Sort and merge gene regions for each chromosome
    for (chr, regions) in raw_regions {
        if regions.is_empty() {
            continue;
        }
        
        // Sort regions by start position
        let mut sorted_regions = regions;
        sorted_regions.sort_by_key(|r| r.start);
        
        // Merge overlapping regions
        let mut merged_regions = Vec::new();
        let mut current = sorted_regions[0].clone();
        
        for region in sorted_regions.iter().skip(1) {
            // If the current region overlaps with the previous one, merge them
            if region.start <= current.end + 1 {
                // Extend the current region
                current.end = current.end.max(region.end);
            } else {
                // No overlap, add the current region to the results and start a new one
                merged_regions.push(current);
                current = region.clone();
            }
        }
        // Add the last region
        merged_regions.push(current);
        
        eprintln!("Chromosome {}: merged {} raw regions into {} non-overlapping regions", 
                 chr, sorted_regions.len(), merged_regions.len());
        
        gene_regions.insert(chr, merged_regions);
    }
    
    eprintln!("Gene region processing completed: merged regions for {} chromosomes", gene_regions.len());

    Ok(gene_regions)
}

/// Determine if a given position is in an intergenic region
/// 
/// This function uses binary search to quickly find if the position falls between any gene regions
/// If the chromosome has no gene annotations, all positions are considered intergenic regions
pub fn is_in_intergenic(pos: usize, chr: &str, gene_regions: &MergedGeneRegions) -> bool {
    if let Some(regions) = gene_regions.get(chr) {
        // If there are no gene regions, consider it intergenic
        if regions.is_empty() {
            return true;
        }
        
        // Check if before the first gene
        if pos < regions[0].start {
            return true;
        }
        
        // Check if after the last gene
        if pos > regions.last().unwrap().end {
            return true;
        }
        
        // Binary search to find the interval containing pos
        let idx = match regions.binary_search_by(|r| {
            if pos < r.start { std::cmp::Ordering::Greater }
            else if pos > r.end { std::cmp::Ordering::Less }
            else { std::cmp::Ordering::Equal }
        }) {
            Ok(_) => return false, // Position is within a gene region
            Err(i) => i,
        };
        
        // Check if the position is between two gene regions
        // Since we merged overlapping regions, we only need to check if it's after the previous region and before the next one
        if idx > 0 && idx < regions.len() {
            return pos > regions[idx-1].end && pos < regions[idx].start;
        } else if idx == 0 {
            // Before the first gene
            return pos < regions[0].start;
        } else {
            // After the last gene
            return pos > regions.last().unwrap().end;
        }
    } else {
        // If the chromosome has no gene annotations, all positions are intergenic regions
        true
    }
}

/// Find the nearest intergenic region to a given position
/// 
/// Searches forward and backward from the target position to find the nearest intergenic region
/// Returns the position of the found intergenic region
pub fn find_nearest_intergenic_region(
    pos: usize, 
    chr: &str, 
    gene_regions: &MergedGeneRegions,
    max_search_distance: usize
) -> Option<usize> {
    if let Some(regions) = gene_regions.get(chr) {
        if regions.is_empty() {
            return Some(pos); // If there are no genes, the current position is an intergenic region
        }
        
        // If the current position is already an intergenic region, return it directly
        if is_in_intergenic(pos, chr, gene_regions) {
            return Some(pos);
        }
        
        // Search forward and backward
        for offset in 1..=max_search_distance {
            // Search backward
            if pos > offset {
                let backward_pos = pos - offset;
                if is_in_intergenic(backward_pos, chr, gene_regions) {
                    return Some(backward_pos);
                }
            }
            
            // Search forward
            let forward_pos = pos + offset;
            if is_in_intergenic(forward_pos, chr, gene_regions) {
                return Some(forward_pos);
            }
        }
    } else {
        // If the chromosome has no gene annotations, all positions are intergenic regions
        return Some(pos);
    }
    
    None // No suitable intergenic region found
} 
use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

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

        let mut fields = line.split('\t');
        let chr_str = match fields.next() {
            Some(s) => s,
            None => {
                skipped_count += 1;
                continue;
            }
        };
        let _source = fields.next();
        let type_str = match fields.next() {
            Some(s) => s,
            None => {
                skipped_count += 1;
                continue;
            }
        };

        if !type_str.eq_ignore_ascii_case("gene") {
            continue;
        }

        let start_str = match fields.next() {
            Some(s) => s,
            None => {
                skipped_count += 1;
                continue;
            }
        };
        let end_str = match fields.next() {
            Some(s) => s,
            None => {
                skipped_count += 1;
                continue;
            }
        };

        let start: usize = start_str.parse()?;
        let end: usize = end_str.parse()?;

        let chr = chr_str.to_string();

        raw_regions
            .entry(chr)
            .or_insert_with(Vec::new)
            .push(GeneRegion { start, end });

        feature_count += 1;
    }

    eprintln!(
        "Initial scan completed: processed {} features, {} chromosomes, skipped {} lines",
        feature_count,
        raw_regions.len(),
        skipped_count
    );

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

        eprintln!(
            "Chromosome {}: merged {} raw regions into {} non-overlapping regions",
            chr,
            sorted_regions.len(),
            merged_regions.len()
        );

        gene_regions.insert(chr, merged_regions);
    }

    eprintln!(
        "Gene region processing completed: merged regions for {} chromosomes",
        gene_regions.len()
    );

    Ok(gene_regions)
}

/// Determine if a given position is in an intergenic region
///
/// This function uses binary search to quickly find if the position falls between any gene regions
/// If the chromosome has no gene annotations, all positions are considered intergenic regions
pub fn is_in_intergenic(pos_0based: usize, chr: &str, gene_regions: &MergedGeneRegions) -> bool {
    let pos = pos_0based + 1; // Convert 0-based idx to 1-based GTF coord
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
            if pos < r.start {
                std::cmp::Ordering::Greater
            } else if pos > r.end {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        }) {
            Ok(_) => return false, // Position is within a gene region
            Err(i) => i,
        };

        // Check if the position is between two gene regions
        // Since we merged overlapping regions, we only need to check if it's after the previous region and before the next one
        if idx > 0 && idx < regions.len() {
            return pos > regions[idx - 1].end && pos < regions[idx].start;
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

/// Find the nearest intergenic region to a given position by searching LEFTWARD only.
///
/// Searches backward from the target position to find the nearest intergenic region
/// Returns the position of the found intergenic region. By only searching backward,
/// we guarantee that the resulting fragment length will be <= the target_pos limit.
pub fn find_nearest_intergenic_region_leftward(
    pos_0based: usize,
    chr: &str,
    gene_regions: &MergedGeneRegions,
    max_search_distance: usize,
) -> Option<usize> {
    let pos = pos_0based + 1; // 1-based pos
    if let Some(regions) = gene_regions.get(chr) {
        if regions.is_empty() {
            return Some(pos_0based); // 如果没有基因，当前位置就是基因间区域
        }

        // 如果当前位置已经是基因间区域，找到它所在的区间并返回中点
        if is_in_intergenic(pos_0based, chr, gene_regions) {
            // 找到当前位置所在的基因间区域
            if pos < regions[0].start {
                // 在第一个基因之前
                return Some(pos_0based.min((regions[0].start / 2).saturating_sub(1)));
            // 取当前位置和第一个基因起始位置一半的较小值
            } else if pos > regions.last().unwrap().end {
                // 在最后一个基因之后
                return Some(pos_0based); // 保持原位置
            } else {
                // 在两个基因之间
                for i in 0..regions.len() - 1 {
                    if pos > regions[i].end && pos < regions[i + 1].start {
                        // 计算两个基因之间的中点
                        let mid = regions[i].end + (regions[i + 1].start - regions[i].end) / 2;
                        return Some(mid.saturating_sub(1));
                    }
                }
            }
            return Some(pos_0based); // 默认返回原位置
        }

        // 搜索最近的基因间区域 (仅向左/向前搜索以确保不超长)
        let mut best_pos = None;

        // 向前搜索
        for offset in 1..=max_search_distance {
            if pos > offset {
                let backward_pos = pos - offset;
                if is_in_intergenic(backward_pos.saturating_sub(1), chr, gene_regions) {
                    // 找到向前的基因间区域
                    for i in 0..regions.len() - 1 {
                        if backward_pos > regions[i].end && backward_pos < regions[i + 1].start {
                            // 计算两个基因之间的中点
                            let midpoint =
                                regions[i].end + (regions[i + 1].start - regions[i].end) / 2;
                            best_pos = Some(midpoint);
                            break;
                        }
                    }
                    if best_pos.is_none() {
                        // 如果在第一个基因之前
                        if backward_pos < regions[0].start {
                            best_pos = Some(backward_pos.min(regions[0].start / 2));
                        }
                        // 如果在最后一个基因之后
                        else if backward_pos > regions.last().unwrap().end {
                            best_pos = Some(backward_pos);
                        }
                    }

                    if best_pos.is_some() {
                        break;
                    }
                }
            }
        }

        return best_pos.map(|p| p.saturating_sub(1));
    } else {
        // 如果染色体没有基因注释，所有位置都是基因间区域
        return Some(pos_0based);
    }
}

/// Find the nearest intergenic region to a given position by searching RIGHTWARD only.
///
/// Searches forward from the target position to find the nearest intergenic region.
pub fn find_nearest_intergenic_region_rightward(
    pos_0based: usize,
    chr: &str,
    gene_regions: &MergedGeneRegions,
    max_search_distance: usize,
) -> Option<usize> {
    let pos = pos_0based + 1; // 1-based pos
    if let Some(regions) = gene_regions.get(chr) {
        if regions.is_empty() {
            return Some(pos_0based);
        }

        if is_in_intergenic(pos_0based, chr, gene_regions) {
            // ... (same as leftward logic for already in intergenic) ...
            if pos < regions[0].start {
                return Some(pos_0based.min((regions[0].start / 2).saturating_sub(1)));
            } else if pos > regions.last().unwrap().end {
                return Some(pos_0based);
            } else {
                for i in 0..regions.len() - 1 {
                    if pos > regions[i].end && pos < regions[i + 1].start {
                        let mid = regions[i].end + (regions[i + 1].start - regions[i].end) / 2;
                        return Some(mid.saturating_sub(1));
                    }
                }
            }
            return Some(pos_0based);
        }

        let mut best_pos = None;

        // 向后搜索 (Rightward)
        for offset in 1..=max_search_distance {
            let forward_pos = pos + offset;
            if is_in_intergenic(forward_pos.saturating_sub(1), chr, gene_regions) {
                // 找到向后的基因间区域
                for i in 0..regions.len() - 1 {
                    if forward_pos > regions[i].end && forward_pos < regions[i + 1].start {
                        let midpoint = regions[i].end + (regions[i + 1].start - regions[i].end) / 2;
                        best_pos = Some(midpoint);
                        break;
                    }
                }
                if best_pos.is_none() {
                    if forward_pos < regions[0].start {
                        best_pos = Some(forward_pos.min(regions[0].start / 2));
                    } else if forward_pos > regions.last().unwrap().end {
                        best_pos = Some(forward_pos);
                    }
                }

                if best_pos.is_some() {
                    break;
                }
            }
        }

        return best_pos.map(|p| p.saturating_sub(1));
    } else {
        return Some(pos_0based);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_in_intergenic() {
        let mut regions = MergedGeneRegions::new();
        regions.insert(
            "chr1".to_string(),
            vec![
                GeneRegion { start: 10, end: 20 },
                GeneRegion { start: 30, end: 40 },
            ],
        );

        // Before first gene (0-based idx 8 is 1-based coord 9, 9 < 10 matches)
        assert!(is_in_intergenic(8, "chr1", &regions));
        // Inside first gene (0-based idx 9 is 1-based coord 10)
        assert!(!is_in_intergenic(9, "chr1", &regions));
        assert!(!is_in_intergenic(19, "chr1", &regions));
        // Between genes (0-based idx 20 is 1-based coord 21)
        assert!(is_in_intergenic(20, "chr1", &regions));
        assert!(is_in_intergenic(25, "chr1", &regions));
        // Inside second gene (0-based idx 29 is 1-based coord 30)
        assert!(!is_in_intergenic(29, "chr1", &regions));
        assert!(!is_in_intergenic(39, "chr1", &regions));
        // After second gene (0-based idx 40 is 1-based coord 41)
        assert!(is_in_intergenic(40, "chr1", &regions));
        assert!(is_in_intergenic(100, "chr1", &regions));
    }

    #[test]
    fn test_find_nearest_intergenic_leftward() {
        let mut regions = MergedGeneRegions::new();
        regions.insert(
            "chr1".to_string(),
            vec![
                GeneRegion {
                    start: 100,
                    end: 200,
                },
                GeneRegion {
                    start: 300,
                    end: 400,
                },
            ],
        );

        // Target inside intergenic region (idx = 249 -> pos 250)
        let nearest = find_nearest_intergenic_region_leftward(249, "chr1", &regions, 1000);
        assert_eq!(nearest, Some(249));

        // Target inside gene (idx = 149 -> pos 150)
        // Now it ONLY searches leftward. Nearest intergenic before 150 is the region before 100.
        // Midpoint of intergenic before 100 is 100/2 = 50 (1-based), which is 49 (0-based).
        let nearest2 = find_nearest_intergenic_region_leftward(149, "chr1", &regions, 1000);
        let target_pos = nearest2.unwrap();
        assert!(is_in_intergenic(target_pos, "chr1", &regions));
        assert_eq!(target_pos, 49);
    }

    #[test]
    fn test_find_nearest_intergenic_rightward() {
        let mut regions = MergedGeneRegions::new();
        regions.insert(
            "chr1".to_string(),
            vec![
                GeneRegion {
                    start: 100,
                    end: 200,
                },
                GeneRegion {
                    start: 300,
                    end: 400,
                },
            ],
        );

        // Target inside gene (idx = 149 -> pos 150)
        // Now it ONLY searches rightward. Nearest intergenic after 150 is the region between 200 and 300.
        // Midpoint of intergenic is (200 + 300) / 2 = 250 (1-based), which is 249 (0-based).
        let nearest = find_nearest_intergenic_region_rightward(149, "chr1", &regions, 1000);
        let target_pos = nearest.unwrap();
        assert!(is_in_intergenic(target_pos, "chr1", &regions));
        assert_eq!(target_pos, 249);
    }
}

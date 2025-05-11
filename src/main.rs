mod args;
mod fasta;
mod gene;
mod process;
mod utils;

use std::fs;
use std::path::Path;
use anyhow::{Result, anyhow};
use clap::Parser;

use args::Args;
use process::{process_gxf, process_with_cut_sites};
use fasta::process_fasta;

fn main() -> Result<()> {
    let args = Args::parse();

    if args.max_length <= args.min_length {
        return Err(anyhow!("Maximum length must be greater than minimum length"));
    }

    // Delete existing output files
    for ext in &[".fa", ".cutsite.tsv"] {
        let path = format!("{}{}", args.prefix, ext);
        if Path::new(&path).exists() {
            fs::remove_file(path)?;
        }
    }

    // If cut_site is provided, use it directly instead of computing split points
    if let Some(cut_site_file) = &args.cut_site {
        // Process genome using predefined cut sites
        process_with_cut_sites(&args.fa, cut_site_file, &args.prefix)?;
        
        // If GTF/GFF file is provided, process it
        if let Some(gxf_file) = &args.gtf {
            if !gxf_file.to_lowercase().ends_with(".gff") && 
               !gxf_file.to_lowercase().ends_with(".gff3") && 
               !gxf_file.to_lowercase().ends_with(".gtf") {
                return Err(anyhow!("Invalid file extension"));
            }
            
            let out_gxf = format!("{}.{}", args.prefix, gxf_file.split('.').last().unwrap_or("gff"));
            process_gxf(gxf_file, &format!("{}.cutsite.tsv", args.prefix), &out_gxf)?;
        }
        
        return Ok(());
    }

    // Process the FASTA file and split it based on criteria
    process_fasta(&args)?;

    // If GXF file is provided, process it
    if let Some(gxf_file) = &args.gtf {
        if !gxf_file.to_lowercase().ends_with(".gff") && 
           !gxf_file.to_lowercase().ends_with(".gff3") && 
           !gxf_file.to_lowercase().ends_with(".gtf") {
            return Err(anyhow!("Invalid GXF file extension (must be gtf/gff/gff3)"));
        }
        
        let ext = Path::new(gxf_file)
            .extension()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow!("Invalid file extension"))?;
            
        process_gxf(gxf_file, &format!("{}.cutsite.tsv", args.prefix), &format!("{}.{}", args.prefix, ext))?;
    }

    Ok(())
}

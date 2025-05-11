use std::fs;
use std::io::{BufWriter, Write};
use anyhow::Result;
use regex::Regex;

/// Create a regex pattern for N-stretches based on user-specified length
pub fn create_n_pattern(num_n: usize) -> Regex {
    Regex::new(&format!(r"N{{{},}}", num_n)).unwrap()
}

/// Process a single genome fragment and write it to output files
pub fn process_fragment(
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
use anyhow::Result;
use regex::Regex;
use std::io::Write;

/// Create a regex pattern for N-stretches based on user-specified length
pub fn create_n_pattern(num_n: usize) -> Regex {
    Regex::new(&format!(r"N{{{},}}", num_n)).unwrap()
}

pub fn process_fragment(
    id: &str,
    fragment: &str,
    start_position: usize,
    length: usize,
    out_fa: &mut impl Write,
    out_txt: &mut impl Write,
    real_len: usize,
) -> Result<()> {
    let end = start_position + length - 1;
    let start_position = start_position + 1;
    let end = end + 1;

    if real_len == length {
        writeln!(out_fa, ">{}", id)?;
    } else {
        writeln!(out_fa, ">{}_{}_{}", id, start_position, end)?;
    }

    // Output sequence with line breaks every 60 characters
    for chunk in fragment.as_bytes().chunks(60) {
        out_fa.write_all(chunk)?;
        out_fa.write_all(b"\n")?;
    }

    writeln!(out_txt, "{}\t{}\t{}", id, start_position, end)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_n_pattern() {
        let pattern = create_n_pattern(5);
        assert!(pattern.is_match("NNNNN"));
        assert!(pattern.is_match("ATCNNNNNG"));
        assert!(!pattern.is_match("NNNN"));
    }
}

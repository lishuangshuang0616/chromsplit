# ChromSplit - Genome Sequence Splitting Tool

## Overview
ChromSplit is a tool for splitting large genome sequences into smaller fragments. It can split at long N-base stretches or intergenic regions, avoiding disruption of gene annotations. When provided with GFF/GTF annotation files, the tool ensures splits occur only between genes, maintaining gene annotation integrity.

## Features
- Split at long N-base stretches or intergenic regions
- Preserve gene annotation integrity
- Output split FASTA sequence files
- Output split position information files
- Optional annotation file adjustment

## Installation
1. Ensure Rust toolchain is installed
2. Clone this project
3. Build the project:
```bash
cargo build --release
```

## Usage
```bash
./target/release/chromsplit -f input.fa -o output_prefix [options]
```

### Command Line Arguments
| Argument | Description |
|------|------|
| -f, --fasta | Input genome file in FASTA format |
| -g, --gtf | Optional GTF/GFF annotation file |
| -o, --prefix | Output file prefix |
| --min_length | Minimum fragment length (default: 300000000) |
| --max_length | Maximum fragment length (default: 500000000) |

## Output Files
- `output_prefix.fa`: Split FASTA sequences
- `output_prefix.cutsite.tsv`: Split position information
- `output_prefix.[gff/gtf]`: Adjusted annotation file (if GFF/GTF provided)

## Examples
```bash
# Basic usage
./target/release/chromsplit -f genome.fa -o split_genome

# With annotation file
./target/release/chromsplit -f genome.fa -g annotation.gtf -o split_genome
```
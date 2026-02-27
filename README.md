# ChromSplit 🧬✂️

![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)

**ChromSplit** is a high-performance, industry-grade command-line utility written in Rust. It is designed to intelligently split massive genome sequence files (e.g., human, plant genomes) into more manageable fragments for downstream analysis, visualization, and parallel processing. 

Crucially, **ChromSplit** avoids arbitrary cuts that destroy biological context. It dynamically identifies optimal split points based on user-defined criteria—specifically targeting long continuous stretches of unsequenced gaps (`N`s) or long intergenic regions—ensuring that gene components remain completely intact.

## ✨ Key Features

- **Biological Context Preservation**: Smartly targets N-stretches and uses reference GFF/GTF annotations to guarantee splits happen exactly in intergenic regions.
- **Extreme High Performance**: 
  - **Memory-Efficient Streaming**: Implements a $O(N)$ streaming architecture where $N$ is the length of the longest individual chromosome, not the entire genome, effectively bringing memory footprint down by orders of magnitude (e.g., processing a 3GB human genome within <300MB RAM footprint).
  - **Zero-Allocation Parsing**: Utilizes zero-cost abstractions and string slice comparisons to parse multi-gigabyte GTF/GTF annotations with near-instant speeds and minimal heap allocations.
- **Automated Coordinate Adjustment**: Automatically produces a newly shifted `gff`/`gtf` matching your newly split sequences.
- **Explicit Target Cuts**: Supports providing a custom TSV file of exact cut coordinates via `--cut_site`.

---

## 🚀 Installation

Ensure you have the [Rust toolchain](https://rustup.rs/) installed via `rustup`.

Clone the repository and compile in release mode to enable maximum compiler optimizations (LTO, Single Codegen unit).

```bash
git clone https://github.com/lishuangshuang0616/chromsplit.git
cd chromsplit
cargo build --release
```

The optimized binary will be located at `./target/release/chromsplit`.

---

## 📖 Usage

### Command Line Interface

```bash
chromsplit -f <INPUT.fa> -o <OUTPUT_PREFIX> [OPTIONS]
```

### Core Arguments

| Argument | Short | Description | Default |
| :--- | :---: | :--- | :--- |
| `--fasta` | `-f` | **(Required)** Input genome sequence file in FASTA `.fa` format. | |
| `--prefix` | `-o` | **(Required)** Prefix for all generated output files. | |
| `--gtf` | `-g` | Optional GTF/GFF3 annotation file. Required if you want splits to respect gene boundaries and output a shifted annotation file. | `None` |
| `--min_length` | | Minimum length of output scaffold fragments (bp). | `300,000,000` |
| `--max_length` | | Maximum length of output scaffold fragments (bp). | `500,000,000` |
| `--cut_site` | | Optional TSV file containing explicitly predefined split positions (`chr \t start \t end`). Overrides dynamic sliding window splits. | `None` |
| `--Ns` | | Minimum length of consecutive `N` bases to be treated as a valid biological split separator point. (Hidden flag). | `10` |

---

## 🛠️ Examples

### 1. Basic Length-based Splitting (No Annotations)
Split a raw genome strictly based on lengths and sequence gaps.
```bash
./target/release/chromsplit -f hg38.fa -o hg38_split --min_length 100000000 --max_length 200000000
```

### 2. Annotation-Aware Splitting (Recommended)
Provide a GTF to ensure no single gene structure is caught in the crossfire of a split. The model will search up to 10Mb around the target split point to find a safe intergenic gap.
```bash
./target/release/chromsplit \
    -f hg38.fa \
    -g hg38_annotation.gtf \
    -o hg38_split_annotated \
    --min_length 100000000 \
    --max_length 250000000
```

### 3. Splitting with Pre-calculated Custom Targets
If you've identified specific problematic genomic loci where you want an explicit separation, feed it in.
```bash
./target/release/chromsplit -f hg38.fa --cut_site my_custom_sites.tsv -o custom_split
```

---

## 📁 Output Artifacts

Running the tool produces a set of deterministic files sharing the same prefix:

1. **`{prefix}.fa`**: The resulting split sequences in standard multi-FASTA format. Sequences retain exactly 60 structural characters per line. Headers format: `>{chr}_{start}_{end}`.
2. **`{prefix}.cutsite.tsv`**: A log of the actual 1-based start and end fragments mapped back to the origin chromosome (`chr \t start \t end`).
3. **`{prefix}.gtf` / `.gff`**: *Generated only if `-g` was provided.* The structural annotation file translated exactly to the coordinates of the newly split FASTA file.

---

## 🧠 Algorithmic Strategy & Architecture 

1. **Gene Region Accumulator**: `chromsplit` parses your heavy `.gtf` using highly optimized split iterators, grabbing only the `gene` typed objects. Overlapping genes are naturally merged into interval trees to form solid blocks of "Do Not Cut" regions.
2. **Binary Search Intersection**: When processing the `.fa` file streamingly, the engine identifies dense `N` regions or hits the target `max_length`. It then uses $O(\log n)$ binary search upon the interval tree to verify if the candidate cut-site falls perfectly within safe intergenic space.
3. **Memory Streaming**: The FASTA `.fa` is processed one logical block at a time. Extracted fragments are flushed tightly to output disks over buffered writers `BufWriter`, enabling `chromsplit` to process TB-sized assemblies on personal laptops.

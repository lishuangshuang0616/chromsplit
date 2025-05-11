use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    author, 
    version, 
    about = "Split large genome sequences into smaller fragments at N-stretches or intergenic regions",
    long_about = "A tool for splitting large genome sequences into manageable fragments.\n\
    It identifies suitable split points either at long stretches of N bases or in intergenic regions\n\
    to avoid disrupting gene annotations. When a GFF/GTF file is provided, the tool ensures splits\n\
    occur only between genes, maintaining the integrity of gene annotations.\n\
    \n\
    The tool outputs:\n\
    - Split sequences in FASTA format (.fa)\n\
    - Split positions in TSV format (.cutsite.tsv)\n\
    - Adjusted annotation file if GFF/GTF is provided"
)]
pub struct Args {
    /// Input genome sequence file in FASTA format
    #[arg(short = 'f', long = "fasta")]
    pub fa: String,

    /// Optional GTF/GFF annotation file for the genome
    #[arg(short = 'g', long = "gtf")]
    pub gtf: Option<String>,

    /// Prefix for output files (.fa and .cutsite.tsv will be appended)
    #[arg(short = 'o', long = "prefix")]
    pub prefix: String,

    /// Minimum length of consecutive N bases to be used as a split separator
    #[arg(long = "Ns", default_value = "10", hide = true)]
    pub num_n: usize,

    /// Minimum length of output scaffold fragments (in base pairs)
    #[arg(long = "min_length", default_value = "300000000")]
    pub min_length: usize,

    /// Maximum length of output scaffold fragments (in base pairs)
    #[arg(long = "max_length", default_value = "500000000")]
    pub max_length: usize,

    /// Optional cut site file containing predefined split positions
    #[arg(long = "cut_site")]
    pub cut_site: Option<String>,
} 
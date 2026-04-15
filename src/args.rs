use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Path to the Newick tree file
    #[arg(short, long)]
    pub tree_file: PathBuf,

    /// TKF92 lambda (insertion rate)
    #[arg(short, long, default_value_t = 0.1)]
    pub lambda: f64,

    /// TKF92 mu (deletion rate)
    #[arg(short, long, default_value_t = 0.11)]
    pub mu: f64,

    /// TKF92 r (fragmentation parameter)
    #[arg(short, long, default_value_t = 0.8)]
    pub r: f64,

    /// Max insertion length
    #[arg(long, default_value_t = 50)]
    pub max_insertion_length: usize,

    /// Random seed
    #[arg(short, long)]
    pub seed: Option<u64>,

    /// Root sequence length
    #[arg(long)]
    pub root_length: Option<usize>,

    /// Output base directory
    #[arg(short, long, default_value = "data")]
    pub output_dir: PathBuf,

    /// Open results in AliView after simulation
    #[arg(long)]
    pub ali_view: bool,
}

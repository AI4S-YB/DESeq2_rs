use clap::Parser;
use deseq2_rs::data::Contrast;
use deseq2_rs::io;
use deseq2_rs::pipeline::DESeqDataSet;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "deseq2-rs",
    about = "Rust implementation of DESeq2 differential expression analysis"
)]
struct Cli {
    /// Path to tab-separated count matrix (genes x samples)
    #[arg(long)]
    counts: PathBuf,

    /// Path to tab-separated column data (sample metadata)
    #[arg(long)]
    coldata: PathBuf,

    /// Column name in coldata to use for the design formula
    #[arg(long)]
    design: String,

    /// Reference level for the design factor
    #[arg(long)]
    reference: String,

    /// Output path for results TSV
    #[arg(long, default_value = "results.tsv")]
    output: PathBuf,

    /// Significance threshold for independent filtering
    #[arg(long, default_value = "0.1")]
    alpha: f64,

    /// Number of threads for parallel computation
    #[arg(long)]
    threads: Option<usize>,

    /// Directory to export intermediate results (size factors, dispersions, etc.)
    #[arg(long)]
    intermediate_dir: Option<PathBuf>,
}

fn main() {
    let cli = Cli::parse();

    if let Some(n) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .unwrap();
    }

    eprintln!("Loading data...");
    let mut dds =
        DESeqDataSet::from_csv(&cli.counts, &cli.coldata, &cli.design, &cli.reference)
            .expect("Failed to load data");
    eprintln!(
        "Loaded {} genes x {} samples",
        dds.gene_names.len(),
        dds.sample_names.len()
    );

    eprintln!("Running DESeq2 pipeline...");
    dds.run().expect("Pipeline failed");

    if let Some(ref dir) = cli.intermediate_dir {
        eprintln!("Exporting intermediate results to {:?}", dir);
        dds.export_intermediates(dir).expect("Failed to export");
    }

    eprintln!("Extracting results...");
    let results = dds
        .results(Contrast::LastCoefficient)
        .expect("Failed to extract results");

    eprintln!("Writing results to {:?}", cli.output);
    io::write_results(&cli.output, &results).expect("Failed to write");

    let n_sig = results
        .iter()
        .filter(|r| !r.p_adjusted.is_nan() && r.p_adjusted < 0.01 && r.log2_fold_change.abs() > 1.0)
        .count();
    eprintln!(
        "Done. Significant genes (padj<0.01, |log2FC|>1): {}",
        n_sig
    );
}

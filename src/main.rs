use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

const DEFAULT_KMER_LENGTH: u8 = 31;
const DEFAULT_WINDOW_SIZE: u8 = 15;

#[derive(Parser)]
#[command(author, version, about = "Fast minimizer-based coverage analysis for genomic sequences", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Estimate minimizer coverage and depth from reference and reads
    Cov {
        /// Reference FASTA file
        reference: PathBuf,

        /// Reads FASTA/FASTQ file (supports .gz, .zst, .xz compression). Use '-' for stdin, or omit to read from stdin.
        #[arg(default_value = "-")]
        reads: PathBuf,

        /// K-mer length used for minimizer computation (1-61)
        #[arg(short = 'k', long = "kmer-length", default_value_t = DEFAULT_KMER_LENGTH, value_parser = clap::value_parser!(u8).range(1..=61))]
        kmer_length: u8,

        /// Minimizer window size
        #[arg(short = 'w', long = "window-size", default_value_t = DEFAULT_WINDOW_SIZE)]
        window_size: u8,

        /// Number of execution threads (0 = auto)
        #[arg(short = 't', long = "threads", default_value_t = 8)]
        threads: usize,

        /// Path to output file (- for stdout)
        #[arg(short = 'o', long = "output", default_value = "-")]
        output: String,

        /// Suppress progress reporting
        #[arg(short = 'q', long = "quiet", default_value_t = false)]
        quiet: bool,

        /// Output format: table, csv, or json
        #[arg(short = 'f', long = "format", default_value = "table", value_parser = ["table", "csv", "json"])]
        format: String,
    },
}

fn main() -> Result<()> {
    // Check we have either AVX2 or NEON for SIMD acceleration
    #[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
    {
        eprintln!(
            "Warning: SIMD acceleration is unavailable. For best performance, compile with `cargo build --release -C target-cpu=native`"
        );
    }

    let cli = Cli::parse();

    match &cli.command {
        Commands::Cov {
            reference,
            reads,
            kmer_length,
            window_size,
            threads,
            output,
            quiet,
            format,
        } => {
            // Validate k-mer and window size constraints
            let k = *kmer_length as usize;
            let w = *window_size as usize;

            // Check constraints: k <= 61, k+w <= 96, k+w even (ensures k odd and k+w-1 odd)
            if k > 61 || k + w > 96 || (k + w) % 2 != 0 {
                return Err(anyhow::anyhow!(
                    "Invalid k-w combination: k={}, w={}, k+w={} (constraints: k<=61, k+w<=96, k+w even)",
                    k,
                    w,
                    k + w
                ));
            }

            // Configure thread pool if specified (non-zero)
            if *threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .context("Failed to initialize thread pool")?;
            }

            // Parse outfmt
            let output_format = match format.as_str() {
                "table" => grate::OutputFormat::Table,
                "csv" => grate::OutputFormat::Csv,
                "json" => grate::OutputFormat::Json,
                _ => unreachable!("clap should have validated the format"),
            };

            let config = grate::CoverageConfig {
                reference_path: reference.clone(),
                reads_path: reads.clone(),
                kmer_length: *kmer_length,
                window_size: *window_size,
                threads: *threads,
                output_path: if output == "-" {
                    None
                } else {
                    Some(PathBuf::from(output))
                },
                quiet: *quiet,
                output_format,
            };

            config
                .execute()
                .context("Failed to run coverage analysis")?;
        }
    }

    Ok(())
}

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

const DEFAULT_KMER_LENGTH: u8 = 31;
const DEFAULT_WINDOW_SIZE: u8 = 15;

/// Parse sample string with K/M/G/T suffix into bp count
fn parse_sample(s: &str) -> Result<u64> {
    let s = s.trim().to_uppercase();
    let (num_str, multiplier) = if s.ends_with('T') {
        (&s[..s.len()-1], 1_000_000_000_000u64)
    } else if s.ends_with('G') {
        (&s[..s.len()-1], 1_000_000_000u64)
    } else if s.ends_with('M') {
        (&s[..s.len()-1], 1_000_000u64)
    } else if s.ends_with('K') {
        (&s[..s.len()-1], 1_000u64)
    } else {
        (s.as_str(), 1u64)
    };

    let num: u64 = num_str.parse()
        .with_context(|| format!("Invalid sample value: {}", s))?;

    Ok(num * multiplier)
}

#[derive(Parser)]
#[command(author, version, about = "Streaming containment and abundance estimation using minimizers", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Estimate containment and abundance of target sequence(s) in a read file or stream
    Cov {
        /// Path to fasta file containing target sequence record(s)
        targets: PathBuf,

        /// Path to fastx file containing reads (or - for stdin).
        #[arg(default_value = "-")]
        reads: PathBuf,

        /// Minimizer length (1-61)
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

        /// Output format
        #[arg(short = 'f', long = "format", default_value = "table", value_parser = ["table", "csv", "json"])]
        format: String,

        /// Comma-separated abundance thresholds for containment calculation
        #[arg(
            short = 'a',
            long = "abundance-thresholds",
            value_delimiter = ',',
            default_value = "10"
        )]
        abundance_thresholds: Vec<usize>,

        /// Retain only minimizers exclusive to each target
        #[arg(short = 'd', long = "discriminatory", default_value_t = false)]
        discriminatory: bool,

        /// Terminate read processing after approximately this many bases (e.g. 50M, 10G)
        #[arg(short = 's', long = "sample")]
        sample: Option<String>,
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
            targets,
            reads,
            kmer_length,
            window_size,
            threads,
            output,
            quiet,
            format,
            abundance_thresholds,
            discriminatory,
            sample,
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

            // Parse sample if provided
            let sample_bp = if let Some(s) = sample {
                Some(parse_sample(s)?)
            } else {
                None
            };

            let config = grate::CoverageConfig {
                targets_path: targets.clone(),
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
                abundance_thresholds: abundance_thresholds.clone(),
                discriminatory: *discriminatory,
                sample_bp,
            };

            config
                .execute()
                .context("Failed to run coverage analysis")?;
        }
    }

    Ok(())
}

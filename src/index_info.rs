use anyhow::Result;
use std::path::Path;
use std::time::Instant;

use crate::index::load_minimizer_hashes;

/// Show info about an index
pub fn info<P: AsRef<Path>>(index_path: P) -> Result<()> {
    let start_time = Instant::now();

    // Load index file
    let (minimizers, header) = load_minimizer_hashes(&index_path)?;

    // Show index info
    eprintln!("Index information:");
    eprintln!("  Format version: {}", header.format_version);
    eprintln!("  K-mer length (k): {}", header.kmer_length());
    eprintln!("  Window size (w): {}", header.window_length());
    eprintln!("  Distinct minimizer count: {}", minimizers.len());

    let total_time = start_time.elapsed();
    eprintln!("Retrieved index info in {:.2?}", total_time);

    Ok(())
}

use anyhow::Result;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::index::{load_minimizer_hashes, write_minimizers};

/// Combine minimizer indexes (set union)
pub fn union<P: AsRef<Path>>(inputs: &[P], output: Option<&PathBuf>) -> Result<()> {
    let start_time = Instant::now();
    // Check input files
    if inputs.is_empty() {
        return Err(anyhow::anyhow!(
            "No input files provided for union operation"
        ));
    }
    // Load the first file to get header information
    let (mut all_minimizers, header) = load_minimizer_hashes(&inputs[0])?;
    eprintln!(
        "Performing union of indexes (k={}, w={})",
        header.kmer_length(),
        header.window_length()
    );
    eprintln!(
        "Loaded {} minimizers from first index",
        all_minimizers.len()
    );
    // Process remaining files
    for (_i, path) in inputs.iter().skip(1).enumerate() {
        let (minimizers, file_header) = load_minimizer_hashes(path)?;
        // Verify header compat
        if file_header.kmer_length() != header.kmer_length()
            || file_header.window_length() != header.window_length()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index at {:?} has k={}, w={}, but first index has k={}, w={}",
                path.as_ref(),
                file_header.kmer_length(),
                file_header.window_length(),
                header.kmer_length(),
                header.window_length()
            ));
        }
        // Count minimizers before union
        let before_count = all_minimizers.len();
        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);
        eprintln!(
            " Added {} new minimizers, total: {}",
            all_minimizers.len() - before_count,
            all_minimizers.len()
        );
    }

    // Write output
    write_minimizers(&all_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!(
        "United {} indexes with {} total minimizers in {:.2?}",
        inputs.len(),
        all_minimizers.len(),
        total_time
    );

    Ok(())
}

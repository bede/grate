use anyhow::Result;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::index::{load_minimizer_hashes, write_minimizers};

/// Compute the set difference between two minimizer indexes (A - B)
pub fn diff<P: AsRef<Path>>(first: P, second: P, output: Option<&PathBuf>) -> Result<()> {
    let start_time = Instant::now();

    // Load first file
    let (mut first_minimizers, header) = load_minimizer_hashes(&first)?;
    eprintln!(
        "Loaded {} minimizers from first index",
        first_minimizers.len()
    );

    // Load second file
    let (second_minimizers, second_header) = load_minimizer_hashes(&second)?;
    eprintln!(
        "Loaded {} minimizers from second index",
        second_minimizers.len()
    );

    // Validate compatible headers
    if second_header.kmer_length() != header.kmer_length()
        || second_header.window_length() != header.window_length()
    {
        return Err(anyhow::anyhow!(
            "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={}",
            second_header.kmer_length(),
            second_header.window_length(),
            header.kmer_length(),
            header.window_length()
        ));
    }

    // Count minimizers before diff
    let before_count = first_minimizers.len();

    // Perform diff operation - remove all hashes in second_minimizers from main_minimizers
    for hash in &second_minimizers {
        first_minimizers.remove(hash);
    }

    // Report results
    eprintln!(
        "Removed {} minimizers, {} remaining",
        before_count - first_minimizers.len(),
        first_minimizers.len()
    );

    // Write to output
    write_minimizers(&first_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!("Completed difference operation in {:.2?}", total_time);

    Ok(())
}

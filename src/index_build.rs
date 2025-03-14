use crate::index::write_minimizers;
use crate::index_format::IndexHeader;
use crate::minimizers::compute_minimizer_hashes;
use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Build an index of minimizers from a FASTX file
pub fn build<P: AsRef<Path>>(
    input: P,
    kmer_length: usize,
    window_size: usize,
    output: Option<PathBuf>,
) -> Result<()> {
    let start_time = Instant::now();
    let path = input.as_ref();

    // Create a fastx reader using needletail (handles gzip automatically)
    let mut reader = parse_fastx_file(path).context("Failed to open input file")?;

    // Init FxHashSet with 0.5B capacity
    let mut all_minimizers: FxHashSet<u64> =
        FxHashSet::with_capacity_and_hasher(500_000_000, Default::default());

    eprintln!("Indexing (k={}, w={})…", kmer_length, window_size);
    let mut seq_count = 0;
    let mut total_bp = 0;

    while let Some(record_result) = reader.next() {
        let record = record_result.context("Error reading fastx record")?;
        let seq = record.seq();

        // Compute and insert hashes
        all_minimizers.extend(compute_minimizer_hashes(
            seq.as_ref(),
            kmer_length,
            window_size,
        ));

        // Update counters
        seq_count += 1;
        total_bp += seq.len();

        // Report progress
        let id = std::str::from_utf8(record.id()).unwrap_or("unknown");
        eprintln!(
            "  {} ({}bp), total minimizers: {}",
            id,
            seq.len(),
            all_minimizers.len()
        );
    }

    eprintln!(
        "Indexed {} minimizers from {} sequence(s) ({}bp)",
        all_minimizers.len(),
        seq_count,
        total_bp
    );

    // Create header
    let header = IndexHeader::new(kmer_length, window_size);

    // Write to output path or stdout
    write_minimizers(&all_minimizers, &header, output.as_ref())?;

    let total_time = start_time.elapsed();
    eprintln!("Completed in {:.2?}", total_time);

    Ok(())
}

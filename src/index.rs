use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use needletail::parse_fastx_file;

/// Serializable header for the index file
#[derive(Serialize, Deserialize, Debug)]
pub struct IndexHeader {
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
}

impl IndexHeader {
    pub fn new(kmer_length: usize, window_size: usize) -> Self {
        IndexHeader {
            format_version: 2,
            kmer_length: kmer_length as u8,
            window_size: window_size as u8,
        }
    }

    /// Validate header
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != 2 {
            return Err(anyhow::anyhow!(
                "Unsupported index format version: {}",
                self.format_version
            ));
        }

        Ok(())
    }

    /// Get k
    pub fn kmer_length(&self) -> usize {
        self.kmer_length as usize
    }

    /// Get w
    pub fn window_size(&self) -> usize {
        self.window_size as usize
    }
}

/// Load the hashes without spiking memory usage with an extra vec
pub fn load_minimizer_hashes<P: AsRef<Path>>(path: &P) -> Result<(FxHashSet<u64>, IndexHeader)> {
    let file =
        File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    // Deserialize header
    let header: IndexHeader = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialize index header")?;
    header.validate()?;

    // Deserialize the count of minimizers so we can init a FxHashSet with the right capacity
    let count: usize = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialize minimizer count")?;

    // Pre-allocate FxHashSet with correct capacity
    let mut minimizers = FxHashSet::with_capacity_and_hasher(count, Default::default());

    // Populate FxHashSet
    for _ in 0..count {
        let hash: u64 = decode_from_std_read(&mut reader, bincode::config::standard())
            .context("Failed to deserialize minimizer hash")?;
        minimizers.insert(hash);
    }

    Ok((minimizers, header))
}

/// Helper function to write minimizers to output file or stdout
pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    // Create writer based on output path
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        if path.to_string_lossy() == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(
                File::create(path).context("Failed to create output file")?,
            ))
        }
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // Serialize header and minimizers
    let mut writer = BufWriter::new(writer);
    encode_into_std_write(header, &mut writer, bincode::config::standard())
        .context("Failed to serialize index header")?;

    // Serialize the count of minimizers first
    let count = minimizers.len();
    encode_into_std_write(&count, &mut writer, bincode::config::standard())
        .context("Failed to serialize minimizer count")?;

    // Serialize each minimizer directly
    for &hash in minimizers {
        encode_into_std_write(&hash, &mut writer, bincode::config::standard())
            .context("Failed to serialize minimizer hash")?;
    }
    Ok(())
}

/// Build an index of minimizers from a fastx file
pub fn build<P: AsRef<Path>>(
    input: P,
    kmer_length: usize,
    window_size: usize,
    output: Option<PathBuf>,
    capacity_millions: usize,
    threads: usize,
) -> Result<()> {
    let start_time = Instant::now();
    let path = input.as_ref();

    let version: String = env!("CARGO_PKG_VERSION").to_string();
    eprintln!("Deacon v{}", version,);

    // Configure thread pool if specified (non-zero)
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    // Create a fastx reader using needletail (handles gzip automatically)
    let mut reader = parse_fastx_file(path).context("Failed to open input file")?;

    // Init FxHashSet with user-specified capacity
    let capacity = capacity_millions * 1_000_000;
    let mut all_minimizers: FxHashSet<u64> =
        FxHashSet::with_capacity_and_hasher(capacity, Default::default());

    eprintln!("Indexing (k={}, w={})…", kmer_length, window_size);

    let mut seq_count = 0;
    let mut total_bp = 0;

    // Process sequences in batches for parallelization
    let batch_size = 10000;

    loop {
        // Collect a batch of sequences
        let mut batch = Vec::with_capacity(batch_size);
        let mut reached_end = false;

        for _ in 0..batch_size {
            if let Some(record_result) = reader.next() {
                match record_result {
                    Ok(record) => {
                        let seq_data = record.seq().to_vec();
                        let id = record.id().to_vec();
                        batch.push((seq_data, id));
                    }
                    Err(e) => return Err(e.into()),
                }
            } else {
                reached_end = true;
                break;
            }
        }

        if batch.is_empty() {
            break;
        }

        // Process batch in parallel
        let batch_results: Vec<_> = batch
            .par_iter()
            .map(|(seq_data, _id)| {
                // Compute minimizer hashes for this sequence
                crate::minimizers::compute_minimizer_hashes(seq_data, kmer_length, window_size)
            })
            .collect();

        // Merge results sequentially (to avoid concurrent modification of hash set)
        for (i, hashes) in batch_results.iter().enumerate() {
            let (seq_data, id) = &batch[i];

            all_minimizers.extend(hashes.iter());

            seq_count += 1;
            total_bp += seq_data.len();

            let id_str = std::str::from_utf8(id).unwrap_or("unknown");
            eprintln!(
                "  {} ({}bp), total minimizers: {}",
                id_str,
                seq_data.len(),
                all_minimizers.len()
            );
        }

        // Check if we've reached the end of the file
        if reached_end {
            break;
        }
    }

    eprintln!(
        "Indexed {} minimizers from {} sequence(s) ({}bp)",
        all_minimizers.len(),
        seq_count,
        total_bp
    );

    let header = IndexHeader::new(kmer_length, window_size);

    // Write to output path or stdout
    write_minimizers(&all_minimizers, &header, output.as_ref())?;

    let total_time = start_time.elapsed();
    eprintln!("Completed in {:.2?}", total_time);

    Ok(())
}

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
        || second_header.window_size() != header.window_size()
    {
        return Err(anyhow::anyhow!(
            "Incompatible headers: second index has k={}, w={}, but first index has k={}, w={}",
            second_header.kmer_length(),
            second_header.window_size(),
            header.kmer_length(),
            header.window_size()
        ));
    }

    // Count minimizers before diff
    let before_count = first_minimizers.len();

    // Remove all hashes in second_minimizers from main_minimizers
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

/// Show info about an index
pub fn info<P: AsRef<Path>>(index_path: P) -> Result<()> {
    let start_time = Instant::now();

    // Load index file
    let (minimizers, header) = load_minimizer_hashes(&index_path)?;

    // Show index info
    eprintln!("Index information:");
    eprintln!("  Format version: {}", header.format_version);
    eprintln!("  K-mer length (k): {}", header.kmer_length());
    eprintln!("  Window size (w): {}", header.window_size());
    eprintln!("  Distinct minimizer count: {}", minimizers.len());

    let total_time = start_time.elapsed();
    eprintln!("Retrieved index info in {:.2?}", total_time);

    Ok(())
}

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
        header.window_size()
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
            || file_header.window_size() != header.window_size()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index at {:?} has k={}, w={}, but first index has k={}, w={}",
                path.as_ref(),
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_creation() {
        let header = IndexHeader::new(31, 21);

        assert_eq!(header.format_version, 2);
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
    }

    #[test]
    fn test_header_validation() {
        // Valid header
        let valid_header = IndexHeader {
            format_version: 2,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(valid_header.validate().is_ok());

        // Invalid format version
        let invalid_header = IndexHeader {
            format_version: 1, // Unsupported version
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header.validate().is_err());
    }
}

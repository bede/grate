use crate::filter::{Buffers, ProcessingStats};
use crate::minimizers::KmerHasher;
use crate::IndexConfig;
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use paraseq::prelude::{ParallelProcessor, ParallelReader};
use paraseq::Record;
use parking_lot::Mutex;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};
use std::time::Instant;

/// Serialisable header for the index file
#[derive(Serialize, Deserialize, Debug)]
pub struct IndexHeader {
    pub format_version: u8,
    pub kmer_length: u8,
    pub window_size: u8,
}

impl IndexHeader {
    pub fn new(kmer_length: u8, window_size: u8) -> Self {
        IndexHeader {
            format_version: 3,
            kmer_length,
            window_size,
        }
    }

    /// Validate header
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != 2 && self.format_version != 3 {
            return Err(anyhow::anyhow!(
                "Unsupported index format version: {}",
                self.format_version
            ));
        }

        Ok(())
    }

    /// Get k
    pub fn kmer_length(&self) -> u8 {
        self.kmer_length
    }

    /// Get w
    pub fn window_size(&self) -> u8 {
        self.window_size
    }
}

/// Load just the header and count from an index file
pub fn load_header_and_count<P: AsRef<Path>>(path: &P) -> Result<(IndexHeader, usize)> {
    let file =
        File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    // Deserialise header
    let header: IndexHeader = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialise index header")?;
    header.validate()?;

    // Deserialise the count of minimizers
    let count: usize = decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialise minimizer count")?;

    Ok((header, count))
}

static INDEX: OnceLock<(PathBuf, FxHashSet<u64>, IndexHeader)> = OnceLock::new();

pub fn load_minimizer_hashes_cached(
    path: &Path,
) -> Result<(&'static FxHashSet<u64>, &'static IndexHeader)> {
    let (p, minimizers, header) = INDEX.get_or_init(|| {
        let (m, h) = load_minimizer_hashes(path).unwrap();
        (path.to_owned(), m, h)
    });
    assert_eq!(
        p, path,
        "Currently, the server can only have one index loaded."
    );

    Ok((minimizers, header))
}

/// Load the hashes without spiking memory usage with an extra vec
/// This older version uses variable-width integer encoding.
/// Use `deacon index convert` for the new format.
fn load_minimizer_hashes_varint(mut reader: impl std::io::Read) -> Result<FxHashSet<u64>> {
    let config = bincode::config::standard();

    // Deserialise the count of minimizers so we can init a FxHashSet with the right capacity
    let count: usize = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;

    // Populate FxHashSet
    let minimizers: FxHashSet<u64> = (0..count)
        .map(|_| decode_from_std_read(&mut reader, config).unwrap())
        .collect();

    Ok(minimizers)
}

fn load_minimizer_hashes_fixedint(mut reader: impl std::io::Read) -> Result<FxHashSet<u64>> {
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise the count of minimizers so we can init a FxHashSet with the right capacity
    let count: usize = decode_from_std_read(&mut reader, config)
        .context("Failed to deserialise minimizer count")?;

    // Populate FxHashSet
    let mut minimizers = FxHashSet::<u64>::with_capacity_and_hasher(count, Default::default());
    const B: usize = 16 * 1024;
    let mut hashes = vec![0u8; 8 * B];
    for i in (0..count).step_by(B) {
        let batch_count = B.min(count - i);
        let batch = &mut hashes[..8 * batch_count];
        reader.read_exact(batch).unwrap();
        for h in batch.as_chunks::<8>().0 {
            // Read as little-endian.
            minimizers.insert(u64::from_le_bytes(*h));
        }
    }

    Ok(minimizers)
}

/// Load the hashes without spiking memory usage with an extra vec
/// This new version uses fixed-width integer encoding.
/// Use `deacon index convert` to convert from the old format.
pub fn load_minimizer_hashes(path: &Path) -> Result<(FxHashSet<u64>, IndexHeader)> {
    let file = File::open(path).context(format!("Failed to open index file {:?}", path))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let config = bincode::config::standard().with_fixed_int_encoding();

    // Deserialise header
    let header: IndexHeader =
        decode_from_std_read(&mut reader, config).context("Failed to deserialise index header")?;
    header.validate()?;

    let minimizers = match header.format_version {
        2 => load_minimizer_hashes_varint(reader)?,
        3 => load_minimizer_hashes_fixedint(reader)?,
        _ => unreachable!(),
    };

    Ok((minimizers, header))
}

/// Helper function to write minimizers to output file or stdout
pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &IndexHeader,
    output_path: Option<&Path>,
) -> Result<()> {
    // Create writer based on output path
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        if path.as_os_str() == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(
                File::create(path).context("Failed to create output file")?,
            ))
        }
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // Serialise header and minimizers
    let mut writer = BufWriter::new(writer);
    // Use fixed-width little-endian encoding for all values.
    let config = bincode::config::standard().with_fixed_int_encoding();
    encode_into_std_write(header, &mut writer, config)
        .context("Failed to serialise index header")?;

    // Serialise the count of minimizers first
    let count = minimizers.len();
    encode_into_std_write(count, &mut writer, config)
        .context("Failed to serialise minimizer count")?;

    // Serialise each minimizer directly
    for &hash in minimizers {
        encode_into_std_write(hash, &mut writer, config)
            .context("Failed to serialise minimizer hash")?;
    }
    Ok(())
}

fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<paraseq::fastx::Reader<Box<dyn Read + Send>>> {
    let mut reader = paraseq::fastx::Reader::from_optional_path(in_path).unwrap();
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

#[derive(Clone)]
struct BuildIndexProcessor<'c> {
    config: &'c IndexConfig,
    hasher: KmerHasher,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_hashes: FxHashSet<u64>,
    local_record_info: Vec<(Vec<u8>, usize)>, // (record_id, seq_len) for progress output
    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_hashes: Arc<Mutex<FxHashSet<u64>>>,
}

impl<Rf: Record> ParallelProcessor<Rf> for BuildIndexProcessor<'_> {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        crate::minimizers::fill_minimizer_hashes(
            &seq,
            &self.hasher,
            self.config.kmer_length,
            self.config.window_size,
            self.config.entropy_threshold,
            &mut self.buffers,
        );
        self.local_hashes.extend(&self.buffers.hashes);

        // Store record info for progress output (printed sequentially in on_batch_complete)
        if !self.config.quiet {
            self.local_record_info
                .push((record.id().to_vec(), seq.len()));
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local hashes into global hash set and get new total count
        let minimizer_count = {
            let mut global_hashes = self.global_hashes.lock();
            global_hashes.extend(&self.local_hashes);
            self.local_hashes.clear();
            global_hashes.len()
        };

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;
            self.local_stats = ProcessingStats::default();
        }

        // Print progress for each record in this batch (sequentially, after merging)
        if !self.config.quiet {
            for (record_id, seq_len) in &self.local_record_info {
                let id_str = std::str::from_utf8(record_id).unwrap_or("unknown");
                eprintln!(
                    "  {} ({}bp), total minimizers: {}",
                    id_str, seq_len, minimizer_count
                );
            }
            self.local_record_info.clear();
        }

        Ok(())
    }
}

/// Build an index of minimizers from a fastx file
pub fn build(config: &IndexConfig) -> Result<()> {
    let start_time = Instant::now();
    let path = &config.input_path;

    let version: String = env!("CARGO_PKG_VERSION").to_string();

    // Build options string similar to filter
    let mut options = Vec::<String>::new();
    if config.threads > 0 {
        options.push(format!("threads={}", config.threads));
    }

    eprintln!(
        "Deacon v{}; mode: build; input: single; options: {}",
        version,
        options.join(", ")
    );

    // Ensure l = k + w - 1 is odd so that canonicalisation tie breaks work correctly
    let l = config.kmer_length as usize + config.window_size as usize - 1;
    if l % 2 == 0 {
        return Err(anyhow::anyhow!(
            "Constraint violated: k + w - 1 must be odd (k={}, w={})",
            config.kmer_length,
            config.window_size
        ));
    }

    let in_path = if path.as_os_str() == "-" {
        None
    } else {
        Some(path.as_path())
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    eprintln!(
        "Building index (k={}, w={})",
        config.kmer_length, config.window_size
    );

    let mut processor = BuildIndexProcessor {
        config,
        hasher: KmerHasher::new(config.kmer_length as usize),
        local_stats: ProcessingStats::default(),
        buffers: Buffers::default(),
        local_hashes: Default::default(),
        local_record_info: Vec::new(),
        global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
        global_hashes: Arc::new(Mutex::new(FxHashSet::<u64>::default())),
    };
    reader.process_parallel(&mut processor, config.threads)?;

    let all_minimizers = Arc::try_unwrap(processor.global_hashes)
        .unwrap()
        .into_inner();
    let stats = Arc::try_unwrap(processor.global_stats)
        .unwrap()
        .into_inner();

    eprintln!(
        "Indexed {} minimizers from {} record(s) ({}bp)",
        all_minimizers.len(),
        stats.total_seqs,
        stats.total_bp
    );

    let header = IndexHeader::new(config.kmer_length, config.window_size);

    // Write to output path or stdout
    write_minimizers(&all_minimizers, &header, config.output_path.as_deref())?;

    let total_time = start_time.elapsed();
    eprintln!("Completed in {:.2?}", total_time);

    Ok(())
}

#[derive(Clone)]
struct DiffIndexProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_hashes: FxHashSet<u64>,
    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    initial_size: usize,
    global_hashes: Arc<Mutex<FxHashSet<u64>>>,
}

impl<Rf: Record> ParallelProcessor<Rf> for DiffIndexProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        crate::minimizers::fill_minimizer_hashes(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            0.0,
            &mut self.buffers,
        );
        self.local_hashes.extend(&self.buffers.hashes);

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        let len;
        {
            let mut global_hashes = self.global_hashes.lock();
            for h in &self.local_hashes {
                global_hashes.remove(h);
            }
            len = global_hashes.len();
            self.local_hashes.clear();
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                eprintln!(
                    "  Processed {} sequences ({}bp), removed {} minimizers",
                    stats.total_seqs,
                    stats.total_bp,
                    self.initial_size - len
                );
                stats.last_reported = current_gb;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Stream minimizers from a FASTX file or stdin and remove those present in first_minimizers
fn stream_diff_fastx(
    fastx_path: &Path,
    kmer_length: u8,
    window_size: u8,
    first_header: &IndexHeader,
    threads: usize,
    first_minimizers: &mut FxHashSet<u64>,
) -> Result<(usize, usize)> {
    let path = fastx_path;

    // Validate parameters match the first index
    if kmer_length != first_header.kmer_length() || window_size != first_header.window_size() {
        return Err(anyhow::anyhow!(
            "FASTX parameters (k={}, w={}) must match first index (k={}, w={})",
            kmer_length,
            window_size,
            first_header.kmer_length(),
            first_header.window_size()
        ));
    }

    if path.as_os_str() == "-" {
        eprintln!(
            "Second index: processing FASTX from stdin (k={}, w={})…",
            kmer_length, window_size
        );
    } else {
        eprintln!(
            "Second index: processing FASTX from file (k={}, w={})…",
            kmer_length, window_size
        );
    }

    let in_path = if path.as_os_str() == "-" {
        None
    } else {
        Some(path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;

    let mut processor = DiffIndexProcessor {
        kmer_length,
        window_size,
        hasher: KmerHasher::new(kmer_length as usize),
        local_stats: ProcessingStats::default(),
        buffers: Buffers::default(),
        local_hashes: Default::default(),
        global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
        initial_size: first_minimizers.len(),
        global_hashes: Arc::new(Mutex::new(FxHashSet::<u64>::default())),
    };
    reader.process_parallel(&mut processor, threads)?;

    let stats = Arc::try_unwrap(processor.global_stats)
        .unwrap()
        .into_inner();

    eprintln!(
        "Processed {} sequences ({}bp) from FASTX file",
        stats.total_seqs, stats.total_bp
    );

    Ok((stats.total_seqs as usize, stats.total_bp as usize))
}

/// Compute the set difference between two minimizer indexes (A - B)
pub fn diff(
    first: &Path,
    second: &Path,
    kmer_length: Option<u8>,
    window_size: Option<u8>,
    threads: usize,
    output: Option<&Path>,
) -> Result<()> {
    let start_time = Instant::now();

    // Load first file (always an index)
    let (mut first_minimizers, header) = load_minimizer_hashes(first)?;
    eprintln!("First index: loaded {} minimizers", first_minimizers.len());

    // Guess if second file is an index or FASTX file
    let second_minimizers = if let (Some(k), Some(w)) = (kmer_length, window_size) {
        // Second file is a FASTX file - stream diff with provided k, w
        let before_count = first_minimizers.len();
        let (_seq_count, _total_bp) =
            stream_diff_fastx(second, k, w, &header, threads, &mut first_minimizers)?;

        // Report results
        eprintln!(
            "Removed {} minimizers, {} remaining",
            before_count - first_minimizers.len(),
            first_minimizers.len()
        );

        write_minimizers(&first_minimizers, &header, output)?;

        let total_time = start_time.elapsed();
        eprintln!("Completed difference operation in {:.2?}", total_time);

        return Ok(());
    } else {
        // Try to load as index file first
        if let Ok((second_minimizers, second_header)) = load_minimizer_hashes(&second) {
            // Second file is an index file
            eprintln!(
                "Second index: loaded {} minimizers",
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

            second_minimizers
        } else {
            // Second file is not a valid index, treat as FASTX file
            // Use k and w from first index header and do a streaming diff
            let k = header.kmer_length();
            let w = header.window_size();

            // Count minimizers before diff
            let before_count = first_minimizers.len();

            let (_seq_count, _total_bp) =
                stream_diff_fastx(&second, k, w, &header, threads, &mut first_minimizers)?;

            // Report results
            eprintln!(
                "Removed {} minimizers, {} remaining",
                before_count - first_minimizers.len(),
                first_minimizers.len()
            );

            write_minimizers(&first_minimizers, &header, output)?;

            let total_time = start_time.elapsed();
            eprintln!("Completed difference operation in {:.2?}", total_time);

            return Ok(());
        }
    };

    // Handle straightforward index-to-index diffing
    // Count minimizers before diff
    let before_count = first_minimizers.len();

    // Remove all hashes in second_minimizers from first_minimizers
    for hash in &second_minimizers {
        first_minimizers.remove(hash);
    }

    // Report results
    eprintln!(
        "Removed {} minimizers, {} remaining",
        before_count - first_minimizers.len(),
        first_minimizers.len()
    );

    write_minimizers(&first_minimizers, &header, output)?;

    let total_time = start_time.elapsed();
    eprintln!("Completed diff operation in {:.2?}", total_time);

    Ok(())
}

/// Show info about an index
pub fn info(index_path: &Path) -> Result<()> {
    let start_time = Instant::now();

    // Load index file
    let (minimizers, header) = load_minimizer_hashes(index_path)?;

    // Show index info
    eprintln!("Index information:");
    eprintln!("  Format version: {}", header.format_version);
    eprintln!("  K-mer length (k): {}", header.kmer_length());
    eprintln!("  Window size (w): {}", header.window_size());
    eprintln!("  Distinct minimizer count: {}", minimizers.len());

    let total_time = start_time.elapsed();
    eprintln!("Loaded index info in {:.2?}", total_time);

    Ok(())
}

pub fn convert_index(from: &Path, to: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    let (minimizers, mut header) = load_minimizer_hashes(from)?;
    let load_time = start_time.elapsed();
    eprintln!(
        "Loaded index (k={}, w={}) in {:.2?}",
        header.kmer_length, header.window_size, load_time
    );
    assert_eq!(header.format_version, 2);
    header.format_version = 3;
    let start_time = Instant::now();
    write_minimizers(&minimizers, &header, to)?;
    let write_time = start_time.elapsed();
    eprintln!("Converted index in {:.2?}", write_time);

    Ok(())
}

/// Combine minimizer indexes (set union)
pub fn union(inputs: &[PathBuf], output: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    // Check input files
    if inputs.is_empty() {
        return Err(anyhow::anyhow!(
            "No input files provided for union operation"
        ));
    }

    // Read all headers first to determine total capacity needed
    let mut headers_and_counts = Vec::new();

    for path in inputs {
        let (header, count) = load_header_and_count(path)?;
        headers_and_counts.push((header, count));
    }

    // Get header from first file for output
    let header = &headers_and_counts[0].0;

    eprintln!(
        "Performing union of indexes (k={}, w={})",
        header.kmer_length(),
        header.window_size()
    );

    // Verify all headers are compatible
    for (i, (file_header, _)) in headers_and_counts.iter().enumerate() {
        if file_header.kmer_length() != header.kmer_length()
            || file_header.window_size() != header.window_size()
        {
            return Err(anyhow::anyhow!(
                "Incompatible headers: index {} has k={}, w={}, but first index has k={}, w={}",
                i,
                file_header.kmer_length(),
                file_header.window_size(),
                header.kmer_length(),
                header.window_size()
            ));
        }
    }

    // Pre-allocate hash set with total capacity to avoid resizing
    let mut all_minimizers = FxHashSet::<u64>::default();

    // Now load and merge all indexes
    for (i, path) in inputs.iter().enumerate() {
        let (minimizers, _) = load_minimizer_hashes(path)?;
        let before_count = all_minimizers.len();

        // Merge minimizers (set union)
        all_minimizers.extend(minimizers);

        let expected_count = headers_and_counts[i].1;
        eprintln!(
            "Index {}: expected {} minimizers, added {} new, total: {}",
            i + 1,
            expected_count,
            all_minimizers.len() - before_count,
            all_minimizers.len()
        );
    }

    write_minimizers(&all_minimizers, header, output)?;

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

        assert_eq!(header.format_version, 3);
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
    }

    #[test]
    fn test_header_validation() {
        // Valid header
        let valid_header = IndexHeader {
            format_version: 3,
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

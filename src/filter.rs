use crate::index::load_minimizer_hashes;
use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use liblzma::write::XzEncoder;
use packed_seq::SeqVec;
use paraseq::fastx::Reader;
use paraseq::parallel::{
    InterleavedParallelProcessor, InterleavedParallelReader, PairedParallelProcessor,
    PairedParallelReader, ParallelProcessor, ParallelReader,
};
use paraseq::Record;
use parking_lot::Mutex;
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use simd_minimizers;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use xxhash_rust;
use zstd::stream::write::Encoder as ZstdEncoder;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer
const DEFAULT_BUFFER_SIZE: usize = 64 * 1024;

type BoxedWriter = Box<dyn Write + Send>;

/// Create a paraseq reader from optional path (stdin if None or "-")
fn create_paraseq_reader(path: Option<&str>) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    match path {
        None | Some("-") => {
            let stdin_reader = Box::new(std::io::stdin()) as Box<dyn std::io::Read + Send>;
            Reader::new(stdin_reader)
                .map_err(|e| anyhow::anyhow!("Failed to create stdin reader: {}", e))
        }
        Some(p) => {
            let (reader, _format) = niffler::send::from_path(p)
                .map_err(|e| anyhow::anyhow!("Failed to open file {}: {}", p, e))?;
            Reader::new(reader)
                .map_err(|e| anyhow::anyhow!("Failed to create reader for {}: {}", p, e))
        }
    }
}

/// Format a single record into a buffer (FASTA/FASTQ format)
fn format_record_to_buffer<R: Record>(
    record: &R,
    counter: u64,
    rename: bool,
    buffer: &mut Vec<u8>,
) -> Result<()> {
    let is_fasta = record.qual().is_none();

    // Header line
    buffer.write(if is_fasta { b">" } else { b"@" })?;
    if rename {
        buffer.extend_from_slice(counter.to_string().as_bytes());
    } else {
        buffer.extend_from_slice(&record.id());
    }
    buffer.write(b"\n")?;

    // Sequence line
    buffer.extend_from_slice(&record.seq_raw());

    if is_fasta {
        buffer.write(b"\n")?;
    } else {
        // FASTQ: plus line and quality
        buffer.write(b"\n+\n")?;
        if let Some(qual) = record.qual() {
            buffer.extend_from_slice(&qual);
        }
        buffer.write(b"\n")?;
    }
    Ok(())
}

/// Validate compression level for the given format
fn validate_compression_level(level: u8, min: u8, max: u8, format: &str) -> Result<()> {
    if level < min || level > max {
        Err(anyhow::anyhow!(
            "Invalid {} compression level {}. Must be between {} and {}.",
            format,
            level,
            min,
            max
        ))
    } else {
        Ok(())
    }
}

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: &str, compression_level: u8) -> Result<BoxedWriter> {
    if output_path == "-" {
        return Ok(Box::new(BufWriter::with_capacity(
            OUTPUT_BUFFER_SIZE,
            io::stdout(),
        )));
    }

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_path)
        .context(format!("Failed to create output file: {}", output_path))?;

    let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

    match output_path {
        p if p.ends_with(".gz") => {
            validate_compression_level(compression_level, 1, 9, "gzip")?;
            Ok(Box::new(GzEncoder::new(
                buffered_file,
                flate2::Compression::new(compression_level as u32),
            )))
        }
        p if p.ends_with(".zst") => {
            validate_compression_level(compression_level, 1, 22, "zstd")?;
            Ok(Box::new(ZstdEncoder::new(
                buffered_file,
                compression_level as i32,
            )?))
        }
        p if p.ends_with(".xz") => {
            validate_compression_level(compression_level, 0, 9, "xz")?;
            Ok(Box::new(XzEncoder::new(
                buffered_file,
                compression_level as u32,
            )))
        }
        _ => Ok(Box::new(buffered_file)),
    }
}

// JSON summary structure
#[derive(Serialize, Deserialize)]
pub struct FilterSummary {
    version: String,
    index: String,
    input: String,
    input2: Option<String>,
    output: String,
    output2: Option<String>,
    k: u8,
    w: u8,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    seqs_in: u64,
    seqs_out: u64,
    seqs_out_proportion: f64,
    seqs_removed: u64,
    seqs_removed_proportion: f64,
    bp_in: u64,
    bp_out: u64,
    bp_out_proportion: f64,
    bp_removed: u64,
    bp_removed_proportion: f64,
    time: f64,
    seqs_per_second: u64,
    bp_per_second: u64,
}

#[derive(Clone)]
struct FilterProcessor {
    // Minimizer matching parameters
    minimizer_hashes: Arc<FxHashSet<u64>>,
    kmer_length: u8,
    window_size: u8,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    deplete: bool,
    rename: bool,
    debug: bool,

    // Local buffers
    local_buffer: Vec<u8>,
    local_stats: ProcessingStats,

    // Global state
    global_writer: Arc<Mutex<BoxedWriter>>,
    global_writer2: Option<Arc<Mutex<BoxedWriter>>>,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    filtering_start_time: Instant,
}

#[derive(Clone, Default)]
struct ProcessingStats {
    total_seqs: u64,
    filtered_seqs: u64,
    total_bp: u64,
    output_bp: u64,
    filtered_bp: u64,
    output_seq_counter: u64,
}

impl FilterProcessor {
    /// Calculate required hits based on absolute and relative thresholds
    fn calculate_required_hits(&self, total_minimizers: usize) -> usize {
        let abs_required = self.abs_threshold;
        let rel_required = if total_minimizers == 0 {
            0
        } else {
            ((self.rel_threshold * total_minimizers as f64).round() as usize).max(1)
        };
        abs_required.max(rel_required)
    }

    /// Check if sequence meets filtering criteria
    fn meets_filtering_criteria(&self, hit_count: usize, total_minimizers: usize) -> bool {
        let required = self.calculate_required_hits(total_minimizers);
        if self.deplete {
            hit_count < required
        } else {
            hit_count >= required
        }
    }
    fn new(
        minimizer_hashes: Arc<FxHashSet<u64>>,
        kmer_length: u8,
        window_size: u8,
        abs_threshold: usize,
        rel_threshold: f64,
        prefix_length: usize,
        deplete: bool,
        rename: bool,
        debug: bool,
        writer: BoxedWriter,
        writer2: Option<BoxedWriter>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        filtering_start_time: Instant,
    ) -> Self {
        Self {
            minimizer_hashes,
            kmer_length,
            window_size,
            abs_threshold,
            rel_threshold,
            prefix_length,
            deplete,
            rename,
            debug,
            local_buffer: Vec::with_capacity(DEFAULT_BUFFER_SIZE),
            local_stats: ProcessingStats::default(),
            global_writer: Arc::new(Mutex::new(writer)),
            global_writer2: writer2.map(|w| Arc::new(Mutex::new(w))),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            filtering_start_time,
        }
    }

    fn should_keep_sequence(&self, seq: &[u8]) -> (bool, usize, usize, Vec<String>) {
        if seq.len() < self.kmer_length as usize {
            return (self.deplete, 0, 0, Vec::new()); // If too short, keep if in deplete mode
        }

        // Apply prefix length limit if specified
        let effective_seq = if self.prefix_length > 0 && seq.len() > self.prefix_length {
            &seq[..self.prefix_length]
        } else {
            seq
        };

        // Trim the last newline character from `effective_seq` if it has one.
        let effective_seq = effective_seq.strip_suffix(b"\n").unwrap_or(effective_seq);

        // Pack the sequence into 2-bit representation.
        // Any non-ACGT characters are silently converted to 2-bit ACGT as well.
        let packed_seq = packed_seq::PackedSeqVec::from_ascii(effective_seq);

        // TODO: Extract this to some nicer helper function in packed_seq?
        // TODO: Use SIMD?
        // TODO: Should probably add some test for this.
        // +2: one to round up, and one buffer.
        let mut invalid_mask = vec![0u64; packed_seq.len() / 64 + 2];
        for i in (0..effective_seq.len()).step_by(64) {
            let mut mask = 0;
            for (j, b) in effective_seq[i..(i + 64).min(effective_seq.len())]
                .iter()
                .enumerate()
            {
                mask |= ((!matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
                    as u64)
                    << j;
            }

            invalid_mask[i / 64] = mask;
        }

        let mut positions = Vec::new();
        simd_minimizers::canonical_minimizer_positions(
            packed_seq.as_slice(),
            self.kmer_length as usize,
            self.window_size as usize,
            &mut positions,
        );

        assert!(
            self.kmer_length <= 56,
            "Indexing the bitmask of invalid characters requires k<=56, but it is {}",
            self.kmer_length
        );

        // Filter positions to only include k-mers with ACGT bases
        positions.retain(|&pos| {
            // Extract bits pos .. pos+k from the bitmask.

            // mask of k ones in low positions.
            let mask = u64::MAX >> (64 - self.kmer_length);
            let byte = pos as usize / 8;
            let offset = pos as usize % 8;
            // The unaligned u64 read is OK, because we ensure that the underlying `Vec` always
            // has at least 8 bytes of padding at the end.
            let x =
                (unsafe { invalid_mask.as_ptr().byte_add(byte).read_unaligned() } >> offset) & mask;
            x == 0
        });

        // Get hash values for valid positions
        let minimizer_values: Vec<u64> = simd_minimizers::iter_canonical_minimizer_values(
            packed_seq.as_slice(),
            self.kmer_length as usize,
            &positions,
        )
        .map(|kmer| xxhash_rust::xxh3::xxh3_64(&kmer.to_le_bytes()))
        .collect();

        // Count distinct minimizer hits and collect matching k-mers
        let mut seen_hits = FxHashSet::default();
        let mut hit_count = 0;
        let mut hit_kmers = Vec::new();

        for (i, &hash) in minimizer_values.iter().enumerate() {
            if self.minimizer_hashes.contains(&hash) && seen_hits.insert(hash) {
                hit_count += 1;
                // Extract the k-mer sequence at this position
                if self.debug && i < positions.len() {
                    let pos = positions[i] as usize;
                    let kmer = &effective_seq[pos..pos + self.kmer_length as usize];
                    hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
                }
            }
        }

        (
            self.meets_filtering_criteria(hit_count, minimizer_values.len()),
            hit_count,
            minimizer_values.len(),
            hit_kmers,
        )
    }

    fn get_minimizer_hashes_and_positions(&self, seq: &[u8]) -> (Vec<u64>, Vec<u32>) {
        // Canonicalize sequence
        let canonical_seq = seq
            .iter()
            .map(|&b| match b {
                b'A' | b'a' => b'A',
                b'C' | b'c' => b'C',
                b'G' | b'g' => b'G',
                b'T' | b't' => b'T',
                _ => b'C',
            })
            .collect::<Vec<u8>>();

        let mut positions = Vec::new();
        simd_minimizers::canonical_minimizer_positions(
            packed_seq::AsciiSeq(&canonical_seq),
            self.kmer_length as usize,
            self.window_size as usize,
            &mut positions,
        );

        // Filter to valid positions
        let valid_positions: Vec<u32> = positions
            .into_iter()
            .filter(|&pos| {
                let pos_usize = pos as usize;
                if pos_usize + self.kmer_length as usize <= seq.len() {
                    let kmer = &seq[pos_usize..pos_usize + self.kmer_length as usize];
                    kmer.iter().all(|&b| {
                        matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
                    })
                } else {
                    false
                }
            })
            .collect();

        // Get hash values
        let hashes: Vec<u64> = simd_minimizers::iter_canonical_minimizer_values(
            packed_seq::AsciiSeq(&canonical_seq),
            self.kmer_length as usize,
            &valid_positions,
        )
        .map(|kmer| xxhash_rust::xxh3::xxh3_64(&kmer.to_le_bytes()))
        .collect();

        (hashes, valid_positions)
    }

    fn should_keep_pair(&self, seq1: &[u8], seq2: &[u8]) -> (bool, usize, usize, Vec<String>) {
        let mut all_hashes = Vec::new();
        let mut all_positions = Vec::new();
        let mut all_sequences = Vec::new();
        let mut seen_hits_pair = FxHashSet::default();
        let mut pair_hit_count = 0;
        let mut hit_kmers = Vec::new();

        // Process read 1
        if seq1.len() >= self.kmer_length as usize {
            let effective_seq = if self.prefix_length > 0 && seq1.len() > self.prefix_length {
                &seq1[..self.prefix_length]
            } else {
                seq1
            };

            let (hashes, positions) = self.get_minimizer_hashes_and_positions(effective_seq);
            all_hashes.extend(hashes);
            all_positions.extend(positions);
            all_sequences.extend(vec![effective_seq; all_hashes.len()]);
        }

        // Process read 2
        if seq2.len() >= self.kmer_length as usize {
            let effective_seq = if self.prefix_length > 0 && seq2.len() > self.prefix_length {
                &seq2[..self.prefix_length]
            } else {
                seq2
            };

            let (hashes, positions) = self.get_minimizer_hashes_and_positions(effective_seq);
            let start_idx = all_hashes.len();
            all_hashes.extend(hashes);
            all_positions.extend(positions);
            all_sequences.extend(vec![effective_seq; all_hashes.len() - start_idx]);
        }

        // Count hits and collect k-mers
        for (i, &hash) in all_hashes.iter().enumerate() {
            if self.minimizer_hashes.contains(&hash) && seen_hits_pair.insert(hash) {
                pair_hit_count += 1;
                if self.debug && i < all_positions.len() && i < all_sequences.len() {
                    let pos = all_positions[i] as usize;
                    let seq = all_sequences[i];
                    if pos + self.kmer_length as usize <= seq.len() {
                        let kmer = &seq[pos..pos + self.kmer_length as usize];
                        hit_kmers.push(String::from_utf8_lossy(kmer).to_string());
                    }
                }
            }
        }

        let total_minimizers = all_hashes.len();
        (
            self.meets_filtering_criteria(pair_hit_count, total_minimizers),
            pair_hit_count,
            total_minimizers,
            hit_kmers,
        )
    }

    fn write_record<Rf: Record>(&mut self, record: Rf) -> Result<()> {
        self.local_stats.output_seq_counter += 1;
        format_record_to_buffer(
            &record,
            self.local_stats.output_seq_counter,
            self.rename,
            &mut self.local_buffer,
        )
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.filtering_start_time.elapsed();
            let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            let output_seqs = stats.total_seqs - stats.filtered_seqs;
            let output_proportion = if stats.total_seqs > 0 {
                output_seqs as f64 / stats.total_seqs as f64
            } else {
                0.0
            };

            let output_bp_proportion = if stats.total_bp > 0 {
                stats.output_bp as f64 / stats.total_bp as f64
            } else {
                0.0
            };

            spinner.lock().set_message(format!(
                "Retained {}/{} sequences ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                output_seqs,
                stats.total_seqs,
                output_proportion * 100.0,
                stats.output_bp,
                stats.total_bp,
                output_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }
}

impl ParallelProcessor for FilterProcessor {
    fn process_record<Rf: Record>(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq_len = record.seq_raw().len();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq_len as u64;

        let (should_keep, hit_count, total_minimizers, hit_kmers) =
            self.should_keep_sequence(&record.seq_raw());

        // Show debug info for sequences with hits
        if self.debug && hit_count > 0 {
            eprintln!(
                "DEBUG: {} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(&record.id()),
                hit_count,
                total_minimizers,
                should_keep,
                hit_kmers.join(",")
            );
        }

        if should_keep {
            self.local_stats.output_bp += seq_len as u64;
            self.write_record(record)?;
        } else {
            self.local_stats.filtered_seqs += 1;
            self.local_stats.filtered_bp += seq_len as u64;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        if !self.local_buffer.is_empty() {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }

        // Clear buffer after releasing the lock for better performance
        self.local_buffer.clear();

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
            stats.output_seq_counter += self.local_stats.output_seq_counter;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

impl InterleavedParallelProcessor for FilterProcessor {
    fn process_interleaved_pair<Rf: Record>(
        &mut self,
        record1: Rf,
        record2: Rf,
    ) -> paraseq::parallel::Result<()> {
        let seq1_len = record1.seq_raw().len();
        let seq2_len = record2.seq_raw().len();
        self.local_stats.total_seqs += 2;
        self.local_stats.total_bp += (seq1_len + seq2_len) as u64;

        let (should_keep, hit_count, total_minimizers, hit_kmers) =
            self.should_keep_pair(&record1.seq_raw(), &record2.seq_raw());

        // Debug info for interleaved pairs
        if self.debug && hit_count > 0 {
            eprintln!(
                "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(&record1.id()),
                String::from_utf8_lossy(&record2.id()),
                hit_count,
                total_minimizers,
                should_keep,
                hit_kmers.join(",")
            );
        }

        if should_keep {
            self.local_stats.output_bp += (seq1_len + seq2_len) as u64;

            // Write both records to output
            self.write_record(record1)?;
            self.write_record(record2)?;
        } else {
            self.local_stats.filtered_seqs += 2;
            self.local_stats.filtered_bp += (seq1_len + seq2_len) as u64;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        if !self.local_buffer.is_empty() {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }

        // Clear buffer after releasing the lock for better performance
        self.local_buffer.clear();

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
            stats.output_seq_counter += self.local_stats.output_seq_counter;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

impl PairedParallelProcessor for FilterProcessor {
    fn process_record_pair<Rf: Record>(
        &mut self,
        record1: Rf,
        record2: Rf,
    ) -> paraseq::parallel::Result<()> {
        let seq1_len = record1.seq_raw().len();
        let seq2_len = record2.seq_raw().len();
        self.local_stats.total_seqs += 2;
        self.local_stats.total_bp += (seq1_len + seq2_len) as u64;

        let (should_keep, hit_count, total_minimizers, hit_kmers) =
            self.should_keep_pair(&record1.seq_raw(), &record2.seq_raw());

        // Debug info for paired reads
        if self.debug && hit_count > 0 {
            eprintln!(
                "DEBUG: {}/{} hits={}/{} keep={} kmers=[{}]",
                String::from_utf8_lossy(&record1.id()),
                String::from_utf8_lossy(&record2.id()),
                hit_count,
                total_minimizers,
                should_keep,
                hit_kmers.join(",")
            );
        }

        if should_keep {
            self.local_stats.output_bp += (seq1_len + seq2_len) as u64;

            // Write to appropriate writers
            if self.global_writer2.is_some() {
                // Separate output files
                self.write_record(record1)?;

                // Temporarily swap buffers for second writer
                let mut temp_buffer = Vec::with_capacity(DEFAULT_BUFFER_SIZE);
                std::mem::swap(&mut self.local_buffer, &mut temp_buffer);

                self.write_record(record2)?;

                // Write second record to writer2
                if let Some(ref writer2) = self.global_writer2 {
                    let mut global_writer2 = writer2.lock();
                    global_writer2.write_all(&self.local_buffer)?;
                }

                // Restore original buffer with first record
                self.local_buffer = temp_buffer;
            } else {
                // Interleaved output
                self.write_record(record1)?;
                self.write_record(record2)?;
            }
        } else {
            self.local_stats.filtered_seqs += 2;
            self.local_stats.filtered_bp += (seq1_len + seq2_len) as u64;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Write buffer to output
        if !self.local_buffer.is_empty() {
            let mut global_writer = self.global_writer.lock();
            global_writer.write_all(&self.local_buffer)?;
            global_writer.flush()?;
        }

        // Clear buffer after releasing the lock for better performance
        self.local_buffer.clear();

        // Also flush writer2 if it exists
        if let Some(ref writer2) = self.global_writer2 {
            let mut global_writer2 = writer2.lock();
            global_writer2.flush()?;
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.filtered_seqs += self.local_stats.filtered_seqs;
            stats.total_bp += self.local_stats.total_bp;
            stats.output_bp += self.local_stats.output_bp;
            stats.filtered_bp += self.local_stats.filtered_bp;
            stats.output_seq_counter += self.local_stats.output_seq_counter;
        }

        // Update spinner
        self.update_spinner();

        // Reset local stats
        self.local_stats = ProcessingStats::default();

        Ok(())
    }
}

pub fn run<P: AsRef<Path>>(
    minimizers_path: P,
    input_path: &str,
    input2_path: Option<&str>,
    output_path: &str,
    output2_path: Option<&str>,
    abs_threshold: usize,
    rel_threshold: f64,
    prefix_length: usize,
    summary_path: Option<&PathBuf>,
    deplete: bool,
    rename: bool,
    threads: usize,
    compression_level: u8,
    debug: bool,
    quiet: bool,
) -> Result<()> {
    let start_time = Instant::now();
    let version: String = env!("CARGO_PKG_VERSION").to_string();
    let tool_version = format!("deacon {}", version);

    // Enable quiet mode when debug is on to avoid mixed output
    let quiet = quiet || debug;

    // Configure thread pool if nonzero
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("Failed to initialize thread pool")?;
    }

    let mode = if deplete { "deplete" } else { "search" };

    let mut input_type = String::new();
    let mut options = Vec::<String>::new();
    let paired_stdin = input_path == "-" && input2_path.is_some() && input2_path.unwrap() == "-";
    if paired_stdin {
        input_type.push_str("interleaved");
    } else if let Some(_) = input2_path {
        input_type.push_str("paired");
    } else {
        input_type.push_str("single");
    }
    options.push(format!(
        "abs_threshold={}, rel_threshold={}",
        abs_threshold, rel_threshold
    ));
    if prefix_length > 0 {
        options.push(format!("prefix_length={}", prefix_length));
    }
    if rename {
        options.push("rename".to_string());
    }
    if threads > 0 {
        options.push(format!("threads={}", threads));
    }

    if !quiet {
        eprintln!(
            "Deacon v{}; mode: {}; input: {}; options: {}",
            version,
            mode,
            input_type,
            options.join(", ")
        );
    }

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) = load_minimizer_hashes(&minimizers_path)?;
    let minimizer_hashes = Arc::new(minimizer_hashes);

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    if !quiet {
        eprintln!(
            "Loaded index (k={}, w={}) in {:.2?}",
            kmer_length, window_size, load_time
        );
    }

    // Create the appropriate writer(s) based on the output path(s)
    let writer = get_writer(output_path, compression_level)?;
    let writer2 = if let (Some(output2), Some(_)) = (output2_path, input2_path) {
        Some(get_writer(output2, compression_level)?)
    } else {
        None
    };

    // Progress bar setup (only if not quiet)
    let spinner = if !quiet {
        let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{msg}")?,
        );
        pb.set_message("Filtering");
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    // Start timer for filtering rate calculation
    let filtering_start_time = Instant::now();

    // Create processor
    let processor = FilterProcessor::new(
        minimizer_hashes,
        kmer_length,
        window_size,
        abs_threshold,
        rel_threshold,
        prefix_length,
        deplete,
        rename,
        debug,
        writer,
        writer2,
        spinner.clone(),
        filtering_start_time,
    );

    // Process based on input type
    let num_threads = if threads == 0 {
        rayon::current_num_threads()
    } else {
        threads
    };

    if paired_stdin {
        // Interleaved paired from stdin - use native interleaved processor
        let reader = create_paraseq_reader(Some("-"))?;
        reader.process_parallel_interleaved(processor.clone(), num_threads)?;
    } else if let Some(input2_path) = input2_path {
        // Paired files
        let r1_reader = create_paraseq_reader(Some(input_path))?;
        let r2_reader = create_paraseq_reader(Some(input2_path))?;
        r1_reader.process_parallel_paired(r2_reader, processor.clone(), num_threads)?;
    } else {
        // Single file or stdin
        let reader = create_paraseq_reader(Some(input_path))?;
        reader.process_parallel(processor.clone(), num_threads)?;
    }

    // Get final stats
    let final_stats = processor.global_stats.lock();
    let total_seqs = final_stats.total_seqs;
    let filtered_seqs = final_stats.filtered_seqs;
    let total_bp = final_stats.total_bp;
    let output_bp = final_stats.output_bp;
    let filtered_bp = final_stats.filtered_bp;

    drop(final_stats); // Release lock

    // Flush writers - they should auto-flush on drop
    drop(processor.global_writer);
    if let Some(w2) = processor.global_writer2 {
        drop(w2);
    }

    let total_time = start_time.elapsed();
    let seqs_per_sec = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / total_time.as_secs_f64();
    let mbp_per_sec = bp_per_sec / 1_000_000.0;

    // Calculate proportions
    let filtered_proportion = if total_seqs > 0 {
        filtered_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    let filtered_bp_proportion = if total_bp > 0 {
        filtered_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    let output_seqs = total_seqs - filtered_seqs;
    let output_seq_proportion = if total_seqs > 0 {
        output_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    let output_bp_proportion = if total_bp > 0 {
        output_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    if !quiet {
        eprintln!(
            "Retained {}/{} sequences ({:.3}%), {}/{} bp ({:.3}%)",
            output_seqs,
            total_seqs,
            output_seq_proportion * 100.0,
            output_bp,
            total_bp,
            output_bp_proportion * 100.0
        );

        eprintln!(
            "Completed in {:.2?}. Speed: {:.0} seqs/s ({:.1} Mbp/s)",
            total_time, seqs_per_sec, mbp_per_sec
        );
    }

    // Finish and clear spinner - disable it completely
    if let Some(ref spinner) = spinner {
        let pb = spinner.lock();
        pb.set_draw_target(ProgressDrawTarget::hidden());
        pb.finish_and_clear();
    }

    // Build and write JSON summary if path provided
    if let Some(summary_file) = summary_path {
        let seqs_out = total_seqs - filtered_seqs;

        let summary = FilterSummary {
            version: tool_version,
            index: minimizers_path.as_ref().to_string_lossy().to_string(),
            input: input_path.to_string(),
            input2: input2_path.map(|s| s.to_string()),
            output: output_path.to_string(),
            output2: output2_path.map(|s| s.to_string()),
            k: kmer_length,
            w: window_size,
            abs_threshold,
            rel_threshold,
            prefix_length,
            deplete,
            rename,
            seqs_in: total_seqs as u64,
            seqs_out: seqs_out as u64,
            seqs_out_proportion: output_seq_proportion,
            seqs_removed: filtered_seqs as u64,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp as u64,
            bp_out: output_bp as u64,
            bp_out_proportion: output_bp_proportion,
            bp_removed: filtered_bp as u64,
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        let file = File::create(summary_file)
            .context(format!("Failed to create summary: {:?}", summary_file))?;
        let writer = BufWriter::new(file);

        serde_json::to_writer_pretty(writer, &summary).context("Failed to write summary")?;

        if !quiet {
            eprintln!("Summary saved to {:?}", summary_file);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_summary() {
        let summary = FilterSummary {
            version: "deacon 0.1.0".to_string(),
            index: "test.idx".to_string(),
            input: "test.fastq".to_string(),
            input2: Some("test2.fastq".to_string()),
            output: "output.fastq".to_string(),
            output2: Some("output2.fastq".to_string()),
            k: 31,
            w: 21,
            abs_threshold: 1,
            rel_threshold: 0.01,
            prefix_length: 0,
            deplete: false,
            rename: false,
            seqs_in: 100,
            seqs_out: 90,
            seqs_out_proportion: 0.9,
            seqs_removed: 10,
            seqs_removed_proportion: 0.1,
            bp_in: 10000,
            bp_out: 9000,
            bp_out_proportion: 0.9,
            bp_removed: 1000,
            bp_removed_proportion: 0.1,
            time: 1.5,
            seqs_per_second: 66,
            bp_per_second: 6666,
        };

        let json = serde_json::to_string(&summary).unwrap();
        let parsed: FilterSummary = serde_json::from_str(&json).unwrap();

        assert_eq!(parsed.version, "deacon 0.1.0");
        assert_eq!(parsed.seqs_in, 100);
        assert_eq!(parsed.seqs_removed_proportion, 0.1);
        assert_eq!(parsed.seqs_out_proportion, 0.9);
        assert_eq!(parsed.bp_out_proportion, 0.9);
        assert_eq!(parsed.input, "test.fastq");
        assert_eq!(parsed.input2, Some("test2.fastq".to_string()));
        assert_eq!(parsed.output, "output.fastq");
        assert_eq!(parsed.output2, Some("output2.fastq".to_string()));
    }
}

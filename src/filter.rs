use crate::index::load_minimizer_hashes;
use crate::minimizers::fill_minimizer_hashes;
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use needletail::parse_fastx_file;
use needletail::parse_fastx_stdin;
use needletail::parser::{Format, SequenceRecord};
use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use zstd::stream::write::Encoder as ZstdEncoder;

const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // Opt: 8MB output buffer

trait FastxWriter: Write {
    fn flush_all(&mut self) -> io::Result<()>;
}

struct StandardWriter<W: Write>(W);

impl<W: Write> Write for StandardWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write> FastxWriter for StandardWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        self.flush()
    }
}

struct GzipWriter<W: Write>(GzEncoder<W>);

impl<W: Write> Write for GzipWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write> FastxWriter for GzipWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        self.flush()?;
        self.0.try_finish()?;
        Ok(())
    }
}

struct ZstdWriter<W: Write> {
    encoder: Option<ZstdEncoder<'static, W>>,
}

impl<W: Write> ZstdWriter<W> {
    fn new(writer: W, compression_level: i32) -> Result<Self> {
        let encoder =
            ZstdEncoder::new(writer, compression_level).context("Failed to create zstd encoder")?;
        Ok(Self {
            encoder: Some(encoder),
        })
    }
}

impl<W: Write> Write for ZstdWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(encoder) = &mut self.encoder {
            encoder.write(buf)
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(encoder) = &mut self.encoder {
            encoder.flush()
        } else {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "Writer has been closed",
            ))
        }
    }
}

impl<W: Write> FastxWriter for ZstdWriter<W> {
    fn flush_all(&mut self) -> io::Result<()> {
        if let Some(encoder) = self.encoder.take() {
            // Take ownership of the encoder
            encoder.finish()?;
        }
        Ok(())
    }
}

// Return a file writer appropriate for the output path extension
fn get_writer(output_path: &str) -> Result<Box<dyn FastxWriter>> {
    if output_path == "-" {
        // Write to stdout
        let stdout = io::stdout();
        let writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, stdout);
        Ok(Box::new(StandardWriter(writer)))
    } else {
        // Write to file with extension-appropriate encoder
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(output_path)
            .context(format!("Failed to create output file: {}", output_path))?;

        let buffered_file = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

        if output_path.ends_with(".gz") {
            // Use gzip
            let encoder = GzEncoder::new(buffered_file, Compression::fast());
            Ok(Box::new(GzipWriter(encoder)))
        } else if output_path.ends_with(".zst") {
            // Use zstd
            let writer =
                ZstdWriter::new(buffered_file, 1).context("Failed to create zstd encoder")?;
            Ok(Box::new(writer))
        } else {
            // Raw (dog) FASTQ
            Ok(Box::new(StandardWriter(buffered_file)))
        }
    }
}

// JSON report structure
#[derive(Serialize, Deserialize)]
pub struct FilterReport {
    version: String,
    index: String,
    input: String,
    output: String,
    k: usize,
    w: usize,
    m: usize,
    n: usize,
    invert: bool,
    rename: bool,
    seqs_in: u64,
    seqs_out: u64,
    seqs_removed: u64,
    seqs_removed_proportion: f64,
    bp_in: u64,
    bp_out: u64,
    bp_removed: u64,
    bp_removed_proportion: f64,
    time: f64,
    seqs_per_second: u64,
    bp_per_second: u64,
}

pub fn run<P: AsRef<Path>>(
    minimizers_path: P,
    input_path: &str,
    output_path: &str,
    min_matches: usize,
    prefix_length: usize,
    report_path: Option<&PathBuf>,
    invert: bool,
    rename: bool,
) -> Result<()> {
    let start_time = Instant::now();

    // Load minimizers hashes and parse header
    let (minimizer_hashes, header) = load_minimizer_hashes(&minimizers_path)?;

    let kmer_length = header.kmer_length();
    let window_size = header.window_size();

    let load_time = start_time.elapsed();
    eprintln!(
        "Loaded index (k={}, w={}) in {:.2?}",
        kmer_length, window_size, load_time
    );

    if invert {
        eprintln!("Invert mode: keeping sequences with matches");
    }

    if rename {
        eprintln!("Rename mode: replacing sequence headers with sequential numbers");
    }

    // Create the appropriate writer based on the output path
    let mut writer = get_writer(output_path)?;

    // A progress bar would require a denominator, so let's spin
    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::default_spinner()
            .tick_strings(&[".  ", ".. ", "...", " ..", "  .", "   "])
            .template("{msg}{spinner} ")?,
    );
    spinner.set_message("Filtering");
    spinner.enable_steady_tick(Duration::from_millis(100)); // Update once per second

    // Init counters
    let mut total_seqs = 0;
    let mut filtered_seqs = 0;
    let mut total_bp = 0;
    let mut output_bp = 0;
    let mut filtered_bp = 0; // Track filtered base pairs explicitly
    let mut output_seq_counter = 0;

    // Opt: pre-allocate buffers for reuse
    let mut minimizer_buffer = Vec::with_capacity(64);
    let mut output_record_buffer = Vec::with_capacity(1024);
    let mut seen_hits = FxHashSet::default();

    // Parse FASTX from input (file or stdin)
    let mut reader = if input_path == "-" {
        parse_fastx_stdin()?
    } else {
        parse_fastx_file(input_path)?
    };

    while let Some(record) = reader.next() {
        let record = record?;

        let seq = record.seq();
        total_seqs += 1;
        total_bp += seq.len();

        // Check for minimizer hits
        let mut hit_count = 0;

        if seq.len() >= kmer_length {
            // Reuse pre-allocated buffers
            minimizer_buffer.clear();
            seen_hits.clear();

            // Apply prefix length limit if specified (if > 0)
            let effective_seq = if prefix_length > 0 && seq.len() > prefix_length {
                &seq.as_ref()[..prefix_length]
            } else {
                seq.as_ref()
            };

            // Get minimizer hash values using parameters from header
            fill_minimizer_hashes(
                effective_seq,
                kmer_length,
                window_size,
                &mut minimizer_buffer,
            );

            // Count distinct minimizer hits
            for &hash in &minimizer_buffer {
                if minimizer_hashes.contains(&hash) && seen_hits.insert(hash) {
                    hit_count += 1;
                    if hit_count >= min_matches {
                        break;
                    }
                }
            }
        }

        // Determine if we should output this sequence based on hit count and invert flag
        let should_output = if !invert {
            // When not inverted, keep sequences with fewer than min_matches
            hit_count < min_matches
        } else {
            // When inverted, keep sequences with greater than or equal to min_matches
            hit_count >= min_matches
        };

        if should_output {
            // Track output base pairs
            output_bp += seq.len();

            // Format as FASTX directly to a byte buffer
            output_record_buffer.clear();

            // Increment output sequence counter
            output_seq_counter += 1;

            // Write the record based on format, with optional renaming
            output_fastx_record(
                &record,
                &mut output_record_buffer,
                rename,
                output_seq_counter,
            );

            // Write the whole record at once
            writer.write_all(&output_record_buffer)?;
        } else {
            filtered_seqs += 1;
            filtered_bp += seq.len(); // Track filtered base pairs
        }

        // Periodically flush and update spinner message with detailed stats
        if total_seqs % (100_000) == 0 {
            writer.flush()?;

            let elapsed = start_time.elapsed();
            let seqs_per_sec = total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = total_bp as f64 / elapsed.as_secs_f64();
            let mbp_per_sec = bp_per_sec / 1_000_000.0;

            // Calculate filtered proportion directly
            let filtered_proportion = if total_seqs > 0 {
                filtered_seqs as f64 / total_seqs as f64
            } else {
                0.0
            };

            // Calculate filtered base pair proportion
            let filtered_bp_proportion = if total_bp > 0 {
                filtered_bp as f64 / total_bp as f64
            } else {
                0.0
            };

            // Update spinner message with detailed stats
            spinner.set_message(format!(
                "Removed {}/{} seqs ({:.2}%), {}/{} bp ({:.2}%). {:.0} seqs/s ({:.1} Mbp/s)",
                filtered_seqs,
                total_seqs,
                filtered_proportion * 100.0, // Convert to percentage for display
                filtered_bp,
                total_bp,
                filtered_bp_proportion * 100.0,
                seqs_per_sec,
                mbp_per_sec
            ));
        }
    }

    // Ensure everything is flushed
    // What to do about this for client server imp?
    writer.flush_all()?;

    let total_time = start_time.elapsed();
    let seqs_per_sec = total_seqs as f64 / total_time.as_secs_f64();
    let bp_per_sec = total_bp as f64 / total_time.as_secs_f64();
    // Convert to Mbp/s for display only
    let mbp_per_sec = bp_per_sec / 1_000_000.0;

    // Calculate filtered proportion directly
    let filtered_proportion = if total_seqs > 0 {
        filtered_seqs as f64 / total_seqs as f64
    } else {
        0.0
    };

    // Calculate filtered base pair proportion
    let filtered_bp_proportion = if total_bp > 0 {
        filtered_bp as f64 / total_bp as f64
    } else {
        0.0
    };

    // Finish spinner with final message
    spinner.finish_with_message(format!(
        "Removed {}/{} sequences ({:.3}%) , {}/{} bp ({:.3}%)",
        filtered_seqs,
        total_seqs,
        filtered_proportion * 100.0, // Convert to percentage for display
        filtered_bp,
        total_bp,
        filtered_bp_proportion * 100.0
    ));

    // Print completion message without spinner
    eprintln!(
        "Completed in {:.2?}. Speed: {:.0} seqs/s ({:.1} Mbp/s)",
        total_time, seqs_per_sec, mbp_per_sec
    );

    // Build and write a JSON report if path provided
    if let Some(report_file) = report_path {
        // Get number of sequences passing filter
        let seqs_out = total_seqs - filtered_seqs;

        let report = FilterReport {
            version: env!("CARGO_PKG_VERSION").to_string(),
            index: minimizers_path.as_ref().to_string_lossy().to_string(),
            input: input_path.to_string(),
            output: output_path.to_string(),
            k: kmer_length,
            w: window_size,
            m: min_matches,
            n: prefix_length,
            invert,
            rename,
            seqs_in: total_seqs as u64,
            seqs_out: seqs_out as u64,
            seqs_removed: filtered_seqs as u64,
            seqs_removed_proportion: filtered_proportion,
            bp_in: total_bp as u64,
            bp_out: output_bp as u64,
            bp_removed: filtered_bp as u64, // Use directly tracked filtered_bp
            bp_removed_proportion: filtered_bp_proportion,
            time: total_time.as_secs_f64(),
            seqs_per_second: seqs_per_sec as u64,
            bp_per_second: bp_per_sec as u64,
        };

        // Write report file
        let file = File::create(report_file)
            .context(format!("Failed to create report file: {:?}", report_file))?;
        let writer = BufWriter::new(file);

        // Serialize and write the report as JSON
        serde_json::to_writer_pretty(writer, &report).context("Failed to write JSON report")?;

        eprintln!("JSON report saved to: {:?}", report_file);
    }

    Ok(())
}

/// Push FASTA or FASTQ record to output buffer
fn output_fastx_record(
    record: &SequenceRecord,
    buffer: &mut Vec<u8>,
    rename: bool,
    seq_number: u64,
) {
    match record.format() {
        Format::Fasta => {
            buffer.push(b'>');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(record.id());
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(record.seq().as_ref());
            buffer.push(b'\n');
        }
        Format::Fastq => {
            buffer.push(b'@');
            if rename {
                // Use sequential numbering for sequence ID
                buffer.extend_from_slice(seq_number.to_string().as_bytes());
            } else {
                // Use original sequence ID
                buffer.extend_from_slice(record.id());
            }
            buffer.push(b'\n');
            buffer.extend_from_slice(record.seq().as_ref());
            buffer.extend_from_slice(b"\n+\n");
            if let Some(qual) = record.qual() {
                buffer.extend_from_slice(qual);
            }
            buffer.push(b'\n');
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::write_minimizers;
    use crate::index_format::IndexHeader;
    use tempfile::TempDir;

    #[allow(dead_code)] // Suppress unused warnings
    fn create_test_index() -> (PathBuf, IndexHeader, TempDir) {
        // Create a temporary directory that will be automatically cleaned up
        let temp_dir = TempDir::new().unwrap();
        let index_path = temp_dir.path().join("test.idx");

        // Create dummy minimizers
        let minimizers: FxHashSet<u64> = [1, 2, 3, 4, 5].iter().cloned().collect();
        let header = IndexHeader::new(5, 3);

        // Write to file
        write_minimizers(&minimizers, &header, Some(&index_path)).unwrap();

        // Return the TempDir along with the other values so it stays in scope
        (index_path, header, temp_dir)
    }

    #[test]
    fn test_filter_report() {
        // Create a sample report
        let report = FilterReport {
            version: "0.1.0".to_string(),
            index: "test.idx".to_string(),
            input: "test.fastq".to_string(),
            output: "output.fastq".to_string(),
            k: 31,
            w: 21,
            m: 1,
            n: 0,
            invert: false,
            rename: false,
            seqs_in: 100,
            seqs_out: 90,
            seqs_removed: 10,
            seqs_removed_proportion: 0.1,
            bp_in: 10000,
            bp_out: 9000,
            bp_removed: 1000,
            bp_removed_proportion: 0.1,
            time: 1.5,
            seqs_per_second: 66,
            bp_per_second: 6666,
        };

        // Test JSON ser+de
        let json = serde_json::to_string(&report).unwrap();
        let parsed: FilterReport = serde_json::from_str(&json).unwrap();

        // Check values
        assert_eq!(parsed.version, "0.1.0");
        assert_eq!(parsed.seqs_in, 100);
        assert_eq!(parsed.seqs_removed_proportion, 0.1);
        assert_eq!(parsed.input, "test.fastq");
        assert_eq!(parsed.output, "output.fastq");
    }
}

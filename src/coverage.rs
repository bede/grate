use crate::minimizers::{
    fill_minimizers, fill_minimizers_with_positions, Buffers, KmerHasher, MinimizerVec,
};
use anyhow::Result;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use paraseq::fastx::Reader;
use paraseq::parallel::{ParallelProcessor, ParallelReader};
use paraseq::Record;
use parking_lot::Mutex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::BuildHasher;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

/// BuildHasher using rapidhash with fixed seed for fast init
#[derive(Clone, Default)]
pub struct FixedRapidHasher;

impl BuildHasher for FixedRapidHasher {
    type Hasher = rapidhash::fast::RapidHasher<'static>;

    fn build_hasher(&self) -> Self::Hasher {
        rapidhash::fast::SeedableState::fixed().build_hasher()
    }
}

/// RapidHashSet using rapidhash with fixed seed for fast init
pub type RapidHashSet<T> = HashSet<T, FixedRapidHasher>;

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
pub enum MinimizerSet {
    U64(RapidHashSet<u64>),
    U128(RapidHashSet<u128>),
}

impl MinimizerSet {
    pub fn len(&self) -> usize {
        match self {
            MinimizerSet::U64(set) => set.len(),
            MinimizerSet::U128(set) => set.len(),
        }
    }

    pub fn is_u64(&self) -> bool {
        matches!(self, MinimizerSet::U64(_))
    }
}

#[derive(Debug, Clone)]
pub struct TargetInfo {
    pub name: String,
    pub length: usize,
    pub minimizers: MinimizerVec,
    pub minimizer_positions: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageResult {
    pub target: String,
    pub length: usize,
    pub total_minimizers: usize,
    pub covered_minimizers: usize,
    pub coverage_fraction: f64,
    pub mean_depth: f64,
    pub median_depth: f64,
    pub depth_histogram: Vec<(usize, usize)>, // (depth, count)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageReport {
    pub version: String,
    pub reference_file: String,
    pub reads_file: String,
    pub parameters: CoverageParameters,
    pub targets: Vec<CoverageResult>,
    pub overall_stats: OverallStats,
    pub timing: TimingStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageParameters {
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverallStats {
    pub total_targets: usize,
    pub total_minimizers: usize,
    pub total_covered_minimizers: usize,
    pub overall_coverage_fraction: f64,
    pub overall_mean_depth: f64,
    pub overall_median_depth: f64,
    pub total_reads_processed: u64,
    pub total_bp_processed: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingStats {
    pub reference_processing_time: f64,
    pub reads_processing_time: f64,
    pub analysis_time: f64,
    pub total_time: f64,
    pub reads_per_second: f64,
    pub bp_per_second: f64,
}

/// Output format options
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Table,
    Csv,
    Json,
}

pub struct CoverageConfig {
    pub reference_path: PathBuf,
    pub reads_path: PathBuf,
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub output_format: OutputFormat,
}

impl CoverageConfig {
    pub fn execute(&self) -> Result<()> {
        run_coverage_analysis(self)
    }
}

fn reader_with_inferred_batch_size(
    in_path: Option<&Path>,
) -> Result<Reader<Box<dyn std::io::Read + Send>>> {
    let mut reader = paraseq::fastx::Reader::from_optional_path(in_path)?;
    reader.update_batch_size_in_bp(256 * 1024)?;
    Ok(reader)
}

/// Processor for collecting reference targets with minimizers
#[derive(Clone)]
struct ReferenceProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    positions: Vec<usize>,
    targets: Arc<Mutex<Vec<TargetInfo>>>,
}

impl ReferenceProcessor {
    fn new(kmer_length: u8, window_size: u8) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        Self {
            kmer_length,
            window_size,
            hasher: KmerHasher::new(kmer_length as usize),
            buffers,
            positions: Vec::new(),
            targets: Arc::new(Mutex::new(Vec::new())),
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ReferenceProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let sequence = record.seq();
        let target_name = String::from_utf8_lossy(record.id()).to_string();

        fill_minimizers_with_positions(
            &sequence,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            &mut self.buffers,
            &mut self.positions,
        );

        // Clone minimizers and positions
        let minimizers = self.buffers.minimizers.clone();
        let minimizer_positions = self.positions.clone();

        self.targets.lock().push(TargetInfo {
            name: target_name,
            length: sequence.len(),
            minimizers,
            minimizer_positions,
        });

        Ok(())
    }
}

fn process_reference_file(
    reference_path: &Path,
    kmer_length: u8,
    window_size: u8,
) -> Result<Vec<TargetInfo>> {

    let in_path = if reference_path.to_string_lossy() == "-" {
        None
    } else {
        Some(reference_path)
    };

    let reader = reader_with_inferred_batch_size(in_path)?;
    let mut processor = ReferenceProcessor::new(kmer_length, window_size);

    // Single thread to preserve order?
    reader.process_parallel(&mut processor, 1)?;

    let targets = Arc::try_unwrap(processor.targets).unwrap().into_inner();

    Ok(targets)
}

#[derive(Clone, Default, Debug)]
struct ProcessingStats {
    total_seqs: u64,
    total_bp: u64,
    last_reported: u64,
}

/// Processor for counting minimizer depths from reads
#[derive(Clone)]
struct ReadsProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    reference_minimizers: Arc<MinimizerSet>,

    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_counts_u64: Option<HashMap<u64, usize>>,
    local_counts_u128: Option<HashMap<u128, usize>>,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_counts_u64: Arc<Mutex<Option<HashMap<u64, usize>>>>,
    global_counts_u128: Arc<Mutex<Option<HashMap<u128, usize>>>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
}

impl ReadsProcessor {
    fn new(
        kmer_length: u8,
        window_size: u8,
        reference_minimizers: Arc<MinimizerSet>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
    ) -> Self {
        let buffers = if kmer_length <= 32 {
            Buffers::new_u64()
        } else {
            Buffers::new_u128()
        };

        let (local_counts_u64, local_counts_u128) = if kmer_length <= 32 {
            (Some(HashMap::new()), None)
        } else {
            (None, Some(HashMap::new()))
        };

        let (global_counts_u64, global_counts_u128) = if kmer_length <= 32 {
            (
                Arc::new(Mutex::new(Some(HashMap::new()))),
                Arc::new(Mutex::new(None)),
            )
        } else {
            (
                Arc::new(Mutex::new(None)),
                Arc::new(Mutex::new(Some(HashMap::new()))),
            )
        };

        Self {
            kmer_length,
            window_size,
            hasher: KmerHasher::new(kmer_length as usize),
            reference_minimizers,
            buffers,
            local_stats: ProcessingStats::default(),
            local_counts_u64,
            local_counts_u128,
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_counts_u64,
            global_counts_u128,
            spinner,
            start_time,
        }
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let reads_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Processing reads: {} reads ({}). {:.0} reads/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                reads_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for ReadsProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let seq = record.seq();
        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += seq.len() as u64;

        fill_minimizers(
            &seq,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            &mut self.buffers,
        );

        // Count minimizers in target ref
        match (&self.buffers.minimizers, &*self.reference_minimizers) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(ref_set)) => {
                let local_counts = self.local_counts_u64.as_mut().unwrap();
                for &minimizer in vec {
                    if ref_set.contains(&minimizer) {
                        *local_counts.entry(minimizer).or_insert(0) += 1;
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(ref_set)) => {
                let local_counts = self.local_counts_u128.as_mut().unwrap();
                for &minimizer in vec {
                    if ref_set.contains(&minimizer) {
                        *local_counts.entry(minimizer).or_insert(0) += 1;
                    }
                }
            }
            _ => panic!("Mismatch between MinimizerVec and MinimizerSet types"),
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Merge local into global counts
        if let Some(local) = &mut self.local_counts_u64 {
            let mut global = self.global_counts_u64.lock();
            let global_map = global.as_mut().unwrap();
            for (&minimizer, &count) in local.iter() {
                *global_map.entry(minimizer).or_insert(0) += count;
            }
            local.clear();
        } else {
            let mut global = self.global_counts_u128.lock();
            let global_map = global.as_mut().unwrap();
            let local = self.local_counts_u128.as_mut().unwrap();
            for (&minimizer, &count) in local.iter() {
                *global_map.entry(minimizer).or_insert(0) += count;
            }
            local.clear();
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            // Update spinner every GB
            let current_gb = stats.total_bp / 1_000_000_000;
            if current_gb > stats.last_reported {
                drop(stats); // Release lock before updating spinner
                self.update_spinner();
                self.global_stats.lock().last_reported = current_gb;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Enum to abstract over u64 and u128 depth maps
enum DepthMap {
    U64(HashMap<u64, usize>),
    U128(HashMap<u128, usize>),
}

fn process_reads_file(
    reads_path: &Path,
    reference_minimizers: Arc<MinimizerSet>,
    kmer_length: u8,
    window_size: u8,
    threads: usize,
    quiet: bool,
) -> Result<(DepthMap, u64, u64)> {
    let in_path = if reads_path.to_string_lossy() == "-" {
        None
    } else {
        Some(reads_path)
    };
    let reader = reader_with_inferred_batch_size(in_path)?;

    // Progress bar
    let spinner = if !quiet {
        let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{msg}")?,
        );
        pb.set_message("Processing reads: 0 reads (0bp)");
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    let start_time = Instant::now();
    let mut processor = ReadsProcessor::new(
        kmer_length,
        window_size,
        reference_minimizers,
        spinner.clone(),
        start_time,
    );

    reader.process_parallel(&mut processor, threads)?;

    // Finish spinner
    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    let stats = processor.global_stats.lock().clone();
    let depth_map = if kmer_length <= 32 {
        let map = Arc::try_unwrap(processor.global_counts_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        DepthMap::U64(map)
    } else {
        let map = Arc::try_unwrap(processor.global_counts_u128)
            .unwrap()
            .into_inner()
            .unwrap();
        DepthMap::U128(map)
    };

    if !quiet {
        let elapsed = start_time.elapsed();
        let unique_minimizers = match &depth_map {
            DepthMap::U64(m) => m.len(),
            DepthMap::U128(m) => m.len(),
        };
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
        eprintln!(
            "Reads: {} records, {}, {} unique minimizers ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            unique_minimizers,
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((depth_map, stats.total_seqs, stats.total_bp))
}

fn calculate_coverage_statistics(
    targets: &[TargetInfo],
    depth_map: &DepthMap,
) -> Vec<CoverageResult> {
    targets
        .iter()
        .map(|target| {
            let total_minimizers = target.minimizers.len();
            let mut depths = Vec::new();
            let mut covered_count = 0;

            // Collect depths for all minimizers in this target
            match (&target.minimizers, depth_map) {
                (MinimizerVec::U64(vec), DepthMap::U64(map)) => {
                    for &minimizer in vec {
                        let depth = map.get(&minimizer).copied().unwrap_or(0);
                        depths.push(depth);
                        if depth > 0 {
                            covered_count += 1;
                        }
                    }
                }
                (MinimizerVec::U128(vec), DepthMap::U128(map)) => {
                    for &minimizer in vec {
                        let depth = map.get(&minimizer).copied().unwrap_or(0);
                        depths.push(depth);
                        if depth > 0 {
                            covered_count += 1;
                        }
                    }
                }
                _ => panic!("Mismatch between MinimizerVec and DepthMap types"),
            }

            let coverage_fraction = if total_minimizers > 0 {
                covered_count as f64 / total_minimizers as f64
            } else {
                0.0
            };

            let mean_depth = if !depths.is_empty() {
                depths.iter().sum::<usize>() as f64 / depths.len() as f64
            } else {
                0.0
            };

            // Calculate median depth
            let mut sorted_depths = depths.clone();
            sorted_depths.sort_unstable();
            let median_depth = if sorted_depths.is_empty() {
                0.0
            } else if sorted_depths.len() % 2 == 0 {
                let mid = sorted_depths.len() / 2;
                (sorted_depths[mid - 1] + sorted_depths[mid]) as f64 / 2.0
            } else {
                sorted_depths[sorted_depths.len() / 2] as f64
            };

            // Depth hist
            let mut depth_counts: HashMap<usize, usize> = HashMap::new();
            for depth in depths {
                *depth_counts.entry(depth).or_insert(0) += 1;
            }
            let mut depth_histogram: Vec<(usize, usize)> = depth_counts.into_iter().collect();
            depth_histogram.sort_by_key(|(depth, _)| *depth);

            CoverageResult {
                target: target.name.clone(),
                length: target.length,
                total_minimizers,
                covered_minimizers: covered_count,
                coverage_fraction,
                mean_depth,
                median_depth,
                depth_histogram,
            }
        })
        .collect()
}

fn output_results(
    report: &CoverageReport,
    output_path: Option<&PathBuf>,
    output_format: OutputFormat,
) -> Result<()> {
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = writer;

    match output_format {
        OutputFormat::Json => {
            // Output JSON
            serde_json::to_writer_pretty(&mut writer, report)?;
            writeln!(writer)?;
        }
        OutputFormat::Csv => {
            // Output CSV with headers
            writeln!(
                writer,
                "target,length_bp,total_minimizers,covered_minimizers,coverage_fraction,mean_depth,median_depth"
            )?;

            for result in &report.targets {
                writeln!(
                    writer,
                    "{},{},{},{},{:.4},{:.2},{:.2}",
                    result.target,
                    result.length,
                    result.total_minimizers,
                    result.covered_minimizers,
                    result.coverage_fraction,
                    result.mean_depth,
                    result.median_depth,
                )?;
            }

            // Add overall statistics as a final row
            let total_length: usize = report.targets.iter().map(|r| r.length).sum();

            writeln!(
                writer,
                "ALL,{},{},{},{:.4},{:.2},{:.2}",
                total_length,
                report.overall_stats.total_minimizers,
                report.overall_stats.total_covered_minimizers,
                report.overall_stats.overall_coverage_fraction,
                report.overall_stats.overall_mean_depth,
                report.overall_stats.overall_median_depth,
            )?;
        }
        OutputFormat::Table => {
            // Output human-readable table with headers
            writeln!(writer)?;
            writeln!(
                writer,
                "{:<30} | {:>10} | {:>12} | {:>12} | {:>12} | {:>12}",
                "Target", "Length", "Minimizers", "Coverage", "Mean depth", "Median depth"
            )?;
            writeln!(
                writer,
                "{}",
                "-".repeat(20 + 3 + 10 + 3 + 21 + 3 + 15 + 3 + 17 + 3 + 19)
            )?;

            for result in &report.targets {
                writeln!(
                    writer,
                    "{:<30} | {:>10} | {:>12} | {:>11.1}% | {:>11.1}x | {:>11.1}x",
                    truncate_string(&result.target, 30),
                    format_bp(result.length),
                    result.total_minimizers,
                    result.coverage_fraction * 100.0,
                    result.mean_depth,
                    result.median_depth
                )?;
            }

            writeln!(writer)?;
            writeln!(
                writer,
                "Overall: Coverage: {:.1}%, Median depth: {:.1}x, Mean depth: {:.1}x",
                report.overall_stats.overall_coverage_fraction * 100.0,
                report.overall_stats.overall_median_depth,
                report.overall_stats.overall_mean_depth
            )?;
        }
    }

    Ok(())
}

fn truncate_string(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len - 3])
    }
}

fn format_bp(bp: usize) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.1}Gbp", bp as f64 / 1_000_000_000.0)
    } else if bp >= 1_000_000 {
        format!("{:.1}Mbp", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}Kbp", bp as f64 / 1_000.0)
    } else {
        format!("{}bp", bp)
    }
}

fn format_bp_per_sec(bp_per_sec: f64) -> String {
    if bp_per_sec >= 1_000_000_000.0 {
        format!("{:.1} Gbp/s", bp_per_sec / 1_000_000_000.0)
    } else if bp_per_sec >= 1_000_000.0 {
        format!("{:.1} Mbp/s", bp_per_sec / 1_000_000.0)
    } else if bp_per_sec >= 1_000.0 {
        format!("{:.1} Kbp/s", bp_per_sec / 1_000.0)
    } else {
        format!("{:.0} bp/s", bp_per_sec)
    }
}

pub fn run_coverage_analysis(config: &CoverageConfig) -> Result<()> {
    let start_time = Instant::now();
    let version = env!("CARGO_PKG_VERSION").to_string();

    if !config.quiet {
        eprintln!(
            "Grate v{}; mode: cov; options: k={}, w={}, threads={}",
            version, config.kmer_length, config.window_size, config.threads
        );
    }

    // Process reference file
    let ref_start = Instant::now();
    let targets = process_reference_file(
        &config.reference_path,
        config.kmer_length,
        config.window_size,
    )?;
    let ref_time = ref_start.elapsed();

    if !config.quiet {
        let total_minimizers: usize = targets.iter().map(|t| t.minimizers.len()).sum();
        let total_bp: usize = targets.iter().map(|t| t.length).sum();
        eprintln!(
            "References: {} records, {}, {} minimizers",
            targets.len(),
            format_bp(total_bp),
            total_minimizers
        );
    }

    // Build set of target minimizers
    let reference_minimizers = if config.kmer_length <= 32 {
        let mut set = RapidHashSet::default();
        for target in &targets {
            if let MinimizerVec::U64(vec) = &target.minimizers {
                set.extend(vec.iter());
            }
        }
        MinimizerSet::U64(set)
    } else {
        let mut set = RapidHashSet::default();
        for target in &targets {
            if let MinimizerVec::U128(vec) = &target.minimizers {
                set.extend(vec.iter());
            }
        }
        MinimizerSet::U128(set)
    };

    let reference_minimizers = Arc::new(reference_minimizers);

    // Process reads file
    let reads_start = Instant::now();
    let (depth_map, total_reads, total_bp) = process_reads_file(
        &config.reads_path,
        reference_minimizers,
        config.kmer_length,
        config.window_size,
        config.threads,
        config.quiet,
    )?;
    let reads_time = reads_start.elapsed();

    // Coverage stats
    let analysis_start = Instant::now();
    let coverage_results = calculate_coverage_statistics(&targets, &depth_map);
    let analysis_time = analysis_start.elapsed();

    let total_time = start_time.elapsed();

    // Overall stats
    let total_minimizers: usize = coverage_results.iter().map(|r| r.total_minimizers).sum();
    let total_covered_minimizers: usize =
        coverage_results.iter().map(|r| r.covered_minimizers).sum();
    let overall_coverage_fraction = if total_minimizers > 0 {
        total_covered_minimizers as f64 / total_minimizers as f64
    } else {
        0.0
    };

    let overall_mean_depth = if total_minimizers > 0 {
        coverage_results
            .iter()
            .map(|r| r.mean_depth * r.total_minimizers as f64)
            .sum::<f64>()
            / total_minimizers as f64
    } else {
        0.0
    };

    // Calculate overall median depth from all median depths weighted by minimizer count
    let overall_median_depth = if total_minimizers > 0 {
        coverage_results
            .iter()
            .map(|r| r.median_depth * r.total_minimizers as f64)
            .sum::<f64>()
            / total_minimizers as f64
    } else {
        0.0
    };

    let reads_per_second = total_reads as f64 / reads_time.as_secs_f64();
    let bp_per_second = total_bp as f64 / reads_time.as_secs_f64();

    // Report
    let report = CoverageReport {
        version: format!("grate {}", version),
        reference_file: config.reference_path.to_string_lossy().to_string(),
        reads_file: if config.reads_path.to_string_lossy() == "-" {
            "stdin".to_string()
        } else {
            config.reads_path.to_string_lossy().to_string()
        },
        parameters: CoverageParameters {
            kmer_length: config.kmer_length,
            window_size: config.window_size,
            threads: config.threads,
        },
        targets: coverage_results,
        overall_stats: OverallStats {
            total_targets: targets.len(),
            total_minimizers,
            total_covered_minimizers,
            overall_coverage_fraction,
            overall_mean_depth,
            overall_median_depth,
            total_reads_processed: total_reads,
            total_bp_processed: total_bp,
        },
        timing: TimingStats {
            reference_processing_time: ref_time.as_secs_f64(),
            reads_processing_time: reads_time.as_secs_f64(),
            analysis_time: analysis_time.as_secs_f64(),
            total_time: total_time.as_secs_f64(),
            reads_per_second,
            bp_per_second,
        },
    };

    // Output results
    output_results(&report, config.output_path.as_ref(), config.output_format)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_bp() {
        assert_eq!(format_bp(500), "500bp");
        assert_eq!(format_bp(1500), "1.5Kbp");
        assert_eq!(format_bp(1500000), "1.5Mbp");
        assert_eq!(format_bp(1500000000), "1.5Gbp");
    }

    #[test]
    fn test_format_bp_per_sec() {
        assert_eq!(format_bp_per_sec(500.0), "500 bp/s");
        assert_eq!(format_bp_per_sec(1500.0), "1.5 Kbp/s");
        assert_eq!(format_bp_per_sec(1500000.0), "1.5 Mbp/s");
        assert_eq!(format_bp_per_sec(1500000000.0), "1.5 Gbp/s");
    }

    #[test]
    fn test_truncate_string() {
        assert_eq!(truncate_string("short", 10), "short");
        assert_eq!(truncate_string("verylongstring", 8), "veryl...");
    }
}

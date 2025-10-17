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

/// Sort order for results
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Original,     // Original order from reference file
    Alphabetical, // Alphabetical by target name
    Containment,  // Descending by containment (highest first)
}

/// Zero-cost (hopefully?) abstraction over u64 and u128 minimizer sets
#[derive(Debug, Clone)]
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
    pub minimizers: MinimizerSet,
    pub minimizer_positions: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageResult {
    pub target: String,
    pub length: usize,
    pub total_minimizers: usize,
    pub contained_minimizers: usize,
    pub containment: f64,
    pub median_abundance: f64,
    pub abundance_histogram: Vec<(u16, usize)>, // (abundance, count)
    pub containment_at_threshold: HashMap<usize, f64>, // threshold -> containment
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>, // Only used in multi-sample mode
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageParameters {
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
    pub abundance_thresholds: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverallStats {
    pub total_targets: usize,
    pub total_minimizers: usize,
    pub total_contained_minimizers: usize,
    pub overall_containment: f64,
    pub overall_median_abundance: f64,
    pub total_reads_processed: u64,
    pub total_bp_processed: u64,
    pub overall_containment_at_threshold: HashMap<usize, f64>, // threshold -> overall containment
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

/// Results for a single sample in multi-sample mode
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleResults {
    pub sample_name: String,
    pub reads_file: String,
    pub targets: Vec<CoverageResult>,
    pub overall_stats: OverallStats,
    pub timing: TimingStats,
}

/// Report containing results for one or more samples
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Report {
    pub version: String,
    pub targets_file: String,
    pub parameters: CoverageParameters,
    pub samples: Vec<SampleResults>,
    pub total_timing: TimingStats,
}

/// Output format options
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Table,
    Csv,
    Json,
}

pub struct CoverageConfig {
    pub targets_path: PathBuf,
    pub reads_paths: Vec<PathBuf>,
    pub sample_names: Vec<String>,
    pub kmer_length: u8,
    pub window_size: u8,
    pub threads: usize,
    pub output_path: Option<PathBuf>,
    pub quiet: bool,
    pub output_format: OutputFormat,
    pub abundance_thresholds: Vec<usize>,
    pub discriminatory: bool,
    pub limit_bp: Option<u64>,
    pub sort_order: SortOrder,
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

/// Processor for collecting target sequence records with minimizers
#[derive(Clone)]
struct TargetsProcessor {
    kmer_length: u8,
    window_size: u8,
    hasher: KmerHasher,
    buffers: Buffers,
    positions: Vec<usize>,
    targets: Arc<Mutex<Vec<TargetInfo>>>,

    // Progress tracking
    local_stats: ProcessingStats,
    global_stats: Arc<Mutex<ProcessingStats>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: std::time::Instant,
}

impl TargetsProcessor {
    fn new(
        kmer_length: u8,
        window_size: u8,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: std::time::Instant,
    ) -> Self {
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
            local_stats: ProcessingStats::default(),
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            spinner,
            start_time,
        }
    }

    fn update_spinner(&self) {
        if let Some(ref spinner) = self.spinner {
            let stats = self.global_stats.lock();
            let elapsed = self.start_time.elapsed();
            let seqs_per_sec = stats.total_seqs as f64 / elapsed.as_secs_f64();
            let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();

            spinner.lock().set_message(format!(
                "Processing targets: {} seqs ({}). {:.0} seqs/s ({})",
                stats.total_seqs,
                format_bp(stats.total_bp as usize),
                seqs_per_sec,
                format_bp_per_sec(bp_per_sec)
            ));
        }
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for TargetsProcessor {
    fn process_record(&mut self, record: Rf) -> paraseq::parallel::Result<()> {
        let sequence = record.seq();
        let target_name = String::from_utf8_lossy(record.id()).to_string();

        self.local_stats.total_seqs += 1;
        self.local_stats.total_bp += sequence.len() as u64;

        fill_minimizers_with_positions(
            &sequence,
            &self.hasher,
            self.kmer_length,
            self.window_size,
            &mut self.buffers,
            &mut self.positions,
        );

        // Build unique minimizer set for this target
        let minimizers = match &self.buffers.minimizers {
            MinimizerVec::U64(vec) => {
                let set: RapidHashSet<u64> = vec.iter().copied().collect();
                MinimizerSet::U64(set)
            }
            MinimizerVec::U128(vec) => {
                let set: RapidHashSet<u128> = vec.iter().copied().collect();
                MinimizerSet::U128(set)
            }
        };

        let minimizer_positions = self.positions.clone();

        self.targets.lock().push(TargetInfo {
            name: target_name,
            length: sequence.len(),
            minimizers,
            minimizer_positions,
        });

        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            // Update spinner every 0.1 Gbp
            let current_progress = stats.total_bp / 100_000_000; // 0.1 Gbp increments
            if current_progress > stats.last_reported {
                drop(stats); // Release lock before updating spinner
                self.update_spinner();
                self.global_stats.lock().last_reported = current_progress;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

fn process_targets_file(
    targets_path: &Path,
    kmer_length: u8,
    window_size: u8,
    quiet: bool,
) -> Result<Vec<TargetInfo>> {
    let in_path = if targets_path.to_string_lossy() == "-" {
        None
    } else {
        Some(targets_path)
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
        pb.set_message("Processing targets: 0 seqs (0bp)");
        Some(Arc::new(Mutex::new(pb)))
    } else {
        None
    };

    let start_time = std::time::Instant::now();
    let mut processor =
        TargetsProcessor::new(kmer_length, window_size, spinner.clone(), start_time);

    // Single thread to preserve order
    reader.process_parallel(&mut processor, 1)?;

    // Finish spinner and clear
    if let Some(ref pb) = spinner {
        pb.lock().finish_and_clear();
    }

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
    targets_minimizers: Arc<MinimizerSet>,

    // Local buffers
    buffers: Buffers,
    local_stats: ProcessingStats,
    local_counts_u64: Option<HashMap<u64, u16>>,
    local_counts_u128: Option<HashMap<u128, u16>>,

    // Global state
    global_stats: Arc<Mutex<ProcessingStats>>,
    global_counts_u64: Arc<Mutex<Option<HashMap<u64, u16>>>>,
    global_counts_u128: Arc<Mutex<Option<HashMap<u128, u16>>>>,
    spinner: Option<Arc<Mutex<ProgressBar>>>,
    start_time: Instant,
    limit_bp: Option<u64>,
}

impl ReadsProcessor {
    fn new(
        kmer_length: u8,
        window_size: u8,
        targets_minimizers: Arc<MinimizerSet>,
        spinner: Option<Arc<Mutex<ProgressBar>>>,
        start_time: Instant,
        limit_bp: Option<u64>,
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
            targets_minimizers,
            buffers,
            local_stats: ProcessingStats::default(),
            local_counts_u64,
            local_counts_u128,
            global_stats: Arc::new(Mutex::new(ProcessingStats::default())),
            global_counts_u64,
            global_counts_u128,
            spinner,
            start_time,
            limit_bp,
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
        // Check if we've reached sample limit
        if let Some(limit) = self.limit_bp {
            let global_bp = self.global_stats.lock().total_bp;
            if global_bp >= limit {
                // Stop processing - return error to halt pipeline
                return Err(paraseq::parallel::ProcessError::IoError(
                    std::io::Error::new(std::io::ErrorKind::Interrupted, "Sample limit reached"),
                ));
            }
        }

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

        // Count minimizers present in targets
        match (&self.buffers.minimizers, &*self.targets_minimizers) {
            (MinimizerVec::U64(vec), MinimizerSet::U64(targets_set)) => {
                let local_counts = self.local_counts_u64.as_mut().unwrap();
                for &minimizer in vec {
                    if targets_set.contains(&minimizer) {
                        local_counts
                            .entry(minimizer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
                    }
                }
            }
            (MinimizerVec::U128(vec), MinimizerSet::U128(targets_set)) => {
                let local_counts = self.local_counts_u128.as_mut().unwrap();
                for &minimizer in vec {
                    if targets_set.contains(&minimizer) {
                        local_counts
                            .entry(minimizer)
                            .and_modify(|e| *e = e.saturating_add(1))
                            .or_insert(1);
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
                global_map
                    .entry(minimizer)
                    .and_modify(|e| *e = e.saturating_add(count))
                    .or_insert(count);
            }
            local.clear();
        } else {
            let mut global = self.global_counts_u128.lock();
            let global_map = global.as_mut().unwrap();
            let local = self.local_counts_u128.as_mut().unwrap();
            for (&minimizer, &count) in local.iter() {
                global_map
                    .entry(minimizer)
                    .and_modify(|e| *e = e.saturating_add(count))
                    .or_insert(count);
            }
            local.clear();
        }

        // Update global stats
        {
            let mut stats = self.global_stats.lock();
            stats.total_seqs += self.local_stats.total_seqs;
            stats.total_bp += self.local_stats.total_bp;

            // Update spinner every 0.1 Gbp
            let current_progress = stats.total_bp / 100_000_000; // 0.1 Gbp increments
            if current_progress > stats.last_reported {
                drop(stats); // Release lock before updating spinner
                self.update_spinner();
                self.global_stats.lock().last_reported = current_progress;
            }

            self.local_stats = ProcessingStats::default();
        }

        Ok(())
    }
}

/// Enum to abstract over u64 and u128 abundance maps
enum AbundanceMap {
    U64(HashMap<u64, u16>),
    U128(HashMap<u128, u16>),
}

fn process_reads_file(
    reads_path: &Path,
    targets_minimizers: Arc<MinimizerSet>,
    kmer_length: u8,
    window_size: u8,
    threads: usize,
    quiet: bool,
    limit_bp: Option<u64>,
) -> Result<(AbundanceMap, u64, u64)> {
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

    let total_target_minimizers = targets_minimizers.len();

    let start_time = Instant::now();
    let mut processor = ReadsProcessor::new(
        kmer_length,
        window_size,
        targets_minimizers,
        spinner.clone(),
        start_time,
        limit_bp,
    );

    // Process reads - may terminate early if sample limit reached
    let process_result = reader.process_parallel(&mut processor, threads);

    // Check if we stopped due to sampling
    let stopped_early = if let Err(ref e) = process_result {
        e.to_string().contains("Sample limit reached")
    } else {
        false
    };

    // If it's not a sample stop, propagate the error
    if !stopped_early {
        process_result?;
    }

    // Finish spinner
    if let Some(ref pb) = spinner {
        pb.lock().finish_with_message("");
    }

    let stats = processor.global_stats.lock().clone();
    let abundance_map = if kmer_length <= 32 {
        let map = Arc::try_unwrap(processor.global_counts_u64)
            .unwrap()
            .into_inner()
            .unwrap();
        AbundanceMap::U64(map)
    } else {
        let map = Arc::try_unwrap(processor.global_counts_u128)
            .unwrap()
            .into_inner()
            .unwrap();
        AbundanceMap::U128(map)
    };

    if !quiet {
        let elapsed = start_time.elapsed();
        let unique_minimizers = match &abundance_map {
            AbundanceMap::U64(m) => m.len(),
            AbundanceMap::U128(m) => m.len(),
        };
        let bp_per_sec = stats.total_bp as f64 / elapsed.as_secs_f64();
        eprintln!(
            "Reads: {} records ({}), found {} of {} target minimizers ({})",
            stats.total_seqs,
            format_bp(stats.total_bp as usize),
            unique_minimizers,
            total_target_minimizers,
            format_bp_per_sec(bp_per_sec)
        );
    }

    Ok((abundance_map, stats.total_seqs, stats.total_bp))
}

fn calculate_containment_statistics(
    targets: &[TargetInfo],
    abundance_map: &AbundanceMap,
    abundance_thresholds: &[usize],
    sample_name: Option<String>,
) -> Vec<CoverageResult> {
    targets
        .iter()
        .map(|target| {
            let total_minimizers = target.minimizers.len();
            let mut abundances: Vec<u16> = Vec::new();
            let mut contained_count = 0;

            // Collect abundances for all unique minimizers in this target
            match (&target.minimizers, abundance_map) {
                (MinimizerSet::U64(set), AbundanceMap::U64(map)) => {
                    for &minimizer in set {
                        let abundance = map.get(&minimizer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                (MinimizerSet::U128(set), AbundanceMap::U128(map)) => {
                    for &minimizer in set {
                        let abundance = map.get(&minimizer).copied().unwrap_or(0);
                        abundances.push(abundance);
                        if abundance > 0 {
                            contained_count += 1;
                        }
                    }
                }
                _ => panic!("Mismatch between MinimizerSet and AbundanceMap types"),
            }

            let containment = if total_minimizers > 0 {
                contained_count as f64 / total_minimizers as f64
            } else {
                0.0
            };

            // Ignore zero-abundance minimizers for median calc
            let non_zero_abundances: Vec<u16> =
                abundances.iter().copied().filter(|&a| a > 0).collect();

            // Calculate median w/o zeros
            let mut sorted_abundances = non_zero_abundances;
            sorted_abundances.sort_unstable();
            let median_abundance = if sorted_abundances.is_empty() {
                0.0
            } else if sorted_abundances.len() % 2 == 0 {
                let mid = sorted_abundances.len() / 2;
                (sorted_abundances[mid - 1] + sorted_abundances[mid]) as f64 / 2.0
            } else {
                sorted_abundances[sorted_abundances.len() / 2] as f64
            };

            // Abundance hist
            let mut abundance_counts: HashMap<u16, usize> = HashMap::new();
            for abundance in &abundances {
                *abundance_counts.entry(*abundance).or_insert(0) += 1;
            }
            let mut abundance_histogram: Vec<(u16, usize)> = abundance_counts.into_iter().collect();
            abundance_histogram.sort_by_key(|(abundance, _)| *abundance);

            // Calculate containment at threshold
            let mut containment_at_threshold = HashMap::new();
            for &threshold in abundance_thresholds {
                let count_at_threshold = abundances
                    .iter()
                    .filter(|&&a| a >= threshold as u16)
                    .count();
                let containment_value = if total_minimizers > 0 {
                    count_at_threshold as f64 / total_minimizers as f64
                } else {
                    0.0
                };
                containment_at_threshold.insert(threshold, containment_value);
            }

            CoverageResult {
                target: target.name.clone(),
                length: target.length,
                total_minimizers,
                contained_minimizers: contained_count,
                containment,
                median_abundance,
                abundance_histogram,
                containment_at_threshold,
                sample_name: sample_name.clone(),
            }
        })
        .collect()
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
        let mut options = format!(
            "k={}, w={}, threads={}",
            config.kmer_length, config.window_size, config.threads
        );

        if config.reads_paths.len() > 1 {
            options.push_str(&format!(", samples={}", config.reads_paths.len()));
        }

        if !config.abundance_thresholds.is_empty() {
            options.push_str(&format!(
                ", abundance-thresholds={}",
                config
                    .abundance_thresholds
                    .iter()
                    .map(|t| t.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }

        if config.discriminatory {
            options.push_str(", discriminatory");
        }

        if let Some(limit) = config.limit_bp {
            options.push_str(&format!(", limit={}", format_bp(limit as usize)));
        }

        options.push_str(&format!(
            ", format={}",
            match config.output_format {
                OutputFormat::Table => "table",
                OutputFormat::Csv => "csv",
                OutputFormat::Json => "json",
            }
        ));

        eprintln!("Grate v{}; mode: cov; options: {}", version, options);
    }

    // Process targets file
    let targets_start = Instant::now();
    let mut targets = process_targets_file(
        &config.targets_path,
        config.kmer_length,
        config.window_size,
        config.quiet,
    )?;
    let targets_time = targets_start.elapsed();

    // Count minimizers shared by multiple targets
    let mut minimizer_target_counts: HashMap<u64, usize> = HashMap::new();
    let mut minimizer_target_counts_u128: HashMap<u128, usize> = HashMap::new();

    for target in &targets {
        match &target.minimizers {
            MinimizerSet::U64(set) => {
                for &minimizer in set {
                    *minimizer_target_counts.entry(minimizer).or_insert(0) += 1;
                }
            }
            MinimizerSet::U128(set) => {
                for &minimizer in set {
                    *minimizer_target_counts_u128.entry(minimizer).or_insert(0) += 1;
                }
            }
        }
    }

    let shared_minimizers = if config.kmer_length <= 32 {
        minimizer_target_counts
            .values()
            .filter(|&&count| count > 1)
            .count()
    } else {
        minimizer_target_counts_u128
            .values()
            .filter(|&&count| count > 1)
            .count()
    };

    let unique_across_all = if config.kmer_length <= 32 {
        minimizer_target_counts.len()
    } else {
        minimizer_target_counts_u128.len()
    };

    // Apply discriminatory filtering if enabled
    if config.discriminatory {
        for target in &mut targets {
            match &mut target.minimizers {
                MinimizerSet::U64(set) => {
                    set.retain(|minimizer| {
                        minimizer_target_counts
                            .get(minimizer)
                            .map_or(true, |&count| count == 1)
                    });
                }
                MinimizerSet::U128(set) => {
                    set.retain(|minimizer| {
                        minimizer_target_counts_u128
                            .get(minimizer)
                            .map_or(true, |&count| count == 1)
                    });
                }
            }
        }
    }

    if !config.quiet {
        let total_unique_minimizers: usize = targets.iter().map(|t| t.minimizers.len()).sum();
        let total_bp: usize = targets.iter().map(|t| t.length).sum();

        if config.discriminatory {
            eprintln!(
                "Targets: {} records ({}), {} discriminatory minimizers ({} shared minimizers dropped)",
                targets.len(),
                format_bp(total_bp),
                total_unique_minimizers,
                shared_minimizers
            );
        } else {
            let shared_pct = if unique_across_all > 0 {
                shared_minimizers as f64 / unique_across_all as f64 * 100.0
            } else {
                0.0
            };
            eprintln!(
                "Targets: {} records ({}), {} minimizers, of which {} ({:.1}%) shared by multiple targets",
                targets.len(),
                format_bp(total_bp),
                total_unique_minimizers,
                shared_minimizers,
                shared_pct
            );
        }
    }

    // Build set of all unique minimizers across targets
    let targets_minimizers = if config.kmer_length <= 32 {
        let mut set = RapidHashSet::default();
        for target in &targets {
            if let MinimizerSet::U64(target_set) = &target.minimizers {
                set.extend(target_set.iter());
            }
        }
        MinimizerSet::U64(set)
    } else {
        let mut set = RapidHashSet::default();
        for target in &targets {
            if let MinimizerSet::U128(target_set) = &target.minimizers {
                set.extend(target_set.iter());
            }
        }
        MinimizerSet::U128(set)
    };

    let targets_minimizers = Arc::new(targets_minimizers);

    // Process each sample
    let mut sample_results = Vec::new();
    let mut total_reads_time = 0.0;
    let mut total_analysis_time = 0.0;

    for (idx, (reads_path, sample_name)) in config
        .reads_paths
        .iter()
        .zip(&config.sample_names)
        .enumerate()
    {
        if !config.quiet {
            eprintln!(
                "Processing sample '{}' ({}/{}):",
                sample_name,
                idx + 1,
                config.reads_paths.len()
            );
        }

        // Process reads for this sample
        let reads_start = Instant::now();
        let (abundance_map, total_reads, total_bp) = process_reads_file(
            reads_path,
            Arc::clone(&targets_minimizers),
            config.kmer_length,
            config.window_size,
            config.threads,
            config.quiet,
            config.limit_bp,
        )?;
        let reads_time = reads_start.elapsed();
        total_reads_time += reads_time.as_secs_f64();

        // Calculate coverage statistics for this sample
        let analysis_start = Instant::now();
        let coverage_results = calculate_containment_statistics(
            &targets,
            &abundance_map,
            &config.abundance_thresholds,
            Some(sample_name.clone()),
        );
        let analysis_time = analysis_start.elapsed();
        total_analysis_time += analysis_time.as_secs_f64();

        // Overall stats for this sample
        let total_minimizers: usize = coverage_results.iter().map(|r| r.total_minimizers).sum();
        let total_contained_minimizers: usize = coverage_results
            .iter()
            .map(|r| r.contained_minimizers)
            .sum();
        let overall_containment = if total_minimizers > 0 {
            total_contained_minimizers as f64 / total_minimizers as f64
        } else {
            0.0
        };

        let overall_median_abundance = if total_minimizers > 0 {
            coverage_results
                .iter()
                .map(|r| r.median_abundance * r.total_minimizers as f64)
                .sum::<f64>()
                / total_minimizers as f64
        } else {
            0.0
        };

        let reads_per_second = total_reads as f64 / reads_time.as_secs_f64();
        let bp_per_second = total_bp as f64 / reads_time.as_secs_f64();

        // Calculate overall containment at each threshold
        let mut overall_containment_at_threshold = HashMap::new();
        for &threshold in &config.abundance_thresholds {
            let total_at_threshold: usize = coverage_results
                .iter()
                .map(|r| {
                    let containment = r.containment_at_threshold.get(&threshold).unwrap_or(&0.0);
                    (containment * r.total_minimizers as f64) as usize
                })
                .sum();
            let overall_containment_value = if total_minimizers > 0 {
                total_at_threshold as f64 / total_minimizers as f64
            } else {
                0.0
            };
            overall_containment_at_threshold.insert(threshold, overall_containment_value);
        }

        sample_results.push(SampleResults {
            sample_name: sample_name.clone(),
            reads_file: if reads_path.to_string_lossy() == "-" {
                "stdin".to_string()
            } else {
                reads_path.to_string_lossy().to_string()
            },
            targets: coverage_results,
            overall_stats: OverallStats {
                total_targets: targets.len(),
                total_minimizers,
                total_contained_minimizers,
                overall_containment,
                overall_median_abundance,
                total_reads_processed: total_reads,
                total_bp_processed: total_bp,
                overall_containment_at_threshold,
            },
            timing: TimingStats {
                reference_processing_time: 0.0, // Not per-sample
                reads_processing_time: reads_time.as_secs_f64(),
                analysis_time: analysis_time.as_secs_f64(),
                total_time: reads_time.as_secs_f64() + analysis_time.as_secs_f64(),
                reads_per_second,
                bp_per_second,
            },
        });
    }

    let total_time = start_time.elapsed();

    // Create report
    let report = Report {
        version: format!("grate {}", version),
        targets_file: config.targets_path.to_string_lossy().to_string(),
        parameters: CoverageParameters {
            kmer_length: config.kmer_length,
            window_size: config.window_size,
            threads: config.threads,
            abundance_thresholds: config.abundance_thresholds.clone(),
        },
        samples: sample_results,
        total_timing: TimingStats {
            reference_processing_time: targets_time.as_secs_f64(),
            reads_processing_time: total_reads_time,
            analysis_time: total_analysis_time,
            total_time: total_time.as_secs_f64(),
            reads_per_second: 0.0, // Not meaningful across samples
            bp_per_second: 0.0,    // Not meaningful across samples
        },
    };

    // Output results
    output_results(
        &report,
        config.output_path.as_ref(),
        config.output_format,
        config.sort_order,
    )?;

    Ok(())
}

/// Sort coverage results based on the specified sort order
fn sort_results(results: &mut [CoverageResult], sort_order: SortOrder) {
    match sort_order {
        SortOrder::Original => {
            // No sorting needed - keep original order
        }
        SortOrder::Alphabetical => {
            results.sort_by(|a, b| a.target.cmp(&b.target));
        }
        SortOrder::Containment => {
            // Sort by containment descending (highest first)
            results.sort_by(|a, b| {
                b.containment
                    .partial_cmp(&a.containment)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
    }
}

fn output_results(
    report: &Report,
    output_path: Option<&PathBuf>,
    output_format: OutputFormat,
    sort_order: SortOrder,
) -> Result<()> {
    // Create a sorted copy of the report if needed
    let mut sorted_report = report.clone();
    if sort_order != SortOrder::Original {
        for sample in &mut sorted_report.samples {
            sort_results(&mut sample.targets, sort_order);
        }
    }

    let writer: Box<dyn Write> = if let Some(path) = output_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    let mut writer = writer;

    match output_format {
        OutputFormat::Json => {
            // Output JSON
            serde_json::to_writer_pretty(&mut writer, &sorted_report)?;
            writeln!(writer)?;
        }
        OutputFormat::Csv => {
            output_csv(&mut writer, &sorted_report)?;
        }
        OutputFormat::Table => {
            output_table(&mut writer, &sorted_report)?;
        }
    }

    Ok(())
}

fn output_csv(writer: &mut dyn Write, report: &Report) -> Result<()> {
    let mut thresholds = report.parameters.abundance_thresholds.clone();
    thresholds.sort_unstable();

    // Build header with sample column first
    let mut header = "sample,target,containment".to_string();
    for threshold in &thresholds {
        header.push_str(&format!(",containment{}", threshold));
    }
    header.push_str(",length_bp,total_minimizers,contained_minimizers,median_abundance");
    writeln!(writer, "{}", header)?;

    // Output data rows for all samples
    for sample in &report.samples {
        for result in &sample.targets {
            let mut row = format!(
                "{},{},{:.5}",
                sample.sample_name, result.target, result.containment
            );
            for threshold in &thresholds {
                let containment = result
                    .containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                row.push_str(&format!(",{:.5}", containment));
            }
            row.push_str(&format!(
                ",{},{},{},{:.0}",
                result.length,
                result.total_minimizers,
                result.contained_minimizers,
                result.median_abundance,
            ));
            writeln!(writer, "{}", row)?;
        }
    }

    Ok(())
}

fn output_table(writer: &mut dyn Write, report: &Report) -> Result<()> {
    let mut thresholds = report.parameters.abundance_thresholds.clone();
    thresholds.sort_unstable();

    // Build header with sample column first
    let mut header = format!(
        "{:<20} | {:<30} | {:>12} ",
        "sample", "target", "containment"
    );
    for threshold in &thresholds {
        header.push_str(&format!("| {:>14} ", format!("containment{}", threshold)));
    }
    header.push_str(&format!(
        "| {:>10} | {:>12} | {:>17}",
        "length", "minimizers", "median_abundance"
    ));

    let separator = "-".repeat(header.len());

    // Output header
    writeln!(writer)?;
    writeln!(writer, "{}", header)?;
    writeln!(writer, "{}", separator)?;

    // Output data rows for all samples
    for sample in &report.samples {
        for result in &sample.targets {
            let mut row = format!(
                "{:<20} | {:<30} | {:>11.1}% ",
                truncate_string(&sample.sample_name, 20),
                truncate_string(&result.target, 30),
                result.containment * 100.0
            );
            for threshold in &thresholds {
                let containment = result
                    .containment_at_threshold
                    .get(threshold)
                    .unwrap_or(&0.0);
                row.push_str(&format!("| {:>13.1}% ", containment * 100.0));
            }
            row.push_str(&format!(
                "| {:>10} | {:>12} | {:>17.0}",
                format_bp(result.length),
                result.total_minimizers,
                result.median_abundance,
            ));
            writeln!(writer, "{}", row)?;
        }
    }

    // Overall statistics per sample
    writeln!(writer)?;
    writeln!(writer, "Overall statistics per sample:")?;
    for sample in &report.samples {
        let mut overall_msg = format!(
            "  {}: containment={:.1}%",
            sample.sample_name,
            sample.overall_stats.overall_containment * 100.0
        );
        for threshold in &thresholds {
            let containment = sample
                .overall_stats
                .overall_containment_at_threshold
                .get(threshold)
                .unwrap_or(&0.0);
            overall_msg.push_str(&format!(
                ", containment{}={:.1}%",
                threshold,
                containment * 100.0
            ));
        }
        overall_msg.push_str(&format!(
            ", median_abundance={:.0}",
            sample.overall_stats.overall_median_abundance,
        ));
        writeln!(writer, "{}", overall_msg)?;
    }

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

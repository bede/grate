//! # Deacon
//!
//! A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format,
//! built for efficient host depletion (*deacon*-tamination).
//!
//! This crate provides both a library and a binary for filtering nucleotide sequences.
//!
#![doc = include_str!("../README.md")]

// Re-export public functionality
pub mod filter;
pub mod index;
pub mod minimizers;

// Re-export the important structures and functions for library users
pub use filter::{FilterSummary, run as run_filter};
pub use index::{
    IndexHeader, build as build_index, diff as diff_index, info as index_info, union as union_index,
};
pub use minimizers::{
    DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE, compute_minimizer_hashes, fill_minimizer_hashes,
};

use anyhow::Result;
use rustc_hash::FxHashSet;
use std::path::{Path, PathBuf};

/// Configuration for filtering operations
pub struct FilterConfig {
    /// Minimizer index file path
    pub minimizers_path: PathBuf,

    /// Path to input fastx file (or - for stdin)
    pub input_path: String,

    /// Path to optional second paired fastx file (or - for interleaved stdin)
    pub input2_path: Option<String>,

    /// Path to output fastx file (or - for stdout; detects .gz and .zst)
    pub output_path: String,

    /// Path to optional second output fastx file for paired reads (detects .gz and .zst)
    pub output2_path: Option<String>,

    /// Minimum number of minimizer matches per query sequence (pair)
    pub min_matches: usize,

    /// Consider only the first N nucleotides per sequence (0 = entire sequence)
    pub prefix_length: usize,

    /// Path to JSON summary file
    pub summary_path: Option<PathBuf>,

    /// Invert filtering (keep sequences WITH matches rather than those WITHOUT)
    pub invert: bool,

    /// Replace sequence headers with sequential numbers (1, 2, 3...)
    pub rename: bool,

    /// Number of execution threads (0 = auto)
    pub threads: usize,
}

impl FilterConfig {
    /// Create a new filtering configuration with the specified minimizers path
    pub fn new<P: AsRef<Path>>(minimizers_path: P) -> Self {
        Self {
            minimizers_path: minimizers_path.as_ref().to_path_buf(),
            input_path: "-".to_string(),
            input2_path: None,
            output_path: "-".to_string(),
            output2_path: None,
            min_matches: 2,
            prefix_length: 0,
            summary_path: None,
            invert: false,
            rename: false,
            threads: 0, // Use all available threads by default
        }
    }

    /// Set the input path
    pub fn with_input<S: Into<String>>(mut self, input_path: S) -> Self {
        self.input_path = input_path.into();
        self
    }

    /// Set the second input path for paired reads
    pub fn with_input2<S: Into<String>>(mut self, input2_path: S) -> Self {
        self.input2_path = Some(input2_path.into());
        self
    }

    /// Set the output path
    pub fn with_output<S: Into<String>>(mut self, output_path: S) -> Self {
        self.output_path = output_path.into();
        self
    }

    /// Set the second output path for paired reads
    pub fn with_output2<S: Into<String>>(mut self, output2_path: S) -> Self {
        self.output2_path = Some(output2_path.into());
        self
    }

    /// Set the minimum matches
    pub fn with_min_matches(mut self, min_matches: usize) -> Self {
        self.min_matches = min_matches;
        self
    }

    /// Set the prefix length
    pub fn with_prefix_length(mut self, prefix_length: usize) -> Self {
        self.prefix_length = prefix_length;
        self
    }

    /// Set the summary path
    pub fn with_summary<P: AsRef<Path>>(mut self, summary_path: P) -> Self {
        self.summary_path = Some(summary_path.as_ref().to_path_buf());
        self
    }

    /// Set invert filtering
    pub fn with_invert(mut self, invert: bool) -> Self {
        self.invert = invert;
        self
    }

    /// Set rename option
    pub fn with_rename(mut self, rename: bool) -> Self {
        self.rename = rename;
        self
    }

    /// Set the num threads
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Execute the filtering operation with this configuration
    pub fn execute(&self) -> Result<()> {
        filter::run(
            &self.minimizers_path,
            &self.input_path,
            self.input2_path.as_deref(),
            &self.output_path,
            self.output2_path.as_deref(),
            self.min_matches,
            self.prefix_length,
            self.summary_path.as_ref(),
            self.invert,
            self.rename,
            self.threads,
        )
    }
}

pub struct IndexConfig {
    /// Path to input fastx file
    pub input_path: PathBuf,

    /// K-mer length used for indexing
    pub kmer_length: usize,

    /// Minimizer window size used for indexing
    pub window_size: usize,

    /// Path to output file (None for stdout)
    pub output_path: Option<PathBuf>,

    /// Hash table pre-allocation capacity in millions
    pub capacity_millions: usize,

    /// Number of execution threads (0 = auto)
    pub threads: usize,
}

impl IndexConfig {
    /// Create a new index configuration with the specified input path
    pub fn new<P: AsRef<Path>>(input_path: P) -> Self {
        Self {
            input_path: input_path.as_ref().to_path_buf(),
            kmer_length: DEFAULT_KMER_LENGTH,
            window_size: DEFAULT_WINDOW_SIZE,
            output_path: None,
            capacity_millions: 500, // Default 500M capacity
            threads: 0,             // Use all available threads by default
        }
    }

    /// Set k-mer length
    pub fn with_kmer_length(mut self, kmer_length: usize) -> Self {
        self.kmer_length = kmer_length;
        self
    }

    /// Set window size
    pub fn with_window_size(mut self, window_size: usize) -> Self {
        self.window_size = window_size;
        self
    }

    /// Set output path
    pub fn with_output<P: AsRef<Path>>(mut self, output_path: P) -> Self {
        self.output_path = Some(output_path.as_ref().to_path_buf());
        self
    }

    /// Set hash table capacity in millions
    pub fn with_capacity_millions(mut self, capacity_millions: usize) -> Self {
        self.capacity_millions = capacity_millions;
        self
    }

    /// Set threads
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Execute index build with this configuration
    pub fn execute(&self) -> Result<()> {
        build_index(
            &self.input_path,
            self.kmer_length,
            self.window_size,
            self.output_path.clone(),
            self.capacity_millions,
            self.threads,
        )
    }
}

pub fn load_minimizers<P: AsRef<Path>>(path: P) -> Result<(FxHashSet<u64>, index::IndexHeader)> {
    index::load_minimizer_hashes(&path)
}

pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &index::IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    index::write_minimizers(minimizers, header, output_path)
}

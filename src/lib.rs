//! # Grate
//!
//! A fast minimizer-based coverage analysis tool for genomic sequences.
//!
//! Grate analyzes the coverage of minimizers from reference sequences in read datasets,
//! providing detailed statistics on coverage depth and breadth per target sequence.
//!

pub mod coverage;
pub mod minimizers;

// Re-export the main functionality
pub use coverage::{
    run_coverage_analysis, CoverageConfig, CoverageParameters, CoverageReport, CoverageResult,
    OutputFormat, OverallStats, TargetInfo, TimingStats,
};

pub use minimizers::{
    decode_u128, decode_u64, fill_minimizers, fill_minimizers_with_positions, Buffers, KmerHasher,
    MinimizerVec, DEFAULT_KMER_LENGTH, DEFAULT_WINDOW_SIZE,
};

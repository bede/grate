use packed_seq::AsciiSeq;
use rustc_hash::FxHashSet;
use xxhash_rust::xxh3;

pub const DEFAULT_KMER_LENGTH: u8 = 31;
pub const DEFAULT_WINDOW_SIZE: u16 = 15;

/// Check if nucleotide is valid ACGT (case insensitive)
#[inline]
fn is_valid_acgt(nucleotide: u8) -> bool {
    matches!(
        nucleotide,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
    )
}

/// Check if k-mer contains only ACGT nucleotides
#[inline]
fn kmer_contains_only_acgt(kmer: &[u8]) -> bool {
    kmer.iter().all(|&b| is_valid_acgt(b))
}

/// Canonicalise IUPAC ambiguous nucleotides to ACGT
#[inline]
fn canonicalise_nucleotide(nucleotide: u8) -> u8 {
    match nucleotide {
        b'A' | b'a' => b'A',
        b'C' | b'c' => b'C',
        b'G' | b'g' => b'G',
        b'T' | b't' => b'T',
        b'R' | b'r' => b'G',
        b'Y' | b'y' => b'C',
        b'S' | b's' => b'G',
        b'W' | b'w' => b'A',
        b'K' | b'k' => b'G',
        b'M' | b'm' => b'C',
        b'B' | b'b' => b'C',
        b'D' | b'd' => b'G',
        b'H' | b'h' => b'C',
        b'V' | b'v' => b'G',
        b'N' | b'n' => b'C',
        _ => b'C',
    }
}

/// Canonicalise a sequence
fn canonicalise_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&nucleotide| canonicalise_nucleotide(nucleotide))
        .collect()
}

/// Returns vector of all minimizer hashes for a sequence
pub fn compute_minimizer_hashes(
    seq: &[u8],
    kmer_length: u8,
    window_size: u16,
    information_threshold: Option<f32>,
) -> Vec<u64> {
    let mut hashes = Vec::new();
    fill_minimizer_hashes(
        seq,
        kmer_length,
        window_size,
        &mut hashes,
        information_threshold,
    );
    hashes
}

/// Calculate linguistic complexity for a k-mer using vocabulary richness approach
/// Returns value between 0.0 (maximally repetitive) and 1.0 (maximally complex)
fn calculate_linguistic_complexity(kmer: &[u8]) -> f32 {
    let n = kmer.len();

    // K-mers less than 4 bases long always pass through (never filtered)
    if n < 4 {
        return 1.0;
    }

    // Always use max word length of 4, regardless of k-mer length
    let window_size = 4;
    let mut complexity = 1.0;

    // Calculate vocabulary usage for each word length from 1 to window_size
    for word_len in 1..=window_size {
        if word_len > n {
            break;
        }

        // Count unique subwords of this length
        let mut vocab = FxHashSet::default();
        for i in 0..=(n - word_len) {
            vocab.insert(&kmer[i..i + word_len]);
        }

        let observed_vocab = vocab.len();
        let max_possible_vocab = 4_usize.pow(word_len as u32).min(n - word_len + 1);
        let vocab_usage = observed_vocab as f32 / max_possible_vocab as f32;

        complexity *= vocab_usage;
    }

    complexity
}

/// Fill a vector with minimizer hashes, skipping k-mers with non-ACGT bases
/// and optionally filtering by linguistic complexity
pub fn fill_minimizer_hashes(
    seq: &[u8],
    kmer_length: u8,
    window_size: u16,
    hashes: &mut Vec<u64>,
    information_threshold: Option<f32>,
) {
    hashes.clear();

    // Skip if sequence is too short
    if seq.len() < kmer_length as usize {
        return;
    }

    let canonical_seq = canonicalise_sequence(seq);

    // Get minimizer positions using simd-minimizers
    let mut positions = Vec::new();
    simd_minimizers::canonical_minimizer_positions(
        AsciiSeq(&canonical_seq),
        kmer_length as usize,
        window_size as usize,
        &mut positions,
    );

    // Filter positions to only include k-mers with ACGT bases and sufficient complexity
    let valid_positions: Vec<u32> = positions
        .into_iter()
        .filter(|&pos| {
            let pos_usize = pos as usize;
            let kmer = &seq[pos_usize..pos_usize + kmer_length as usize];

            // First check ACGT constraint
            if !kmer_contains_only_acgt(kmer) {
                return false;
            }

            // Then check complexity constraint if threshold is specified
            if let Some(threshold) = information_threshold {
                let complexity = calculate_linguistic_complexity(kmer);
                complexity >= threshold
            } else {
                true
            }
        })
        .collect();

    hashes.extend(
        simd_minimizers::iter_canonical_minimizer_values(
            AsciiSeq(&canonical_seq),
            kmer_length as usize,
            &valid_positions,
        )
        .map(|kmer| xxh3::xxh3_64(&kmer.to_le_bytes())),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonicalise_nucleotide() {
        // Test basic nucleotides
        assert_eq!(canonicalise_nucleotide(b'A'), b'A');
        assert_eq!(canonicalise_nucleotide(b'C'), b'C');
        assert_eq!(canonicalise_nucleotide(b'G'), b'G');
        assert_eq!(canonicalise_nucleotide(b'T'), b'T');

        // Test lowercase
        assert_eq!(canonicalise_nucleotide(b'a'), b'A');
        assert_eq!(canonicalise_nucleotide(b'c'), b'C');

        // Test ambiguous nucleotides
        assert_eq!(canonicalise_nucleotide(b'R'), b'G');
        assert_eq!(canonicalise_nucleotide(b'Y'), b'C');
        assert_eq!(canonicalise_nucleotide(b'S'), b'G');
        assert_eq!(canonicalise_nucleotide(b'W'), b'A');
        assert_eq!(canonicalise_nucleotide(b'K'), b'G');
        assert_eq!(canonicalise_nucleotide(b'M'), b'C');
        assert_eq!(canonicalise_nucleotide(b'B'), b'C');
        assert_eq!(canonicalise_nucleotide(b'D'), b'G');
        assert_eq!(canonicalise_nucleotide(b'H'), b'C');
        assert_eq!(canonicalise_nucleotide(b'V'), b'G');
        assert_eq!(canonicalise_nucleotide(b'N'), b'C');
    }

    #[test]
    fn test_canonicalise_sequence() {
        let seq = b"ACGTN";
        let canonical = canonicalise_sequence(seq);
        assert_eq!(canonical, b"ACGTC"); // N becomes C

        let seq = b"RYMKWS";
        let canonical = canonicalise_sequence(seq);
        assert_eq!(canonical, b"GCCGAG"); // R->G, Y->C, M->C, K->G, W->A, S->G
    }

    #[test]
    fn test_compute_minimizer_hashes() {
        // Simple sequence test
        let seq = b"ACGTACGTACGT";
        let k = 5;
        let w = 3;

        let hashes = compute_minimizer_hashes(seq, k, w, None);

        // We should have at least one minimizer hash
        assert!(!hashes.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        let short_hashes = compute_minimizer_hashes(short_seq, k, w, None);
        assert!(short_hashes.is_empty());
    }

    #[test]
    fn test_calculate_linguistic_complexity() {
        // Test maximum complexity (all different nucleotides)
        let max_complexity_kmer = b"ACGT";
        let complexity = calculate_linguistic_complexity(max_complexity_kmer);
        assert!(
            (complexity - 1.0).abs() < 0.01,
            "Expected ~1.0, got {}",
            complexity
        );

        // Test minimum complexity (homopolymer)
        let min_complexity_kmer = b"AAAA";
        let complexity = calculate_linguistic_complexity(min_complexity_kmer);
        assert!(
            complexity < 0.1,
            "Expected very low complexity, got {}",
            complexity
        );

        // Test moderate complexity (alternating pattern)
        let alt_complexity_kmer = b"ATATAT";
        let complexity = calculate_linguistic_complexity(alt_complexity_kmer);
        assert!(
            complexity > 0.05 && complexity < 0.5,
            "Expected moderate complexity, got {}",
            complexity
        );

        // Test empty k-mer (should also never be filtered)
        let empty_kmer = b"";
        let complexity = calculate_linguistic_complexity(empty_kmer);
        assert_eq!(complexity, 1.0);

        // Test k-mers less than 4 bases (should return 1.0 to never be filtered)
        let single_kmer = b"A";
        let complexity = calculate_linguistic_complexity(single_kmer);
        assert_eq!(complexity, 1.0);

        let dinucleotide_kmer = b"AT";
        let complexity = calculate_linguistic_complexity(dinucleotide_kmer);
        assert_eq!(complexity, 1.0);

        let trinucleotide_kmer = b"ACG";
        let complexity = calculate_linguistic_complexity(trinucleotide_kmer);
        assert_eq!(complexity, 1.0);
    }

    #[test]
    fn test_complexity_filtering() {
        let seq = b"AAAAAAAAAACGTACGTACGT"; // Low complexity start, high complexity end
        let k = 10;
        let w = 4; // k + w - 1 = 13 (odd)

        // Without filtering
        let hashes_no_filter = compute_minimizer_hashes(seq, k, w, None);

        // With high threshold (should filter out low complexity k-mers)
        let hashes_filtered = compute_minimizer_hashes(seq, k, w, Some(0.8));

        // Filtered should have fewer or equal hashes
        assert!(hashes_filtered.len() <= hashes_no_filter.len());
    }
}

use packed_seq::AsciiSeq;
use xxhash_rust::xxh3;

pub const DEFAULT_KMER_LENGTH: usize = 31;
pub const DEFAULT_WINDOW_SIZE: usize = 15;

/// Check if nucleotide is valid ACGT (case insensitive)
#[inline]
fn is_valid_acgt(nucleotide: u8) -> bool {
    matches!(nucleotide, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
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
pub fn compute_minimizer_hashes(seq: &[u8], kmer_length: usize, window_size: usize) -> Vec<u64> {
    let mut hashes = Vec::new();
    fill_minimizer_hashes(seq, kmer_length, window_size, &mut hashes);
    hashes
}

/// Fill a vector with minimizer hashes, skipping k-mers with non-ACGT bases
pub fn fill_minimizer_hashes(
    seq: &[u8],
    kmer_length: usize,
    window_size: usize,
    hashes: &mut Vec<u64>,
) {
    hashes.clear();

    // Skip if sequence is too short
    if seq.len() < kmer_length {
        return;
    }

    let canonical_seq = canonicalise_sequence(seq);

    // Get minimizer positions using simd-minimizers
    let mut positions = Vec::new();
    simd_minimizers::canonical_minimizer_positions(
        AsciiSeq(&canonical_seq),
        kmer_length,
        window_size,
        &mut positions,
    );

    // Filter positions to only include k-mers with ACGT bases in original sequence
    let valid_positions: Vec<u32> = positions
        .into_iter()
        .filter(|&pos| {
            let pos_usize = pos as usize;
            let kmer = &seq[pos_usize..pos_usize + kmer_length];
            kmer_contains_only_acgt(kmer)
        })
        .collect();

    hashes.extend(
        simd_minimizers::iter_canonical_minimizer_values(
            AsciiSeq(&canonical_seq),
            kmer_length,
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

        let hashes = compute_minimizer_hashes(seq, k, w);

        // We should have at least one minimizer hash
        assert!(!hashes.is_empty());

        // Test with a sequence shorter than k
        let short_seq = b"ACGT";
        let short_hashes = compute_minimizer_hashes(short_seq, k, w);
        assert!(short_hashes.is_empty());
    }
}

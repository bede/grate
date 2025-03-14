use serde::{Deserialize, Serialize};

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
            format_version: 1,
            kmer_length: kmer_length as u8,
            window_size: window_size as u8,
        }
    }

    /// Validate header
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.format_version != 1 {
            return Err(anyhow::anyhow!(
                "Unsupported file format version: {}",
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_creation() {
        let header = IndexHeader::new(31, 21);

        assert_eq!(header.format_version, 1);
        assert_eq!(header.kmer_length(), 31);
        assert_eq!(header.window_size(), 21);
    }

    #[test]
    fn test_header_validation() {
        // Valid header
        let valid_header = IndexHeader {
            format_version: 1,
            kmer_length: 31,
            window_size: 21,
        };
        assert!(valid_header.validate().is_ok());

        // Invalid format version
        let invalid_header = IndexHeader {
            format_version: 2, // Unsupported version
            kmer_length: 31,
            window_size: 21,
        };
        assert!(invalid_header.validate().is_err());
    }
}

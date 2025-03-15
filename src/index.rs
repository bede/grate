use crate::index_format::IndexHeader;
use anyhow::{Context, Result};
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

// Re-export index operations
pub use crate::index_build::build;
pub use crate::index_diff::diff;
pub use crate::index_info::info;
pub use crate::index_union::union;

/// Load the hashes without spiking memory usage with an extra vec
pub fn load_minimizer_hashes<P: AsRef<Path>>(path: &P) -> Result<(FxHashSet<u64>, IndexHeader)> {
    let file =
        File::open(path).context(format!("Failed to open index file {:?}", path.as_ref()))?;
    let mut reader = BufReader::new(file);

    // Deserialize header
    let header: IndexHeader =
    decode_from_std_read(&mut reader, bincode::config::standard())
        .context("Failed to deserialize index header")?;
    header.validate()?;

    // Deserialize the count of minimizers so we can init a FxHashSet with the right capacity
    let count: usize =
        decode_from_std_read(&mut reader, bincode::config::standard())
            .context("Failed to deserialize minimizer count")?;

    // Pre-allocate the FxHashSet with the right capacity
    let mut minimizers = FxHashSet::with_capacity_and_hasher(count, Default::default());

    // Populate the FxHashSet
    for _ in 0..count {
        let hash: u64 = decode_from_std_read(&mut reader, bincode::config::standard())
            .context("Failed to deserialize minimizer hash")?;
        minimizers.insert(hash);
    }

    Ok((minimizers, header))
}

/// Helper function to write minimizers to output file or stdout
pub fn write_minimizers(
    minimizers: &FxHashSet<u64>,
    header: &IndexHeader,
    output_path: Option<&PathBuf>,
) -> Result<()> {
    // Create writer based on output path
    let writer: Box<dyn Write> = if let Some(path) = output_path {
        if path.to_string_lossy() == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(
                File::create(path).context("Failed to create output file")?,
            ))
        }
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // Serialize header and minimizers
    let mut writer = BufWriter::new(writer);
    encode_into_std_write(header, &mut writer, bincode::config::standard())
        .context("Failed to serialize index header")?;

    // Serialize the count of minimizers first
    let count = minimizers.len();
    encode_into_std_write(&count, &mut writer, bincode::config::standard())
        .context("Failed to serialize minimizer count")?;

    // Serialize each minimizer directly
    for &hash in minimizers {
        encode_into_std_write(&hash, &mut writer, bincode::config::standard())
            .context("Failed to serialize minimizer hash")?;
    }
    Ok(())
}

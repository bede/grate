# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.0] - 2025-07-08

### Added
- `deacon index diff` optionally accepts a fastx file or stream in place of a second index. This enables index masking using massive sequence collections without the need to first index them.

### Changed
- Deacon uses the recently added `simd-minimizers::iter_canonical_minimizer_values()`, increasing filtering speed by up to 50% on Linux/x86_64 systems. Speeds exceeding 500Mbp/s are now possible.
  - Index format is now version 2. Existing indexes must be rebuilt for use with this version. A new version of the panhuman-1 index is available from Zenodo and object storage. Attempting to load an incompatible index throws an error.
- Position-dependent IUPAC ambiguous base canonicalisation was replaced with a simpler and faster fixed mapping, meaning that records containing ambiguous IUPAC bases may be classified differently to before.
- `deacon index union` now automatically preallocates the required hash table capacity, eleminating slowdowns when combining indexes.
- Compatible minimizer _k_ and _w_ is now validated (k+w-1 must be odd) prior to indexing.
- Default index capacity is now 400M (Was 500).




## [0.6.0] - 2025-06-25

### Added

- Support for the .xz compression format via liblzma.
- Adjustable filter output compression level with `--compression-level`.
- Report fields `seqs_out_proportion` and `bp_out_proportion`.

### Changed

- Use zlib-rs for much faster gzip decompression.
- Displays number and proportion of _retained_ reads and base pairs during filtering.

## [0.5.0] - 2025-06-11

### Added

- `--deplete` (`-d`) flag to remove index matches.
- Support for relative thresholds (floats between 0.0 and 1.0) for required minimizer hits to `--matches`.
- `-O` short argument name for `--output2`.
- Tests.

### Changed

- Default filtering behaviour now _passes_ index matches. Use `--deplete` (`-d`) to remove matches.
- Renamed `--nucleotides` (`-n`) to `--prefix-length` (`-p`).
- Renamed `--report` to `--summary`.

### Removed

- `--invert` argument has been removed.

## [0.4.0] - 2025-05-23

### Added

- Non-interleaved paired output file support (`--output2`).

### Changed

- Faster indexing.
- Renamed `--log` argument to `--report`.
- Filter stats are now always sent to stderr (a json report can be written wherever one chooses).
- For paired input sequences, identical minimizer hits in both mates of a read pair are now counted only once.

## [0.3.0] - 2025-05-09

### Added

- Parallel filtering.
  - Up to 10x faster from initial testing.
  - Configurable with new `--threads` parameter.
  - Uses available CPU cores by default (`--threads 0`).
- Tests.

### Changed

- Default minimizer parameters changed to k=31 and w=15.

## [0.2.0] - 2025-04-23

### Added

- Paired read support using either third positional argument or interleaved stdin.

### Changed

- Faster indexing.
- More accurate default parameters.
- Optional argument changes.
- Refactored lib.rs.
- Dependency updates.
  - Bincode2.

## [0.1.0] - 2025-03-14

### Added

- Initial experimental release.

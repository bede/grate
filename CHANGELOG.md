# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

# Grate

Fast estimation of target reference DNA sequence containment (~coverage) and median abundance (~depth).

## Install

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/bede/grate
```

## Usage

```bash
# One samples
grate cov refs.fa reads.fastq.gz

# Many samples
grate cov refs.fa reads1.fastq.gz reads2.fastq.gz reads3.fastq.gz…

# CSV output, needed for plotting
grate cov -f csv refs.fa reads.fastq.gz > results.csv

# Plotting (requires uv)
./plot.py results.csv
uv run plot.py results.csv

# See plotting options
./plot.py --help
```

The plotting script is a self-contained uv script that automatically installs its dependencies (pandas, altair, vl-convert-python) when run. You need [uv](https://docs.astral.sh/uv/) installed.

### CLI Reference

```bash
$ grate cov -h
Estimate containment and abundance of target sequence(s) in read file(s) or stream

Usage: grate cov [OPTIONS] <TARGETS> <READS>...

Arguments:
  <TARGETS>   Path to fasta file containing target sequence record(s)
  <READS>...  Path(s) to fastx file(s) containing reads (or - for stdin). Multiple files are treated as separate samples

Options:
      --sample-names <SAMPLE_NAMES>
          Sample names for read files (optional, defaults to filename without extensions)
  -k, --kmer-length <KMER_LENGTH>
          Minimizer length (1-61) [default: 31]
  -w, --window-size <WINDOW_SIZE>
          Minimizer window size [default: 15]
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
  -o, --output <OUTPUT>
          Path to output file (- for stdout) [default: -]
  -q, --quiet
          Suppress progress reporting
  -f, --format <FORMAT>
          Output format [default: table] [possible values: table, csv, json]
  -a, --abundance-thresholds <ABUNDANCE_THRESHOLDS>
          Comma-separated abundance thresholds for containment calculation [default: 10]
  -d, --discriminatory
          Retain only minimizers exclusive to each target
  -l, --limit <LIMIT>
          Terminate read processing after approximately this many bases (e.g. 50M, 10G)
      --sort <SORT>
          Sort order for results: o=original (default), a=alphabetical, c=containment (max first) [default: o] [possible values: o, a, c]
  -h, --help
          Print help
```


# Grate

A fork of [Deacon](https://github.com/bede/deacon).

## Install

```bash
RUSTFLAGS="-C target-cpu=native" cargo install --git https://github.com/bede/grate
```

## Usage

```bash
# One sample
grate con refs.fa reads.fastq.gz

# Many samples
grate con refs.fa reads1.fastq.gz reads2/ reads3.fa.zst…

# Stdin
zstdcat reads3.fq.zst | grate con refs.fa

# CSV output for plotting
grate con -f csv refs.fa reads.fastq.gz > results.csv

# Plotting (requires uv)
uv run plot/con.py results.csv

# View plotting options
uv run ./plot/con.py --help
```

Run the plotting scripts with [uv](https://docs.astral.sh/uv/) to automatically handle dependencies.

### CLI Reference

```bash
$ grate con -h
Calculate minimizer containment & abundance in fastx files or directories thereof

Usage: grate con [OPTIONS] <TARGETS> <SAMPLES>...

Arguments:
  <TARGETS>     Path to fasta file containing target sequence record(s)
  <SAMPLES>...  Path(s) to fastx files/dirs (- for stdin). Each file/dir is treated as a separate sample

Options:
  -k, --kmer-length <KMER_LENGTH>
          Minimizer length (1-61) [default: 31]
  -w, --window-size <WINDOW_SIZE>
          Minimizer window size [default: 31]
  -a, --abundance-thresholds <ABUNDANCE_THRESHOLDS>
          Comma-separated abundance thresholds for containment calculation [default: 10]
  -d, --discriminatory
          Consider only minimizers unique to each target
  -t, --threads <THREADS>
          Number of execution threads (0 = auto) [default: 8]
  -l, --limit <LIMIT>
          Terminate read processing after approximately this many bases (e.g. 50M, 10G)
  -o, --output <OUTPUT>
          Path to output file (- for stdout) [default: -]
  -f, --format <FORMAT>
          Output format [default: table] [possible values: table, csv, json]
  -n, --names <SAMPLE_NAMES>
          Comma-separated sample names (default is file/dir name without extension)
  -s, --sort <SORT>
          Sort displayed results: o=original, t=target, s=sample, c=containment (descending) [default: o] [possible values: o, t, s, c]
  -q, --quiet
          Suppress progress reporting
  -h, --help
          Print help
```


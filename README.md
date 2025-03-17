[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)

# Deacon

*Currently under development*

A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for efficient host depletion. Default behaviour removes query sequences with one or more minimizers present in the index. Filters at 50Mbp/s using a single Apple M1 core and indexes the human genome in under 60s. Peak memory usage is ~2.5GB for a human genome with default parameters.

The sensitivity/specificity/memory tradeoff can be tuned using *k*-mer length (`-k`), minimizer window size (`-w`), and match threshold (`-m`). With long reads, speeds of hundreds of megabases per second are possible by considering only the first `-n` bases of each query sequence. Deacon uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers) enabling fast vectorised minimizer calculation. This project is currently experimental and unstable.



## Install

### Cargo
[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)


```
cargo install deacon
```



## Usage

### Indexing

Supports FASTA[.gz] input files and outputs to stdout or file

```bash
deacon index build chm13v2.fa > human.k31w21.idx
```

```
deacon index build -k 41 -m 27 chm13v2.fa > human.k41w27.idx
```



### Filtering

Supports FASTA or FASTQ input from stdin or file and outputs to stdout or file. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically. To avoid a compression bottleneck with gzip, you should redirect uncompressed FASTA/Q to pigz (see below).

``` bash
deacon filter human.k31w21.idx reads.fq.gz -o filt.fq.gz  # Basic usage
```

```bash
zcat reads.fq.gz | deacon filter -n 1000 human.k31w21.idx | pigz > filt.fq.gz  # Fast
```

```bash
deacon filter -m 3 human.k31w21.idx | pigz > filt.fq.gz  # Greater precision with -m 3
```



### Set operations on indexes

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to nonredundantly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.

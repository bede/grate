[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)

# Deacon

*Currently under development*

A fast minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for efficient host depletion (*deacon*-tamination). Default behaviour removes query sequences with two or more minimizers present in the index. Filters at 50Mbp/s using a single Apple M1 core and indexes the human genome in under 60s. Peak memory usage is ~2.5GB for a human genome with default parameters.

The sensitivity/specificity/memory tradeoff can be tuned using *k*-mer length (`-k`), minimizer window length (`-w`), and match threshold (`-m`). With long reads, speeds of hundreds of megabases per second are possible by considering only the first `-n` bases of each query sequence. Uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers) for fast vectorised minimizer calculation. This project is currently unstable.



## Install

### Cargo
[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)


```
cargo install deacon
```



## Usage

### Indexing

Supports FASTA[.gz] input files and outputs to stdout or file (`-o`).

```bash
deacon index build chm13v2.fa > human.k31w21.idx
deacon index build -k 41 -m 27 chm13v2.fa > human.k41w27.idx  # Custom minimizer k and w
```



### Filtering

Supports FASTA or FASTQ input from stdin or file and outputs to stdout or file. Paired sequences are supported as either separate files or interleaved stdin, and are  written in interleaved format to either stdout or file. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically. Consider piping uncompressed FASTA/Q to pigz to avoid compression bottlenecks when writing gzip output.

``` bash
deacon filter human.idx reads.fq.gz -o filt.fastq  # File input & output
zcat reads.fq.gz | deacon filter human.idx | pigz > filt.fq.gz  # Fast gzip
zcat reads.fq.gz | deacon -n 1000 filter human.idx | pigz > filt.fq.gz  # Faster
deacon filter -m 3 human.idx reads.fq.gz | pigz > filt.fq.gz  # More precise
deacon filter -m 1 human.idx reads.fq.gz | pigz > filt.fq.gz  # More sensitive
deacon filter human.idx r1.fq.gz r2.fq.gz > filt12.fastq  # Paired file input
zcat r12.fq.gz | deacon filter human.idx - - > filt12.fastq  # Interleaved stdin
deacon filter human.idx reads.fq.gz --log log.json > filt.fq  # Log results JSON
```



## Reports

Use `--log results.json` to save filtering results to a JSON file:
```
{
  "version": "0.2.0",
  "index": "data/chm13v2.k31w21.idx",
  "input1": "data/viruses/ont.100MB.fastq.gz",
  "input2": null,
  "output": "-",
  "k": 31,
  "w": 21,
  "m": 3,
  "n": 0,
  "invert": false,
  "rename": false,
  "seqs_in": 190000,
  "seqs_out": 189628,
  "seqs_removed": 372,
  "seqs_removed_proportion": 0.001957894736842105,
  "bp_in": 356014239,
  "bp_out": 354404228,
  "bp_removed": 1610011,
  "bp_removed_proportion": 0.004522321928814763,
  "time": 12.3978525,
  "seqs_per_second": 15325,
  "bp_per_second": 28715798
}
```



### Set operations on indexes

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to nonredundantly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.

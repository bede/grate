[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)

# Deacon

A minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for efficient host decontamination. Default behaviour removes query sequences with one or more minimizers present in the index. Filters at 50Mbp/s using a single Apple M1 core. Indexing the human genome takes 90s. Peak memory usage is ~2.5GB for a human genome with default parameters.

The sensitivity/specificity/memory tradeoff can be tuned using kmer length (`-k`), minimizer window size (`-w`), and match threshold (`-m`). With long reads, speed can be increased to hundreds of megabases per second by considering the first `-n` bases of each query sequence. This project is currently experimental and unstable.



## Install

### Cargo
[![Crates.io](https://img.shields.io/crates/v/deacon.svg)](https://crates.io/crates/deacon)


```
cargo install deacon
```



## Usage

### Indexing

```bash
deacon index build chm13v2.fa > human.k31w21.idx
```

```
deacon index build -k 41 -m 27 chm13v2.fa > human.k41w27.idx
```



### Filtering

``` bash
deacon filter human.k31w21.idx reads.fq.gz -o filt.fq.gz  # Basic usage
```

```bash
zcat reads.fq.gz | deacon filter -n 1000 human.k31w21.idx | pigz > filt.fq.gz  # Fast
```

```bash
deacon filter -m 3 human.k31w21.idx | pigz > filt.fq.gz  # Greater precision with -m 3
```

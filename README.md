[![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon) [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon) ![CI status](https://img.shields.io/github/actions/workflow/status/bede/deacon/test.yml?style=flat-square) [![Downloads](https://img.shields.io/conda/dn/bioconda/deacon.svg?style=flat-square)](https://anaconda.org/bioconda/deacon)


# Deacon

<div align="center"><img src="deacon.png" width="180" alt="Logo"></div>

Fast minimizer-based filtering of nucleotide sequences in FASTA or FASTQ format for search or depletion of files or streams. Default parameters balance sensitivity and specificity for microbial (meta)genomic host depletion, for which a validated prebuilt index is available. Sensitivity, specificity and required memory can be tuned by varying *k*-mer length (`-k`), minimizer window size (`-w`), and the number or proportion of required index matches (`-m`) per query. Minimizer `k` and `w`  are chosen at index time, while the match threshold `m` can be changed at filter time.

Building on [simd-minimizers](https://github.com/rust-seq/simd-minimizers), Deacon is capable of filtering long reads at >250Mbp/s (Apple M4) and indexing a human genome in <30s. Peak memory usage during filtering is 5GB for the default panhuman index. Partial query matching can be used to further increase speed for long queries by considering only the first `-p` bases per query. Sequences can optionally be renamed for privacy and smaller file sizes. Deacon reports filtering performance in real time and optionally writes a JSON summary. Stay tuned for a preprint evaluating performance and further improvements.

## Install

### conda/mamba/pixi  [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

### cargo [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
cargo install deacon
```

## Usage

> [!IMPORTANT]
> Version 0.5.0 introduced breaking CLI changes. Use `--deplete` (`-d`) to remove sequences matching the index.

### Indexing

Build indexes with `deacon index build`. For human host depletion, the prebuilt validated panhuman index is recommended, available for download below. Object storage is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford.

```shell
deacon index build chm13v2.fa > human.k31w15.idx
```

#### Prebuilt indexes

|                           Name/URL                           |                         Composition                          | Minimizers                   | Subtracted minimizers | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ---------------------------- | --------------------- | ----- | ------- |
| [**panhuman-1**](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/panhuman-1.k31w15.idx) | ([HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40)) - bacteria (FDA-ARGOS)  - viruses (RefSeq) | 409,914,298 (*k*=31, *w*=15) | 20,741 (**0.0051%**)  | 3.7GB | 2025-04 |

### Filtering

The command `deacon filter` accepts an index path followed by up to two query FASTA/FASTQ file paths, depending on whether query sequences originate from stdin, a single file, or paired input files. Paired queries are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file, or else to separate paired output files. For paired reads, distinct minimizer hits originating from either mate are counted. By default, query sequences with two or more minimizer hits to the index (`-m 2`) pass the filter. Filtering can be inverted for e.g. host depletion using the `--deplete` (`-d`) flag. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically by file extension. Consider using Zstandard compression rather than Gzip if possible for best performance.

**Examples**

```bash
# Keep only SARS-CoV-2 reads
deacon filter -o sars2.fq.gz MN908947.idx reads.fq.gz 

# Host depletion using the panhuman-1 index
deacon filter -d -o filt.fq.gz panhuman-1.k31w15.idx reads.fq.gz

# Host depletion using stdin and stdout
zcat reads.fq.gz | deacon filter -d panhuman-1.k31w15.idx > filt.fq

# Faster Zstandard compression
deacon filter -d -o filt.fq.zst panhuman-1.k31w15.idx reads.fq.zst 

# More sensitive match threshold of at 1 minimizer hit
deacon filter -d -m 1 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# More specific match threshold of 50% minimizer hits (minimum 1)
deacon filter -d -m 0.5 panhuman-1.k31w15.idx reads.fq.gz > filt.fq.gz

# Deplete paired reads
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fastq
deacon filter -d panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz -o filt.r1.fq.gz -O filt.r2.fq.gz
zcat r12.fq.gz | deacon filter -d panhuman-1.k31w15.idx - - > filt12.fq

# Save summary JSON
deacon filter -d -o filt.fq.gz --summary summary.json panhuman-1.k31w15.idx reads.fq.gz
```

### Composing indexes with set operations

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to succinctly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.

## Filtering summary statistics

Use `--summary summary.json` to save detailed filtering statistics:
```json
{
  "version": "deacon 0.5.0",
  "index": "panhuman-1.k31w15.idx",
  "input1": "HG02334.1m.fastq.gz",
  "input2": null,
  "output": "-",
  "output2": null,
  "k": 31,
  "w": 21,
  "match_threshold": "2",
  "prefix_length": 0,
  "deplete": true,
  "rename": false,
  "seqs_in": 1000000,
  "seqs_out": 13452,
  "seqs_removed": 986548,
  "seqs_removed_proportion": 0.986548,
  "bp_in": 5477122928,
  "bp_out": 5710050,
  "bp_removed": 5471412878,
  "bp_removed_proportion": 0.9989574727324798,
  "time": 125.755103875,
  "seqs_per_second": 7951,
  "bp_per_second": 43553881
}
```

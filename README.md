[![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon) [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

# Deacon

A fast general purpose minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for rapid and accurate host depletion. Default parameters have been selected to maximise balanced accuracy for short and long reads. Sensitivity, specificity and memory use may nevertheless be tuned by varying *k*-mer length (`-k`), minimizer window length (`-w`), and the number of required index matches (`-m`) per query. Minimizer `k` and `w`  are chosen at index time, while the number of required matches `m` can be specified at filter time.

Building on [simd-minimizers](https://github.com/rust-seq/simd-minimizers), Deacon is capable of filtering at >200Mbp/s (Apple M1) and indexing a human genome in under 60s. Peak memory usage is ~4.5GB for the default panhuman index. Partial query matching can be used to further increase speed for long queries by considering only the first `-n` bases per query. Stay tuned for comprehensive validation and benchmarks. This project is currently unstable.

## Install

### conda/mamba/pixi (recommended)  [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

```bash
conda install -c bioconda deacon
```

### cargo [![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon)

```bash
cargo install deacon
```

## Usage

### Indexing

Custom indexes can be built using `deacon index build`. For human host depletion, the prebuilt validated panhuman index is recommended, available for download below. Object storage is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford.

```shell
deacon index build chm13v2.fa > human.k31w15.idx
```

#### Prebuilt indexes

|                           Name/URL                           |                         Composition                          | Minimizers                   | Masked minimizers    | Size  | Date    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | ---------------------------- | -------------------- | ----- | ------- |
| [**panhuman-1**](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/panhuman-1.k31w15.idx) | ([HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) ∪ [CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) ∪ [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40)) - bacteria (FDA-ARGOS)  - viruses (RefSeq) | 409,914,298 (*k*=31, *w*=15) | 20,741 (**0.0051%**) | 3.7GB | 2025-04 |

### Filtering

The command `deacon filter` accepts an index path followed by up to two query FASTA/FASTQ file paths, depending on whether query sequences originate from stdin, a single file, or paired input files. Paired queries are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file. For paired reads, hits in either mate count towards the threshold for filtering that pair. By default, query sequences with fewer than two index hits are sent to stdout. Filtering can be inverted using the `--invert` flag. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically by file extension. Since (de)compression can be rate limiting, consider using Zstandard rather than Gzip for best performance.

```bash
deacon filter panhuman-1.k31w15.idx reads.fq.gz -o filt.fq  # File input & output
zcat reads.fq.gz | deacon filter panhuman-1.k31w15.idx > filt.fq  # Stdin and stdout
deacon filter panhuman-1.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # Parallel gzip
deacon -n 1000 filter panhuman-1.k31w15.idx reads.fq.zst | zstd > filt.fq.zst  # Fastest
deacon filter -m 3 panhuman-1.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # More precise
deacon filter -m 1 panhuman-1.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # More sensitive
deacon filter panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fastq  # Paired file input
deacon filter panhuman-1.k31w15.idx r1.fq.gz r2.fq.gz --output filt.r1.fq.gz --output2 filt.r2.fq.gz  # Paired file input/output
zcat r12.fq.gz | deacon filter panhuman-1.k31w15.idx - - > filt12.fastq  # Interleaved stdin
deacon filter panhuman-1.k31w15.idx reads.fq.gz --report report.json > filt.fq  # Save report JSON
```

## Reports

Use `--report results.json` to save a filtering report to JSON file:
```json
{
  "version": "0.3.0",
  "index": "data/chm13v2.k31w21.idx",
  "input1": "data/HG02334.1m.fastq.gz",
  "input2": null,
  "output": "-",
  "k": 31,
  "w": 21,
  "m": 1,
  "n": 0,
  "invert": false,
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

### Composing indexes with set operations

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to succinctly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.

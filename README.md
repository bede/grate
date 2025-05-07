[![Crates.io Version](https://img.shields.io/crates/v/deacon?style=flat-square)](https://crates.io/crates/deacon) [![Bioconda version](https://anaconda.org/bioconda/deacon/badges/version.svg)](https://anaconda.org/bioconda/deacon)

# Deacon

A minimizer-based filter for nucleotide sequences in FASTA or FASTQ format, built for fast host depletion. Default behaviour removes query sequences with two or more minimizers present in the index. Capable of filtering at >200Mbp/s on Apple M1 and indexing a human genome in under 60s. Peak memory usage is ~4.5GB for the default panhuman index.

The sensitivity/specificity/memory tradeoff can be tuned using indexes built with varying *k*-mer length (`-k`), minimizer window length (`-w`), and match threshold (`-m`). Filtering speed may be increased by considering only the first `-n` bases per query sequence. Uses [simd-minimizers](https://github.com/rust-seq/simd-minimizers) for accelerated minimizer computation. This project is currently unstable, but validation and benchmarks will be published soon.

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

The command `deacon filter` accepts a path to an index and a FASTA/FASTQ query from file or stdin. Prebuilt indexes are available for download below, and custom indexes may be created using `deacon index build`.

### Prebuilt indexes

Use `deacon index fetch panhuman-1m` to fetch the default panhuman index from object storage for immediate use with `deacon filter`. Object storage is provided by the [ModMedMicro research unit](https://www.expmedndm.ox.ac.uk/modernising-medical-microbiology) at the University of Oxford.

|                             Name                             |                         Composition                          | Parameters     | Minimizers  | Size  | Date    | Masked minimizers    |
| :----------------------------------------------------------: | :----------------------------------------------------------: | -------------- | ----------- | ----- | ------- | -------------------- |
| [**panhuman-1m**](https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/panhuman-1m.k31w15.idx) | [HPRC Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies/blob/main/assembly_index/Year1_assemblies_v2_genbank.index) + [CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/11828891) + [GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40), masked with bacteria (`argos988`) and viruses (`rsviruses17900`) | *k*=31, *w*=15 | 409,914,298 | 3.7GB | 2025-05 | 20,741 (**0.0051%**) |

### Filtering

Supports FASTA or FASTQ input from file or stdin and outputs to stdout or file. Paired sequences are supported as either separate files or interleaved stdin, and written interleaved to either stdout or file. Gzip (.gz) and Zstandard (.zst) compression formats are detected automatically by file extension. Since (de)compression can be rate limiting, consider using Zstandard rather than Gzip for best performance on multicore systems.

```bash
deacon filter panhuman-1m.k31w15.idx reads.fq.gz -o filt.fq  # File input & output
zcat reads.fq.gz | deacon filter panhuman-1m.k31w15.idx > filt.fq  # Stdin and stdout
deacon filter panhuman-1m.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # Parallel gzip
deacon -n 1000 filter panhuman-1m.k31w15.idx reads.fq.zst | zstd > filt.fq.zst  # Fastest
deacon filter -m 3 panhuman-1m.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # More precise
deacon filter -m 1 panhuman-1m.k31w15.idx reads.fq.gz | pigz > filt.fq.gz  # More sensitive
deacon filter panhuman-1m.k31w15.idx r1.fq.gz r2.fq.gz > filt12.fastq  # Paired file input
zcat r12.fq.gz | deacon filter panhuman-1m.k31w15.idx - - > filt12.fastq  # Interleaved stdin
deacon filter panhuman-1m.k31w15.idx reads.fq.gz --log log.json > filt.fq  # Log results JSON
```

## Reports

Use `--log results.json` to save a filtering summary to a JSON file:
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

### Building indexes

 Build custom indexes using`deacon index build`. This accepts FASTA[.gz] input files and outputs to stdout or file (`-o`). 

```bash
deacon index build chm13v2.fa > human.k31w21.idx
deacon index build -k 41 -m 27 chm13v2.fa > human.k41w27.idx  # Custom minimizer k and w
```

### Composing indexes with set operations

- Use `deacon index union 1.idx 2.idx > 1+2.idx` to nonredundantly combine two (or more) deacon minimizer indexes.
- Use `deacon index diff 1.idx 2.idx > 1-2.idx` to subtract minimizers in 2.idx from 1.idx. Useful for masking.

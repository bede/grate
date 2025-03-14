# Deacon

A fast general purpose DNA sequence filter built for host decontamination. Accurately filters human sequences at 50Mbp/s on an Apple silicon using a single thread. Indexing the human genome takes ~3 minutes. Peak memory usage is 2.5GB with the default index.

The sensitivity/specificity/memory tradeoff can be tuned using kmer length (`-k`), minimizer window size (`-w`), and match threshold (`-m`). With long reads, speed can be increased to hundreds of megabases per second by considering the first `-n` bases of each query sequence. Deacon is being actively developed, with benchmarks and a preprint to follow.



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
deacon filter -m 3 human.k31w21.idx | pigz > filt.fq.gz  # Greater precision
```


# Deacon

A fast general purpose DNA sequence filter built for host decontamination. Accurately filters human sequences at 40-50Mbp/s on an Apple M1 using a single thread, and much faster when matching only partial queries. Indexing the human genome takes 3 minutes. Benchmarks and preprint coming soon.



## Usage

### Indexing

```bash
deacon index build chm13v2.fa > human.k31w21.idx  # Three minutes or so
```

### Filtering

``` bash
deacon filter human.k31w21.idx reads.fq.gz -o filt.fq.gz  # Basic usage
```

```bash
zcat reads.fq.gz | deacon filter -n 1000 human.k31w21.idx | pigz > filt.fq.gz  # Fast
```


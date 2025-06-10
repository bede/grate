# Host depletion benchmarking

## Fetch datasets

### Human (`chm13v2`)

Baseline human genome (chm13v2)

```shell
wget -O chm13v2.fastq.gz https://zenodo.org/records/15424966/files/chm13v2.fastq.gz
wget -O chm13v2.r1.fastq.gz https://zenodo.org/records/15424966/files/chm13v2.r1.fastq.gz
wget -O chm13v2.r2.fastq.gz https://zenodo.org/records/15424966/files/chm13v2.r2.fastq.gz
```



### Bacteria (`argos988`)

These files are compressed with Zstandard rather than Gzip in order to squeeze inside Zenodo's 50GB limit.

```shell
wget -O argos988.r1.fastq.zst https://zenodo.org/records/15424142/files/argos988.r1.fastq.zst
wget -O argos988.fastq.zst https://zenodo.org/records/15424142/files/argos988.fastq.zst
wget -O argos988.r2.fastq.zst https://zenodo.org/records/15424142/files/argos988.r2.fastq.zst

# Convert zstandard to gzip
zstdcat argos988.fastq.zst | pigz > argos988.fastq.gz
zstdcat argos988.r1.fastq.zst | pigz > argos988.r1.fastq.gz
zstdcat argos988.r2.fastq.zst | pigz > argos988.r2.fastq.gz
```



### Viruses (`rsviruses17900`)

```shell
wget -O rsviruses17900.r1.fastq.gz https://zenodo.org/records/15411280/files/rsviruses17900.r1.fastq.gz
wget -O rsviruses17900.r2.fastq.gz https://zenodo.org/records/15411280/files/rsviruses17900.r2.fastq.gz
wget -O rsviruses17900.fastq.gz https://zenodo.org/records/15411280/files/rsviruses17900.fastq.gz
```



## Fetch references/indexes

```shell
wget -O panhuman-1.k31w15.idx https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/panhuman-1.k31w15.idx
```


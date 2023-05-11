# PREREQUISITE PREPARATION

Following GATK4 best practices workflow: **_[GATK_Germline_Short_Variant_Discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)_**

## Table of Contents

- [PREREQUISITE PREPARATION](#prerequisite-preparation)
  - [Table of Contents](#table-of-contents)
  - [Tools](#tools)
    - [1. bcftools](#1-bcftools)
    - [2. htslib](#2-htslib)
    - [3. GATK](#3-gatk)
  - [Reference genome for human chromosome 21 and necessary databases](#reference-genome-for-human-chromosome-21-and-necessary-databases)
    - [1. Reference genome](#1-reference-genome)
    - [2. dbsnp\_146 database](#2-dbsnp_146-database)
    - [3. 1000G\_omni2.5.hg38 known snps database](#3-1000g_omni25hg38-known-snps-database)
    - [4. 1000G\_phase1 snps with high confidence hg38 database](#4-1000g_phase1-snps-with-high-confidence-hg38-database)
    - [5. Homo sapiens assembly38 known indels database](#5-homo-sapiens-assembly38-known-indels-database)
    - [6. Mills and 1000G gold standard for known indels hg38 database](#6-mills-and-1000g-gold-standard-for-known-indels-hg38-database)
  - [Dataset for practice](#dataset-for-practice)

## Tools

### 1. bcftools

Please download and install **bcftools** from [here](https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2) to manipulate vcf file.

Please also refer to this [htslib_Guide](http://www.htslib.org/download/) or [latest release guide](https://github.com/samtools/bcftools/releases/latest/) on how to build and install appropriately.

### 2. htslib

Please download and install **htslib** from [here](https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2).

Please have a look at this [htslib_latest](https://github.com/samtools/htslib/releases/tag/1.17) guide on how to install it correctly.

### 3. GATK

This variant discovery practice will be using [GATK v4.4.0.0](https://github.com/broadinstitute/gatk/releases).

One can download the zip file for this latest version [here](https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip).

Please have a look at this in depth tutorial in GATK github: [broadinstitute/gatk](https://github.com/broadinstitute/gatk).

1. [Requirements](https://github.com/broadinstitute/gatk#requirements)
2. [Downloading GATK](https://github.com/broadinstitute/gatk#downloading)
3. [Building GATK](https://github.com/broadinstitute/gatk#building-gatk4)
4. [Running GATK](https://github.com/broadinstitute/gatk#running-gatk4)

## Reference genome for human chromosome 21 and necessary databases

Please download the reference genome and databases required for this variant calling practice. One should always download both the **vcf file** and its **tbi indexed file**.

### 1. Reference genome

One can navigate to UCSC sequence data by chromosome using this link: [UCSC_hg38_sequence_data_by_chromosome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/), then choose to download "chr21.fa.gz" reference fasta file.

OR

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
```

### 2. dbsnp_146 database

- VCF file

```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
```

- INDEX file

```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

### 3. 1000G_omni2.5.hg38 known snps database

- VCF file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
```

- INDEX file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
```

### 4. 1000G_phase1 snps with high confidence hg38 database

- VCF file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```

- INDEX file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

### 5. Homo sapiens assembly38 known indels database

- VCF file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
```

- INDEX file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```

### 6. Mills and 1000G gold standard for known indels hg38 database

- VCF file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

- INDEX file

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```

## Dataset for practice

Please download both **bam** file and **bai** file to follow this practice.
These 2 files can be found in [test_data](./test_data/) folder.

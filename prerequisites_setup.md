# VARIANT CALLING
Following GATK4 best practices workflow: **_[GATK_Germline_Short_Variant_Discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)_**

## Tools

### 1. bcftools
Please download and install [bctools](https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2) to manipulate vcf file.

Please also refer to this [htslib_Guide](http://www.htslib.org/download/) or [latest release guide](https://github.com/samtools/bcftools/releases/latest/) on how to build and install appropriately.

### 2. GATK
This variant discovery practice will be using [GATK v4.4.0.0](https://github.com/broadinstitute/gatk/releases).

One can download the zip file for this latest version [here](https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip).

Please have a look at this in depth tutorial in GATK github: [broadinstitute/gatk](https://github.com/broadinstitute/gatk).

1. [Requirements](https://github.com/broadinstitute/gatk#requirements)
2. [Downloading GATK](https://github.com/broadinstitute/gatk#downloading)
3. [Building GATK](https://github.com/broadinstitute/gatk#building-gatk4)
4. [Running GATK](https://github.com/broadinstitute/gatk#running-gatk4)

## Reference genome for human chromosome 21 and necessary databases

Please download the reference genome and databases required for this variant calling practice. One should always download both the **<u>vcf file</u>** and its **<ins>tbi indexed file</ins>**.

### 1. Reference genome
One can navigate to UCSC sequence data by chromosome using this link: [UCSC_hg38_sequence_data_by_chromosome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/), then choose to download "chr21.fa.gz" reference fasta file.

*OR*

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
```

### 2. dbsnp_146 database
1. VCF file
```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz
```
2. INDEX file
```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

### 3. 1000G_omni2.5.hg38 known snps database 
1. VCF file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
```
2. INDEX file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
```

### 4. 1000G_phase1 snps with high confidence hg38 database
1. VCF file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```
2. INDEX file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

### 5. Homo sapiens assembly38 known indels database
1. VCF file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
```
2. INDEX file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```

### 6. Mills and 1000G gold standard for known indels hg38 database
1. VCF file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```
2. INDEX file
```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```
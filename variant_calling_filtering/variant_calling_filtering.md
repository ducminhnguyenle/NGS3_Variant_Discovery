# VARIANT PROCESSING

---

## Table of Contents

- [VARIANT PROCESSING](#variant-processing)
  - [Table of Contents](#table-of-contents)
  - [VARIANT CALLING](#variant-calling)
    - [1. Workflow of variant calling using HaplotypeCaller](#1-workflow-of-variant-calling-using-haplotypecaller)
      - [Setup working directories](#setup-working-directories)
    - [2. FASTA - Reference genome format](#2-fasta---reference-genome-format)
    - [3. HaplotypeCaller](#3-haplotypecaller)
  - [VARIANT FILTERING (QUALITY-BASED FILTERING)](#variant-filtering-quality-based-filtering)
    - [1. CNNScoreVariants](#1-cnnscorevariants)
    - [2. FilterVariantTranches](#2-filtervarianttranches)
      - [Variant Genotype Filters](#variant-genotype-filters)
    - [3. Hard Filtering](#3-hard-filtering)

## VARIANT CALLING

### 1. Workflow of variant calling using HaplotypeCaller

```mermaid
flowchart TD
    subgraph VARIANT CALLING
    A("Analysis-Ready Reads<br>BAM") -- "Call Variants per sample<br><u>HaplotypeCaller</u> in single sample mode" --> B{{"VCF<br>Raw SNPs-Indels"}}
    style A fill:#34495e,stroke:#333,stroke-width:4px,color:#fff
    style B fill:#767076,stroke:#F9F2F4,stroke-width:2px,color:#fff
    end
    C{{"Reference genome preparation<br>FASTA"}} --"Index\nSequenceDictionary"--> B
    style C fill:#5f51fe,stroke:#F9F2F4,stroke-width:2px,color:#fff

```

#### Setup working directories

- Directories absolute path:

```bash
# Set up your own path to both fasta reference, recalibrated bam file and output variant calling folder
reference="path/to/reference/UCSC_hg38/chr21.fa.gz"
recal_reads="path/to/alignment/chr21_tumor_recal.bam"
output_vcf="path/to/variant_calling"
known_snps="path/to/reference/GATKBundle/known_snps"
known_indels="path/to/reference/GATKBundle/known_indels"
```

- Directories structure:

```text
Variant_Discovery/
|-- alignment/
|    |--chr21_tumor_recal.bai
|    `--chr21_tumor_recal.bam
|-- reference/
|    |-- GATKBundle/
|    |    |-- known_indels/
|    |    |    |-- Homo_sapiens_assembly38.known_indels.vcf.gz
|    |    |    |-- Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
|    |    |    |-- Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
|    |    |    `-- Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
|    |    `-- known_snps/
|    |         |-- 1000G_omni2.5.hg38.vcf.gz
|    |         |-- 1000G_omni2.5.hg38.vcf.gz.tbi
|    |         |-- 1000G_phase1.snps.high_confidence.hg38.vcf.gz
|    |         |-- 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
|    |         |-- dbsnp_146.hg38.vcf.gz
|    |         `-- dbsnp_146.hg38.vcf.gz.tbi
|    `-- UCSC_hg38/
|         |-- chr21.dict
|         |-- chr21.fa.gz
|         |-- chr21.fa.gz.fai
|         `-- chr21.fa.gz.gzi
|-- script/
|    `-- variant_calling_filtering.sh
`-- variant_calling/
     |-- chr21_tumor_raw_variants.vcf
     |-- chr21_tumor_raw_variants.vcf.idz
     |-- chr21_tumor_raw_snps.vcf
     |-- chr21_tumor_raw_snps.vcf.idx
     |-- chr21_tumor_raw_indels.vcf
     |-- chr21_tumor_raw_indels.vcf.idx
     |-- chr21_tumor_filtered_snps.vcf
     |-- chr21_tumor_filtered_snps.vcf.idx
     |-- chr21_tumor_filtered_indels.vcf
     |-- chr21_tumor_filtered_indels.vcf.idx
     |-- chr21_tumor_ready_snps.vcf
     |-- chr21_tumor_ready_snps.vcf.idx
     |-- chr21_tumor_ready_snpsGT.vcf
     |-- chr21_tumor_ready_indels.vcf
     |-- chr21_tumor_ready_indels.vcf.idx
     `-- chr21_tumor_ready_indelsGT.vcf
```

### 2. FASTA - Reference genome format

<https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format>

1. Create the FASTA index file (only bgzip can be indexed)

    ```bash
    reference="path/to/reference/UCSC_hg38/chr21.fa.gz"

    gzip -d ${reference}
    bgzip ${reference%.*} 
    samtools faidx ${reference}
    ```

2. Create the FASTA sequence dictionary file

    ```bash
    reference="path/to/reference/UCSC_hg38/chr21.fa.gz"
    gatk CreateSequenceDictionary -R ${reference}
    ```

3. Testing GATK with CountReads to count the number of sequence reads in a BAM file

    ```bash
    gatk --java-options "-Xmx4G" CountReads \
        -R ${reference} \
        -I ${recal_reads} \
        --QUIET
    ```

### 3. HaplotypeCaller

```bash
output_vcf="path/to/variant_calling"

gatk --java-options "-Xmx4G" HaplotypeCaller \
    -R ${reference} \
    -I ${recal_reads} \
    -O "${output_vcf}/chr21_tumor_raw_variants.vcf" \
    -stand-call-conf 30
```

## VARIANT FILTERING (QUALITY-BASED FILTERING)

```mermaid
flowchart TD
    subgraph "VARIANT FILTERING"
    A("Raw SNPs + Indels"<br>VCF)--"<u>CNNScoreVariants<br>FilterVariantsTranches</u>"--> B{{"Analysis-Ready Variants<br>SNVs + Indels"}}
    style A fill:#34495e,stroke:#333,stroke-width:4px,color:#fff
    style B fill:#767076,stroke:#F9F2F4,stroke-width:2px,color:#fff
    end
```

**NOTE:**
In order to run [CNNSCoreVariants](#1-cnnscorevariants) and [FilterVariantsTranches](#2-filtervarianttranches) tools, one need to set up [Python Dependencies](https://github.com/broadinstitute/gatk#python) environment.

### 1. CNNScoreVariants

```bash
gatk --java-options "-Xmx4G" CNNScoreVariants \
    -R ${reference} \
    -V "${output_vcf}/chr21_tumor_raw_variants.vcf" \
    -O "${output_vcf}/chr21_tumor_cnn.vcf"
```

### 2. FilterVariantTranches

```bash
gatk --java-options "-Xmx4G" FilterVariantTranches \
    -V "${output_vcf}/chr21_tumor_cnn.vcf" \
    --resource "${known_snps}/dbsnp_146.hg38.vcf.gz" \
    --resource "${known_snps}/1000G_omni2.5.hg38.vcf.gz" \
    --resource "${known_snps}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
    --resource "${known_indels}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
    --resource "${known_indels}/Homo_sapiens_assembly38.known_indels.vcf.gz" \
    -O "${output_vcf}/chr21_tumor_filtered.vcf" \
    --snp-tranche 99.95 \
    --indel-tranche 99.4 \
    --info-key CNN_1D
```

#### Variant Genotype Filters

- Remove genotype DP < 30
- Removal of sites with low genotype quality score: GQ < 20

```bash
bcftools filter -i ' FORMAT/GQ>=20 & FORMAT/DP>=30' ${output_vcf}/chr21_tumor_filtered.vcf > ${output_vcf}/chr21_tumor_filteredGT.vcf
```

```mermaid
flowchart TD
    subgraph "VARIANT FILTERING"
    A("Raw SNPs + Indels"<br>VCF)--"Hard Filtering"--> B{{"Analysis-Ready Variants<br>SNVs + Indels"}}
    style A fill:#34495e,stroke:#333,stroke-width:4px,color:#fff
    style B fill:#767076,stroke:#F9F2F4,stroke-width:2px,color:#fff
    end
```

### 3. Hard Filtering

If one can not run both [CNNScoreVariants](#1-cnnscorevariants) and [FilterVariantTranches](#2-filtervarianttranches) tools, please follow  this [Hard Filtering](#3-hard-filtering) option.

**Note**:
Variant filtering using ***Hard filtering*** option will be run after variant calling step with HaplotypeCaller, which is the same as CNNScore filtering process.

If one has already run variant filtering with **CNNScoreVariant** option, there is no need to filter with **Hard Filtering** option anymore.

- #### Extract SNPs and INDELS

```bash
gatk SelectVariants \
    -R ${reference} \
    -V "${output_vcf}/chr21_tumor_raw_variants.vcf" \
    --select-type SNP \
    -O ${output_vcf}/chr21_tumor_raw_snps.vcf
```

```bash
gatk SelectVariants \
    -R ${reference} \
    -V "${output_vcf}/chr21_tumor_raw_variants.vcf" \
    --select-type INDEL \
    -O ${output_vcf}/chr21_tumor_raw_indels.vcf
```

- #### Filter SNPs

```bash
gatk VariantFiltration \
    -R ${reference} \
    -V "${output_vcf}/chr21_tumor_raw_snps.vcf" \
    -O "${output_vcf}/chr21_tumor_filtered_snps.vcf" \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 3.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-expression "DP < 30" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 20" \
    -genotype-filter-name "GQ_filter"
```

- #### Filter INDELs

```bash
gatk VariantFiltration \
    -R ${reference} \
    -V "${output_vcf}/chr21_tumor_raw_indels.vcf" \
    -O "${output_vcf}/chr21_tumor_filtered_indels.vcf" \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0" \
    -genotype-filter-expression "DP < 30" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 20" \
    -genotype-filter-name "GQ_filter"
```

- #### Select Variants that "PASS" filters

```bash
gatk SelectVariants \
    --exclude-filtered \
    -V "${output_vcf}/chr21_tumor_filtered_snps.vcf" \
    -O ${output_vcf}/chr21_tumor_ready_snps.vcf
```

```bash
gatk SelectVariants \
    --exclude-filtered \
    -V "${output_vcf}/chr21_tumor_filtered_indels.vcf" \
    -O ${output_vcf}/chr21_tumor_ready_indels.vcf
```

- #### Exclude variants that failed genotype filters

```bash
cat "${output_vcf}/chr21_tumor_ready_snps.vcf" | grep -v -E "DP_filter|GQ_filter" > ${output_vcf}/chr21_tumor_ready_snpsGT.vcf
```

```bash
cat "${output_vcf}/chr21_tumor_ready_indels.vcf" | grep -v -E "DP_filter|GQ_filter" > ${output_vcf}/chr21_tumor_ready_indelsGT.vcf
```

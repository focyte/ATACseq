# ATAC-seq
Pipeline for analysing ATAC-seq data from raw fastq files to called genomic peaks

---

## Table of Contents
1. [Overview](#overview)
2. [Pipeline Stages](#pipeline-stages)
    - [Set up Conda Environment](#set-up-conda-environment)
    - [Quality Control](#quality-control)
    - [Mapping](#mapping)
    - [Peak Calling](#peak-calling)
    - [Results](#results)
3. [How to Run](#how-to-run)
4. [Results](#results)

---

## Overview

## Pipeline Stages

## Set up Conda Environment

### Load the yml file
```console
conda env create -f atacseq.yml
```

## Quality Control

### Run FastQC
```console
fastqc -t 8 *fastq.gz -o ./fastqc_initial
```

### Run MultiQC
```console
cd ./fastqc_initial/
multiqc .
firefox multiqc_report.html
```

### Trim Reads
```console
cd ../trimming
./trim-reads-paired.sh
```

### Run FastQC on Trimmed Reads
```console
fastqc -t 8 *fq.gz -o fastqc_trimming/
cd ./fastqc_trimming/
multiqc .
firefox multiqc_report.html
```

## Mapping

### Index Reference Genome
Download index files from https://benlangmead.github.io/aws-indexes/bowtie and add them to the mapping folder.

```console
cd ../../../mapping
```

### Run Bowtie2 Alignment
```console
./bowtie2-mapping-paired.sh
```

### Convert SAM to BAM, Sort, and Index BAM
```console
./sam-to-bam.sh
```

## Peak Calling

### Mark Duplicates and Filter
```console
cd ../peak_call
./filter.sh
```

### Use MACS2 to Call Peaks
```console
./peak_call.sh
```


### Open IGV and View Peaks
```console
igv
```
In the IGV GUI choose file > Load from File

## How To Run

## Results

![Figure describing results generated using this pipeline](https://github.com/focyte/ATACseq/blob/main/ExampleFIgure.PNG)

#!/bin/bash

conda update -n base -c defaults conda

conda create -n atacseq

conda activate atacseq

conda install bioconda::bowtie2 bioconda::trim-galore bioconda::cutadapt bioconda::fastqc bioconda::multiqc bioconda::sra-tools bioconda::bedtools conda-forge::firefox bioconda::samtools bioconda::igv bioconda::picard

done
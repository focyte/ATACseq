#!/bin/bash

# Set the input directory
INPUT_DIR="../"

# Loop through each pair of fastq.gz files in the input directory
for fastq_file_1 in "$INPUT_DIR"/*_1.fastq.gz; do
    fastq_file_2="${fastq_file_1/_1.fastq.gz/_2.fastq.gz}"
    
    unique_id=$(basename "$fastq_file_1" | cut -d'_' -f1)

    trim_galore --paired --cores 2 "$fastq_file_1" "$fastq_file_2"

    if [ $? -eq 0 ]; then
        echo "Trimmed $unique_id"
    else
        echo "Error trimming $unique_id"
    fi
done

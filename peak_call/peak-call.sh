#!/bin/bash

# Process all sorted BAM files in the current directory ending with filtered.bam
for BAM_FILE in *filtered.bam; do
    # Skip if BAM_FILE is a directory or does not end with filtered.bam
    if [ ! -f "$BAM_FILE" ]; then
        continue
    fi

    # Extract the base name of the input file
    BASENAME=$(basename "$BAM_FILE")
    FILENAME="${BASENAME%.bam}"

    # Peak Calling with MACS2
    echo "Peak calling for $BAM_FILE..."

    # Run MACS2 with different output file names based on the input file name
    # macs2 callpeak -t "$BAM_FILE" -f BAMPE -n "${FILENAME}_ATACSeq_peak" -g 3100000000 --keep-dup al
    macs2 callpeak --broad -t "$BAM_FILE" -f BAMPE -n "${FILENAME}_ATACSeq_peak" -g 3100000000 --keep-dup all -B

    if [ $? -eq 0 ]; then
        echo "Peak calling for $BAM_FILE completed successfully."
    else
        echo "Error in peak calling for $BAM_FILE."
    fi
done

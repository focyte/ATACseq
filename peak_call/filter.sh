#!/bin/bash

# Process all sorted BAM files in the ../mapping directory
for BAM_FILE in ../mapping/*sorted.bam; do
    # Skip if BAM_FILE is a directory or does not end with sorted.bam
    if [ ! -f "$BAM_FILE" ]; then
        continue
    fi

    echo "Processing $BAM_FILE..."

    # Extract the file name without the directory and extension
    BAM_BASENAME=$(basename "$BAM_FILE")
    BAM_NAME="${BAM_BASENAME%.bam}"

    # Mark Duplicates
    picard MarkDuplicates \
          I="$BAM_FILE" \
          O="./${BAM_NAME}_marked_duplicates.bam" \
          M="./${BAM_NAME}_marked_dup_metrics.txt"

    if [ $? -eq 0 ]; then
        echo "Marking duplicates for $BAM_FILE completed successfully."
    else
        echo "Error in marking duplicates for $BAM_FILE."
        continue
    fi

    # Remove reads mapping to chromosomes "chrM" and "chrUn" and index the resulting BAM file
    samtools view -h "./${BAM_NAME}_marked_duplicates.bam" | awk '{if($3 != "chrM" && $3 != "chrUn"){print $0}}' | samtools view -Shb - > "./${BAM_NAME}_filtered.bam"

    if [ $? -eq 0 ]; then
        echo "Filtering for $BAM_FILE completed successfully."
    else
        echo "Error in filtering for $BAM_FILE."
        continue
    fi

    # Index Sorted BAM
    FILTERED_BAM_FILE="./${BAM_NAME}_filtered.bam"
    if [[ $FILTERED_BAM_FILE == *_filtered.bam ]]; then
        echo "Indexing sorted BAM for $FILTERED_BAM_FILE..."
        samtools index "$FILTERED_BAM_FILE"

        if [ $? -eq 0 ]; then
            echo "Indexing for $FILTERED_BAM_FILE completed successfully."
        else
            echo "Error in indexing sorted BAM for $FILTERED_BAM_FILE."
        fi
    fi
done

echo "All sorted BAM files processed."

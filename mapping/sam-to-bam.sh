#!/bin/bash

# Process all SAM files in the current directory
for SAM_FILE in *.sam; do
    # Skip if SAM_FILE is a directory
    if [ -d "$SAM_FILE" ]; then
        continue
    fi

    echo "Processing $SAM_FILE..."

    # Convert SAM to BAM
    echo "Converting $SAM_FILE to BAM..."
    samtools view -h -S -b -o "${SAM_FILE%.sam}.bam" "$SAM_FILE"

    if [ $? -eq 0 ]; then
        echo "Conversion to BAM for $SAM_FILE completed successfully."
    else
        echo "Error in converting $SAM_FILE to BAM."
        continue
    fi

    # Sort BAM
    echo "Sorting BAM for $SAM_FILE..."
    samtools sort "${SAM_FILE%.sam}.bam" -o "${SAM_FILE%.sam}_sorted.bam"

    if [ $? -eq 0 ]; then
        echo "Sorting for $SAM_FILE completed successfully."
    else
        echo "Error in sorting BAM for $SAM_FILE."
        continue
    fi

    # Index Sorted BAM
    echo "Indexing sorted BAM for $SAM_FILE..."
    samtools index "${SAM_FILE%.sam}_sorted.bam"

    if [ $? -eq 0 ]; then
        echo "Indexing for $SAM_FILE completed successfully."
    else
        echo "Error in indexing sorted BAM for $SAM_FILE."
    fi
done

echo "All SAM files processed."

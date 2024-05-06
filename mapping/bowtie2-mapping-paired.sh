#!/bin/bash

# Set the input directory
INPUT_DIR="../raw_data/trimming"

# Generate paired-end file names (file2)
# input files (file1) for correct pairing in the subsequent loop iterations
for file1 in "$INPUT_DIR"/*1_val_1.fq.gz; do
    file2="${file1/1_val_1/2_val_2}"

    echo "Processing raw paired-end files: $file1 and $file2"

    # Extract the sample name from the base filename of file1
    sample_name=$(basename "$file1" | cut -d '_' -f 1)

    # Run bowtie2 with the specified options.
    bowtie2 -x ./GRCh38_noalt_as -1 "$file1" -2 "$file2"  -p 8 -S "./${sample_name}.sam" -t
done

#!/bin/bash
#for multiply countries after trimming 


# Export the path to BWA
export PATH="/Users/eugenianikonorova/bwa:$PATH"

# Path to the reference genome index
REFERENCE_INDEX="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/DUW_3CX"

# Base directories
SRA_BASE_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/SRA"
MAPPED_BASE_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/mapped"

# Iterate over each country directory
for COUNTRY_DIR in $SRA_BASE_DIR/*; do
    # Ensure it is a directory
    if [ -d "$COUNTRY_DIR" ]; then
        COUNTRY_NAME=$(basename "$COUNTRY_DIR")

        # Create the mapped directory for the country if it doesn't exist
        MAPPED_COUNTRY_DIR="$MAPPED_BASE_DIR/$COUNTRY_NAME"
        mkdir -p "$MAPPED_COUNTRY_DIR"

        # Iterate over each fastq.gz file pair in the Trimmed subfolder
        for F in "$COUNTRY_DIR/Trimmed/"*_1_trimmed.fastq.gz; do
            I=${F%%_1_trimmed.fastq.gz}

            # Perform bwa mem alignment
            bwa mem -t 8 "$REFERENCE_INDEX" "${I}_1_trimmed.fastq.gz" "${I}_2_trimmed.fastq.gz" > "${I}.sam"
            echo "${I} mapping is finished"

            # Convert SAM to BAM
            samtools view -bS -F 0x4 "${I}.sam" -o "${I}.bam"
            echo "${I}.bam is received"

            # Sort BAM file
            samtools sort "${I}.bam" -o "${I}.sorted.bam"
            echo "${I}.bam is sorted"

            # Move the BAM files to the mapped directory
            mv "${I}.bam" "${MAPPED_COUNTRY_DIR}/"
            mv "${I}.sorted.bam" "${MAPPED_COUNTRY_DIR}/"
            echo "${I} was moved to $MAPPED_COUNTRY_DIR"

            # Clean up SAM file
            rm "${I}.sam"
            echo "${I}.sam was removed"
        done
    fi
done

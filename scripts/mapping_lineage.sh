#!/bin/bash

# Paths to the reference genome indexes
REFERENCE_INDEX_F_SW4="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/F_SW4"
REFERENCE_INDEX_L2_434_BU="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/L2_434_Bu"

# Base directories
SRA_BASE_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/test"
MAPPED_BASE_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/mapped"
LOG_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/logs"
RESULTS_FILE="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/mapping_results.txt"

# Path to the lineage text file
LINEAGE_FILE="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/mapping_data.txt"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Initialize results file with header
echo -e "Sample Name\tCountry\tLineage\tMapped\tReference" > "$RESULTS_FILE"

# Initialize the arrays
lineages=()
sample_names=()

# Read the lineage and sample name pairs into two arrays
while IFS=$'\t' read -r lineage sample_name; do
    # Skip the header line
    if [[ "$lineage" == "lineage" && "$sample_name" == "sample_name" ]]; then
        continue
    fi
    lineages+=("$lineage")
    sample_names+=("$sample_name")
done < <(sed 's/[[:space:]]*$//' "$LINEAGE_FILE")

# Function to find the lineage for a given sample name
get_lineage() {
    local sample=$1
    for i in "${!sample_names[@]}"; do
        if [[ "${sample_names[$i]}" == "$sample" ]]; then
            echo "${lineages[$i]}"
            return
        fi
    done
    echo "Unknown"
}

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
            SAMPLE_NAME=$(basename "$I")
            LOG_FILE="$LOG_DIR/${SAMPLE_NAME}_mapping.log"

            # Get the lineage of the sample using the get_lineage function
            LINEAGE=$(get_lineage "$SAMPLE_NAME")

            # Debugging output for lineage
            echo "Sample: $SAMPLE_NAME, Lineage found: $LINEAGE" >> "$LOG_FILE"

            # Determine which reference genome to use based on lineage
            if [[ "$LINEAGE" == "T1" || "$LINEAGE" == "T2" || "$LINEAGE" == "T1-T2" ]]; then
                REFERENCE_INDEX="$REFERENCE_INDEX_F_SW4"
                REFERENCE_NAME="F_SW4"
            elif [[ "$LINEAGE" == "LGV" || "$LINEAGE" == "L1" || "$LINEAGE" == "L2" || "$LINEAGE" == "L3" ]]; then
                REFERENCE_INDEX="$REFERENCE_INDEX_L2_434_BU"
                REFERENCE_NAME="L2_434_BU"
            else
                echo "Lineage for $SAMPLE_NAME not recognized, skipping..." >> "$LOG_FILE"
                continue
            fi

            # Perform bwa mem alignment and log the output
            bwa mem -t 8 "$REFERENCE_INDEX" "${I}_1_trimmed.fastq.gz" "${I}_2_trimmed.fastq.gz" > "${I}.sam" 2> "$LOG_FILE"

            # Check if mapping was successful by inspecting the exit status
            if [ $? -eq 0 ]; then
                MAPPED_STATUS="yes"
                echo "${SAMPLE_NAME} mapping is finished" >> "$LOG_FILE"

                # Convert SAM to BAM
                samtools view -bS -F 0x4 "${I}.sam" -o "${I}.bam" >> "$LOG_FILE" 2>&1
                echo "${SAMPLE_NAME}.bam is received" >> "$LOG_FILE"

                # Sort BAM file
                samtools sort "${I}.bam" -o "${I}.sorted.bam" >> "$LOG_FILE" 2>&1
                echo "${SAMPLE_NAME}.bam is sorted" >> "$LOG_FILE"

                # Move the sorted.BAM files to the mapped directory
                mv "${I}.sorted.bam" "$MAPPED_COUNTRY_DIR/"
                echo "${SAMPLE_NAME} was moved to $MAPPED_COUNTRY_DIR" >> "$LOG_FILE"
            else
                MAPPED_STATUS="no"
                echo "${SAMPLE_NAME} mapping failed" >> "$LOG_FILE"
            fi

            # Clean up SAM and BAM files
            rm "${I}.sam"
            rm "${I}.bam"
            echo "${SAMPLE_NAME}.bam and ${SAMPLE_NAME}.sam were removed" >> "$LOG_FILE"

            # Save the result to the results file
            echo -e "${SAMPLE_NAME}\t${COUNTRY_NAME}\t${LINEAGE}\t${MAPPED_STATUS}\t${REFERENCE_NAME}" >> "$RESULTS_FILE"

        done
    fi
done

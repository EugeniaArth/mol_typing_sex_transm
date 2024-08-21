#!/bin/bash
#It use two references for snp calling and  annotation  depends on the lineage associated with exact sample 
#(it use txt file mapping_data.txt with sample names and lineages). 
#If the sample name  matched with lineage T1, T2 or T1-T2 the .vcf it should uses reference genome F_SW4; 
#if the sample name matched with lineage LGV, L1, L2, or L3 -  reference genome L2_434_Bu. 
#There are also two databases for annovar - cpdb_F_SW4 and cpdb_L2_434_Bu.

# Update paths for BWA and GATK
export PATH="/Users/eugenianikonorova/gatk-4.3.0.0/:$PATH"

# Define reference files and directories
REFERENCE_F_SW4="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/F_SW4.fasta"
REF_ANNOTATION_F_SW4="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/F_SW4.gtf"
REFERENCE_L2_434_BU="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/L2_434_Bu.fasta"
REF_ANNOTATION_L2_434_BU="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/reference/L2_434_Bu.gtf"
ANNOVAR_DB_F_SW4="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/annovar_db/cpdb_F_SW4/"
ANNOVAR_DB_L2_434_BU="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/annovar_db/cpdb_L2_434_Bu/"
RESULTS_FILE="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/calling_results.txt"

# Base directories
MAPPED_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/test/"
CALLING_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/calling/"
LINEAGE_FILE="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/mapping_data.txt"
LOG_DIR="/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/calling/logs/"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Initialize results file with header
echo -e "Sample Name\tCountry\tLineage\tCalling\tReference" > "$RESULTS_FILE"

# Create sequence dictionary and index for both references
#samtools faidx $REFERENCE_F_SW4
#samtools faidx $REFERENCE_L2_434_BU
#java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar CreateSequenceDictionary R=$REFERENCE_F_SW4 O=${REFERENCE_F_SW4%.fasta}.dict
#java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar CreateSequenceDictionary R=$REFERENCE_L2_434_BU O=${REFERENCE_L2_434_BU%.fasta}.dict

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


# Process each country directory in the mapped directory
for COUNTRY_DIR in $MAPPED_DIR*/; do
    COUNTRY=$(basename "$COUNTRY_DIR")
    echo "Processing country: $COUNTRY"

    # Create necessary directories for this country
    BAM_DIR=${CALLING_DIR}/${COUNTRY}/bam
    VCF_DIR=${CALLING_DIR}/${COUNTRY}/vcf
    ANNOTATED_VCF_DIR=${CALLING_DIR}/${COUNTRY}/vcf_annotated
    mkdir -p $BAM_DIR $VCF_DIR $ANNOTATED_VCF_DIR

    # Process each sorted BAM file in the country's directory
    for BAM_FILE in ${COUNTRY_DIR}*.sorted.bam; do
        echo "Processing file: $BAM_FILE"

        # Extract base name from BAM file
        BASE_NAME=$(basename "$BAM_FILE" .sorted.bam)
        LOG_FILE="$LOG_DIR/${BASE_NAME}_calling.log"

        # Get the lineage of the sample using the get_lineage function
        SAMPLE_LINEAGE=$(get_lineage "$BASE_NAME")


        # Determine the reference genome and annotation database based on lineage
        if [[ "$SAMPLE_LINEAGE" == "T1" || "$SAMPLE_LINEAGE" == "T2" || "$SAMPLE_LINEAGE" == "T1-T2" ]]; then
            REFERENCE=$REFERENCE_F_SW4
            REF_ANNOTATION=$REF_ANNOTATION_F_SW4
            ANNOVAR_DB=$ANNOVAR_DB_F_SW4
            BUILDVER="F_SW4"
            REFERENCE_NAME="F_SW4"
        elif [[ "$SAMPLE_LINEAGE" == "LGV" || "$SAMPLE_LINEAGE" == "L1" || "$SAMPLE_LINEAGE" == "L2" || "$SAMPLE_LINEAGE" == "L3" ]]; then
            REFERENCE=$REFERENCE_L2_434_BU
            REF_ANNOTATION=$REF_ANNOTATION_L2_434_BU
            ANNOVAR_DB=$ANNOVAR_DB_L2_434_BU
            BUILDVER="L2_434_Bu"
            REFERENCE_NAME="L2_434_Bu"
        else
            echo "Lineage for $BASE_NAME not recognized, skipping..."
            continue
        fi

        # Add debug statements to verify correct reference assignment

          echo "Sample: $BASE_NAME"
          echo "Lineage: $SAMPLE_LINEAGE"
          echo "Selected Reference: $REFERENCE"
          echo "Reference Annotation: $REF_ANNOTATION"
          echo "ANNOVAR Database: $ANNOVAR_DB"
          echo "Build Version: $BUILDVER"

        # Ensure that the correct reference is being used
            echo "Using reference: $REFERENCE for sample: $BASE_NAME with lineage: $SAMPLE_LINEAGE"

        # Add read groups, sort, and index
        java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar AddOrReplaceReadGroups \
            -I $BAM_FILE \
            -O $BAM_DIR/${BASE_NAME}.rgs.bam \
            -RGID $BASE_NAME -RGLB miseq -RGSM $BASE_NAME -RGPL illumina -RGPU unit1
        samtools index $BAM_DIR/${BASE_NAME}.rgs.bam

        # Sort and index the BAM file with read groups
        samtools sort $BAM_DIR/${BASE_NAME}.rgs.bam -o $BAM_DIR/${BASE_NAME}.rgs.sorted.bam
        samtools index $BAM_DIR/${BASE_NAME}.rgs.sorted.bam

        # Mark duplicates and index
        java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar MarkDuplicates \
            -I $BAM_DIR/${BASE_NAME}.rgs.sorted.bam \
            -O $BAM_DIR/${BASE_NAME}.marked_dup.bam \
            -M $BAM_DIR/${BASE_NAME}.dup_metrics.txt
        samtools index $BAM_DIR/${BASE_NAME}.marked_dup.bam

        # Call variants using GATK HaplotypeCaller
        gatk HaplotypeCaller \
            -ploidy 1 \
            -R $REFERENCE \
            -I $BAM_DIR/${BASE_NAME}.marked_dup.bam  \
            -O $VCF_DIR/${BASE_NAME}.vcf

        # Filter variants based on quality metrics
        gatk --java-options "-Xmx4G" VariantFiltration -R $REFERENCE -V $VCF_DIR/${BASE_NAME}.vcf \
            --filter-expression 'QD<2.0 || FS>60.0 || SOR>3.0 || MQ<40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 ||  AF < 0.25  || DP <10' \
            --filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 --missing-values-evaluate-as-failing \
            -O $VCF_DIR/${BASE_NAME}.filtered.vcf

        echo filtration is performed for ${BASE_NAME}

     # убираем отфильтрованные / Exclude Filtered Variants
       gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V $VCF_DIR/${BASE_NAME}.filtered.vcf   -O $VCF_DIR/${BASE_NAME}.passed.vcf
       echo variants filtered for ${BASE_NAME}


    # Extract SNPs & Indels
    gatk SelectVariants \
            -R $REFERENCE \
            -V $VCF_DIR/${BASE_NAME}.passed.vcf \
            --select-type SNP \
            -O $VCF_DIR/${BASE_NAME}.snp.vcf
    echo SNP extracted V for ${BASE_NAME}

    gatk SelectVariants \
            -R $REFERENCE \
            -V $VCF_DIR/${BASE_NAME}.passed.vcf \
            --select-type INDEL \
            -O $VCF_DIR/${BASE_NAME}.indel.vcf
    echo INDEL extracted V for ${BASE_NAME}


    # name genes by simply intersecting with reference
    bedtools intersect -wb -a $VCF_DIR/${BASE_NAME}.snp.vcf -b $REF_ANNOTATION  -header >  $VCF_DIR/${BASE_NAME}.snp.named.vcf
    bedtools intersect -wb -a $VCF_DIR/${BASE_NAME}.indel.vcf -b $REF_ANNOTATION -header >  $VCF_DIR/${BASE_NAME}.indel.named.vcf
    echo bedtools intersect is finished for ${BASE_NAME}


    # Next step is annotation of variants using Annovar

    ## Need to change directory as annovar will create output in the working directory
    cd $VCF_DIR

    #Annotate SNPs and Predict Effects
    perl  ~/annovar/table_annovar.pl  $VCF_DIR/${BASE_NAME}.snp.named.vcf $ANNOVAR_DB --outfile $ANNOTATED_VCF_DIR/${BASE_NAME}.snp --protocol \
    refGene --operation g --buildver  $BUILDVER --vcfinput
    echo SNP annotation is performed for ${BASE_NAME}

    perl  ~/annovar/table_annovar.pl  $VCF_DIR/${BASE_NAME}.indel.named.vcf $ANNOVAR_DB --outfile $ANNOTATED_VCF_DIR/${BASE_NAME}.indel --protocol \
    refGene --operation g --buildver  $BUILDVER --vcfinput
    echo INDELS annotation is performed for ${BASE_NAME}


# After all processing steps, append the results to the calling_results.txt
    echo -e "${BASE_NAME}\t${COUNTRY}\t${SAMPLE_LINEAGE}\tcompleted\t${REFERENCE_NAME}" >> "$RESULTS_FILE"


      done
  done

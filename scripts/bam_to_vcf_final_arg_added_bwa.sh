#!bin/bash
# to run
# chmod +x bam_1st_vcf.sh
# начиная с sorted.bam файлов и до vcf

#Подготовка индексов референса к работе.
#samtools faidx $REF
#picard CreateSequenceDictionary R=$REF O=GCF_000008725.1_ASM872v1_genomic.dict

#запускаем в папке GATK4

export PATH="/Users/eugenianikonorova/bwa:$PATH"
export PATH="/Users/eugenianikonorova/gatk-4.3.0.0/:$PATH"

REFERENCE=~/Genomes_mol_typing/Chlamydiales/reference/GCF_000008725.1_ASM872v1_genomic.fna
#create bwa index for reference
cd ~/Genomes_mol_typing/Chlamydiales/reference/
#bwa index -p bwa_index $REFERENCE

REF_ANNOTATION=~/Genomes_mol_typing/Chlamydiales/reference/GCF_000008725.1_ASM872v1_genomic.gff

#create folders for resulted files
cd ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test
mkdir -p bam vcf vcf_annotated

for read1 in ~/Genomes_mol_typing/Chlamydiales/SRA/test/*_1.fastq.gz # for each sample
  do
    echo "working with file $read1"

    base=$(basename "$read1" _1.fastq.gz)
    echo "base name is $base"

    read1=~/Genomes_mol_typing/Chlamydiales/SRA/test/${base}_1.fastq.gz
    read2=~/Genomes_mol_typing/Chlamydiales/SRA/test/${base}_2.fastq.gz
    bwa_index=~/Genomes_mol_typing/Chlamydiales/reference/bwa_index
    sam=~/Genomes_mol_typing/Chlamydiales/GATK4/mapped/test/${base}.sam
    bam=~/Genomes_mol_typing/Chlamydiales/GATK4/mapped/test/${base}.bam
    sorted_bam=~/Genomes_mol_typing/Chlamydiales/GATK4/mapped/test/${base}.sorted.bam
    rgs_bam=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/bam/${base}.rgs.bam
    rgs_sorted_bam=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/bam/${base}.rgs.sorted.bam
    marked_dup_bam=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/bam/${base}.marked_dup.bam
    dup_metrics=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/bam/${base}.dup_metrics.txt
    variants=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.vcf
    variants_filtered=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.filtered.vcf
    variants_passed=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.passed.vcf
    snp_vcf=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.snp.vcf
    indel_vcf=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.indel.vcf
    snp_named_vcf=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.snp.named.vcf
    indel_named_vcf=~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf/${base}.indel.named.vcf
    annovar_db=~/Genomes_mol_typing/Chlamydiales/cpdb/


bwa mem -t 8 $bwa_index $read1 $read2 > $sam
echo  mapping is finished for ${base}
samtools view -bS -F 0x4 $sam -o $bam
echo ${base}.bam is recieved
samtools sort $bam -o $sorted_bam
echo ${i}.bam is sorted

echo "Mapping is finished"

#Next step is SNV calling using GATK4

echo "Running GATK..."
#  Add read groups, coordinate sort and index using AddOrReplaceReadGroups
  #добавляем read groups
  java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar \
AddOrReplaceReadGroups -I $sorted_bam -O $rgs_bam -RGID ${base} -RGLB miseq -RGSM ${base} -RGPL illumina -RGPU unit1

  samtools index $rgs_bam

echo AddOrReplaceReadGroups is performed for ${base}

  samtools sort $rgs_bam  -o $rgs_sorted_bam
  samtools index $rgs_sorted_bam


#mark duplicates
  java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar MarkDuplicates -I $rgs_sorted_bam \
        -O $marked_dup_bam \
        -M $dup_metrics

#samtools sort ${i}.marked_dup.bam  -o ${i}.marked_dup.sorted.bam
  samtools index $marked_dup_bam
echo  MarkDuplicatesis finished for ${base}

# удаляем лишнее сразу

#HaplotypeCaller
  gatk HaplotypeCaller -ploidy 1 -R $REFERENCE  -I $marked_dup_bam -O $variants
  echo gatk HaplotypeCaller is performed for ${base}

# Filter SNPs
  gatk --java-options "-Xmx4G" VariantFiltration -R $REFERENCE -V $variants \
 --filter-expression 'QD<2.0 || FS>60.0 || SOR>3.0 || MQ<40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 ||  AF < 0.25' \
 --filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 --missing-values-evaluate-as-failing \
 -O $variants_filtered

 echo filtration is performed for ${base}


 # убираем отфильтрованные / Exclude Filtered Variants
   gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V $variants_filtered   -O $variants_passed
   echo variants filtered for ${base}


# Extract SNPs & Indels
gatk SelectVariants \
        -R $REFERENCE \
        -V $variants_passed \
        --select-type SNP \
        -O $snp_vcf
echo SNP extracted V for ${base}

gatk SelectVariants \
        -R $REFERENCE \
        -V $variants_passed \
        --select-type INDEL \
        -O $indel_vcf
echo INDEL extracted V for ${base}


# name genes by simply intersecting with reference
bedtools intersect -wb -a $snp_vcf -b $REF_ANNOTATION  -header >  $snp_named_vcf
bedtools intersect -wb -a $indel_vcf -b $REF_ANNOTATION -header >  $indel_named_vcf
echo bedtools intersect is finished for ${base}


# Next step is annotation of variants using Annovar

## Need to change directory as annovar will create output in the working directory
cd ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/test/vcf_annotated/

#Annotate SNPs and Predict Effects
perl  ~/annovar/table_annovar.pl  $snp_named_vcf $annovar_db --outfile ${base}.snp --protocol refGene --operation g --buildver  Chl.p --vcfinput
echo SNP annotation is performed for ${base}

perl  ~/annovar/table_annovar.pl  $indel_named_vcf $annovar_db --outfile ${base}.indel --protocol refGene --operation g --buildver  Chl.p --vcfinput
echo INDELS annotation is performed for ${base}


# убираем лишние файлы
#rm -f ${i}.marked_dup.bam ${i}.snp.passed.vcf ${i}.indels.passed.vcf ${i}.SNP.vcf ${i}.INDEL.vcf \
  #${i}.snp.filtered.vcf ${i}.indels.filtered.vcf ${i}.recal_reads.bam ${i}.named.snp.vcf ${i}.named.indels.vcf \
#  ${i}.filtered_snps_final.vcf  ${i}.filtered_indels_final.vcf ${i}.raw_snps_recal.vcf ${i}.raw_indels_recal.vcf \
#  ${i}.raw_variants_recal.vcf ${i}.vcf
#rm -f *.idx *.bai *.orig *.log *.fa *.avinput *.txt

# Compile Statistics
#parse_metrics.sh sample_id > sample_id_report.csv

done

#!bin/bash
# to run
# chmod +x bam_1st_vcf.sh
# начиная с sorted.bam файлов и до vcf

#Подготовка индексов референса к работе.
#samtools faidx $REF
#picard CreateSequenceDictionary R=$REF O=GCF_000008725.1_ASM872v1_genomic.dict


for f in ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/Italy_noBQSR/*.sorted.bam # for each sample
  do
    i=${f%%.sorted.bam}

#  Add read groups, coordinate sort and index using AddOrReplaceReadGroups
  #добавляем read groups
java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar \
AddOrReplaceReadGroups -I ${i}.sorted.bam -O ${i}.rgs.bam -RGID ${i} -RGLB miseq -RGSM ${i} -RGPL illumina -RGPU unit1

samtools index ${i}.rgs.bam

echo AddOrReplaceReadGroups is performed for ${i}

samtools sort ${i}.rgs.bam  -o ${i}.rgs.sorted.bam
samtools index ${i}.rgs.sorted.bam


#mark duplicates
java -jar /Users/eugenianikonorova/gatk-4.3.0.0/picard.jar MarkDuplicates \
        -I ${i}.rgs.sorted.bam \
        -O ${i}.marked_dup.bam \
        -M ${i}.dup_metrics.txt

#samtools sort ${i}.marked_dup.bam  -o ${i}.marked_dup.sorted.bam
samtools index ${i}.marked_dup.bam
echo  MarkDuplicatesis finished for ${i}

# удаляем лишнее сразу
rm ${i}.rgs.bam
rm ${i}.rgs.bam.bai
rm ${i}.rgs.sorted.bam
rm ${i}.rgs.sorted.bam.bai

#HaplotypeCaller
  gatk HaplotypeCaller -ploidy 1 -R ~/Genomes_mol_typing/Chlamydiales/reference/GCF_000008725.1_ASM872v1_genomic.fna \
   -I ${i}.marked_dup.bam -O ${i}.vcf
  echo gatk HaplotypeCaller is performed for ${i}


# Extract SNPs & Indels
gatk SelectVariants \
        -R Reference.fna \
        -V ${i}.vcf \
        --select-type SNP \
        -O ${i}.SNP.vcf
echo SNP extracted V for ${i}

gatk SelectVariants \
        -R Reference.fna\
        -V ${i}.vcf \
        --select-type INDEL \
        -O ${i}.INDEL.vcf
echo INDEL extracted V for ${i}

# Filter SNPs
  gatk --java-options "-Xmx4G" VariantFiltration -R Reference.fna -V ${i}.SNP.vcf \
 --filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 || DP < 10 ' \
 --filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 --missing-values-evaluate-as-failing \
  -O ${i}.snp.filtered.vcf
  echo SNPs filtration №1 is performed for ${i}

#Filter Indels
gatk --java-options "-Xmx4G" VariantFiltration -R Reference.fna -V ${i}.INDEL.vcf \
--filter-expression 'QUAL < 30.0 || QD < 2.0 || FS > 60.0 || DP < 10 ' \
--filter-name lowQualFilter --cluster-window-size 10 --cluster-size 3 --missing-values-evaluate-as-failing \
        -O ${i}.indels.filtered.vcf
echo INDELs filtration №1 is performed for ${i}

# убираем отфильтрованные / Exclude Filtered Variants
  gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V ${i}.snp.filtered.vcf  -O ${i}.snp.passed.vcf
  echo SNPs filtered for ${i}
  gatk --java-options "-Xmx4G" SelectVariants --exclude-filtered -V ${i}.indels.filtered.vcf  -O ${i}.indels.passed.vcf
  echo INDELs filtered for ${i}




# пересекаем с референсом
bedtools intersect -wb -a ${i}.snp.passed.vcf   \
-b ~/Genomes_mol_typing/Chlamydiales/reference/GCF_000008725.1_ASM872v1_genomic.gff  \
-header >  ${i}.named.snp.vcf
bedtools intersect -wb -a ${i}.indels.passed.vcf \
-b ~/Genomes_mol_typing/Chlamydiales/reference/GCF_000008725.1_ASM872v1_genomic.gff  \
-header >  ${i}.named.indels.vcf
echo bedtools intersect is finished for ${i}



#Annotate SNPs and Predict Effects
perl  ~/annovar/table_annovar.pl  ${i}.named.snp.vcf ~/Genomes_mol_typing/Chlamydiales/cpdb/ --outfile ${i}.snp --protocol refGene --operation g --buildver  Chl.p --vcfinput
echo SNP annotation is performed for ${i}

perl  ~/annovar/table_annovar.pl  ${i}.named.indels.vcf ~/Genomes_mol_typing/Chlamydiales/cpdb/ --outfile ${i}.indel --protocol refGene --operation g --buildver  Chl.p --vcfinput
echo INDELS annotation is performed for ${i}

#перемещаем итог в другую папку annovar
mv ${i}.snp.Chl.p_multianno.vcf ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/Italy_noBQSR
mv ${i}.indel.Chl.p_multianno.vcf ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/Italy_noBQSR

# убираем лишние файлы
rm -f ${i}.marked_dup.bam ${i}.snp.passed.vcf ${i}.indels.passed.vcf ${i}.SNP.vcf ${i}.INDEL.vcf \
  ${i}.snp.filtered.vcf ${i}.indels.filtered.vcf ${i}.recal_reads.bam ${i}.named.snp.vcf ${i}.named.indels.vcf \
  ${i}.filtered_snps_final.vcf  ${i}.filtered_indels_final.vcf ${i}.raw_snps_recal.vcf ${i}.raw_indels_recal.vcf \
  ${i}.raw_variants_recal.vcf ${i}.vcf
rm -f *.idx *.bai *.orig *.log *.fa *.avinput *.txt

# Compile Statistics
#parse_metrics.sh sample_id > sample_id_report.csv

done

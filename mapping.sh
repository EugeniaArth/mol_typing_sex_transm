#!bin/bash

#bowtie2.sh

for f in ~/Genomes_mol_typing/Chlamydiales/SRA/Fastq/Italy/*_1.fastq.gz # for each sample
do
  i=${f%%_1.fastq.gz}
  bowtie2 -p 8  -1 ${i}_1.fastq.gz  -2 ${i}_2.fastq.gz -x ~/Genomes_mol_typing/Chlamydiales/reference/Chl.trach_index -S ${i}.sam
echo ${i} mapping is finished
  samtools view -bS -F 0x4 ${i}.sam -o ${i}.bam
echo ${i}.bam is recieved
  samtools sort ${i}.bam -o ${i}.sorted.bam
echo ${i}.bam is sorted
#  mv ${i}.sorted.bam ~/Genomes_mol_typing/Chlamydiales/GATK4/mapped ~/Genomes_mol_typing/Chlamydiales/GATK4/calling/Italy
echo ${i} was moved to GATK4/calling folder
  #rm ${i}.bam
  #rm ${i}.sam
  echo ${i} .bam and .sam were removed
done

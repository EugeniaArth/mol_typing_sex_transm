# mol_typing_sex_transm
Scripts and some data on search for the project on SNP search in N.gonorrhoeae, T.pallidum, Ch.trachomatis

**Data preparation**

Download data according available SRA using prefetch and fasterq-dump (fasterq-dump.sh script). Quality control was performed using fastqc and multiqc. Fastq trimming was performed using trimmomatic (trim_fastq.sh) with parameters: PE -threads 8 -phred33 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10  LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:50



# mol_typing_sex_transm
Scripts and some data on search for the project on SNP search in N.gonorrhoeae, T.pallidum, Ch.trachomatis

**Data preparation**

Download data according available SRA using prefetch and fasterq-dump (fasterq-dump.sh script). Quality control was performed using fastqc and multiqc. Fastq trimming was performed using trimmomatic (trim_fastq.sh) with parameters: PE -threads 8 -phred33 ILLUMINACLIP:Sequencing_adaptors.fasta:2:30:10  LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30 MINLEN:50


**Result summary**

1) First, use the script grab_exonic_variant_function_multi.py – it takes the data from files {ID}.....exonic_variant_function and create summary table of every SNP found in sample with Start, End, Gene name, Protein name, Type of SNP, Position,	REF, and	ALT. The resulting files are {ID}.snp_processed_data.csv
2) Then use script named ratio_snp_%.py – it collects data from {ID}.snp_processed_data.csv and create two tables: 1) summary_table.scv where the data about all SNP in all genes of the all samples in folder are collected (contain fields Sample name,	Gene name, Protein name, Start,	End, Ratio,	nonsynSNP, and synSNP) and filtered based on Ratio with the uniqe line for every sample and gene with SNP; 2) stats_table.scv where the data are concateneted again but now the genes are uniqe with number of samples which have SNPs in this gene, mean Ratio +/- SD added.


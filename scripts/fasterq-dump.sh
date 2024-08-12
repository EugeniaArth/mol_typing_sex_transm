#!bin/bash
# to run
# chmod +x fasterq-dump.sh
date  ## echo the date at start

export PATH=$HOME/sratoolkit.3.0.7-mac64/bin:$PATH

cat ~/Genomes_mol_typing/Chlamydiales/sra_id_cleaned.txt | parallel "prefetch  -O ~/Genomes_mol_typing/Chlamydiales/SRA -p  {}"

# Loop through each directory containing .sra files
for dir in */; do
    # Check if the directory contains .sra files
    if compgen -G "${dir}/*.sra" > /dev/null; then
        # Extract the directory name
        dirname=$(basename "$dir")

        # Loop through each .sra file in the directory
        for sra_file in "${dir}"/*.sra; do
            # Use fasterq-dump to convert .sra to .fastq.gz and save in the Fastq directory
            fasterq-dump --split-files -p -O Fastq "$sra_file"

echo fasterq-dump is performed for "$sra_file"
            # Rename the resulting FASTQ files based on the folder name gzip
           gzip "~/Genomes_mol_typing/Chlamydiales/SRA/${dirname}_1.fastq"
           gzip "~/Genomes_mol_typing/Chlamydiales/SRA/${dirname}_2.fastq"

echo unzipped "$sra_file" is gzipped

rm -r "${dir}"
echo ${dir} is removed


        done
    fi
done

echo "FASTQ files have been saved in the Fastq directory."

date  ## echo the date at start

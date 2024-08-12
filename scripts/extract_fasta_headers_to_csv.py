

import csv

def extract_fasta_headers_to_csv(fasta_file, csv_file):
    # Open the FASTA file for reading
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    # List to hold extracted header data
    header_data = []

    # Process each line in the FASTA file
    for line in lines:
        if line.startswith('>'):  # Identify header lines
            # Remove the initial '>' and split the header by '|'
            parts = line[1:].strip().split('|')
            header_data.append(parts)

    # Define the CSV column headers
    csv_headers = ['ID', 'E_Ar', 'Organism', 'Empty1', 'E', 'Country', 'Project', 'Sample', 'Empty2', 'Run']

    # Write the extracted header data to a CSV file
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_headers)  # Write the header row
        writer.writerows(header_data)  # Write the data rows

# Specify the input FASTA file and output CSV file paths
fasta_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/BIGSdb_4164854_8720967734_79623.fas'

csv_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/test.csv'

# Extract headers and write to CSV
extract_fasta_headers_to_csv(fasta_file, csv_file)

print(f"Data extracted to {csv_file} successfully.")

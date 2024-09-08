
import csv
import os

def read_mapping_file(mapping_file):
    mapping_dict = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            lineage, sample_name = line.strip().split()
            mapping_dict[sample_name] = lineage
    return mapping_dict

import pandas as pd  # Added pandas for reading CSV

def read_country_file(country_file_path):
    # Updated to read the CSV file using pandas with semicolon delimiter
    country_dict = {}
    df = pd.read_csv(country_file_path, sep=';')
    for index, row in df.iterrows():
        country_dict[row['run_accession']] = row['country']
    return country_dict

def process_files(input_directory, output_file, mapping_file, country_file):
    # Read the mapping data
    mapping_dict = read_mapping_file('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/data/mapping_data.txt')

    # Read the country data
    country_dict = read_country_file('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/data/country_data.csv')

    # Your existing processing code here...
    # For example, assume you are reading input files and generating a summary table

    # Reading the original summary table (if it exists)
    if os.path.exists(output_file):
        with open(output_file, 'r') as infile:
            reader = csv.reader(infile)
            headers = next(reader)  # Read the header
            rows = [row for row in reader]
    else:
        headers = ["Sample name", "Gene name", "Protein name", 'Protein ID', "Start", "End", "Ratio", "nonsynSNP", "synSNP"]  # Example headers, replace with your actual headers
        rows = []

    # Adding the "Country" column
    headers.insert(1, 'Country')  # Insert 'Country' after the 'Sample Name'

    # Adding the corresponding country to each row
    for row in rows:
        sample_name = row[0]  # Assuming the sample name is the first column
        country = country_dict.get(sample_name, 'Unknown')  # Get the country or default to 'Unknown'
        row.insert(1, country)

    # Writing the updated table back to the CSV
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(headers)
        writer.writerows(rows)

    print(f"Updated summary table saved to {output_file}")

# Example usage of the script
input_directory = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/'  # Replace with the actual directory
output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_F_SW4.csv'
#output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_L2_434_Bu.csv'  # Replace with the actual output file
mapping_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/data/mapping_data.txt'  # Replace with the actual mapping file
country_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/data/country_data.csv'  # Replace with the actual country file

process_files(input_directory, output_file, mapping_file, country_file)

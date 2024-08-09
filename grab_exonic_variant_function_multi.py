import pandas as pd
import re
import os

def extract_gene_name(info):
    match = re.search(r"Name=([^;]+)", info)
    return match.group(1) if match else None

def extract_protein_name(info):
    match = re.search(r"product=([^;]+)", info)
    return match.group(1) if match else None

def process_exonic_variant_function(file_path, output_dir):
    # Extract the file prefix
    filename = os.path.basename(file_path)
    prefix_parts = filename.split('.')
    prefix = '.'.join(prefix_parts[:2])

    # Load the data from the provided file
    data = pd.read_csv(file_path, sep='\t', header=None)

    # Initialize an empty list to hold the processed data
    processed_data = []

    # Iterate over the rows and extract necessary information
    for i in range(len(data)):
        # Check if the current row contains 'gene' in column 29
        if 'ID=gene-' in str(data.iloc[i, 29]):
            start = data.iloc[i, 24]
            end = data.iloc[i, 25]
            gene_info = data.iloc[i, 29]
            gene_name = extract_gene_name(gene_info)

            # Find the corresponding CDS line
            if i + 1 < len(data) and 'ID=cds-' in str(data.iloc[i + 1, 29]):
                protein_info = data.iloc[i + 1, 29]
                protein_name = extract_protein_name(protein_info)

                # Extract SNP data from the relevant columns
                snp_type = data.iloc[i, 1]
                position = data.iloc[i, 5]
                ref = data.iloc[i, 6]
                alt = data.iloc[i, 7]

                # Append the extracted data as a new row
                processed_data.append([start, end, gene_name, protein_name, snp_type, position, ref, alt])

    # Create a DataFrame from the processed data
    columns = ["Start", "End", "Gene name", "Protein name", "Type of SNP", "Position", "REF", "ALT"]
    result_df = pd.DataFrame(processed_data, columns=columns)

    # Save the result to a new CSV file with the prefix
    output_file = os.path.join(output_dir, f"{prefix}_processed_data.csv")
    result_df.to_csv(output_file, index=False)

    print(f"Processing complete for {file_path}. The result is saved to {output_file}.")

# Example usage
input_dir = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/calling/Italy'  # Replace with the actual path to your directory
output_dir = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/Italy'  # Replace with the actual path to your output directory
os.makedirs(output_dir, exist_ok=True)

# Process all files in the input directory
for file in os.listdir(input_dir):
    if file.endswith('.exonic_variant_function'):
        file_path = os.path.join(input_dir, file)
        process_exonic_variant_function(file_path, output_dir)

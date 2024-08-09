import pandas as pd
import os

def extract_gene_name(info):
    match = re.search(r"Name=([^;]+)", info)
    return match.group(1) if match else None

def extract_protein_name(info):
    match = re.search(r"product=([^;]+)", info)
    return match.group(1) if match else None

def calculate_gene_length(start, end):
    return abs(int(end) - int(start))

def process_summary(input_dir, summary_output_file, stats_output_file):
    # Initialize an empty list to hold the summary data
    summary_data = []
    ratio_data = []

    # Process each file in the input directory
    for file in os.listdir(input_dir):
        if file.endswith('.snp_processed_data.csv'):
            # Extract sample name from the file name
            sample_name = file.split('.')[0]

            # Load the data from the file
            file_path = os.path.join(input_dir, file)
            data = pd.read_csv(file_path)

            # Group by gene name, protein name, and start and end coordinates, and calculate the counts of nonsynonymous and synonymous SNVs
            grouped = data.groupby(['Gene name', 'Protein name', 'Start', 'End', 'Type of SNP']).size().unstack(fill_value=0)

            for (gene_name, protein_name, start, end), counts in grouped.iterrows():
                nonsynonymous_count = counts.get('nonsynonymous SNV', 0)
                synonymous_count = counts.get('synonymous SNV', 0)
                gene_length = calculate_gene_length(start, end)
                n_of_nonsynonSNP_to_100bp = round((nonsynonymous_count / gene_length) * 100, 2) if gene_length > 0 else 0
                n_of_synonSNP_to_100bp = round((synonymous_count / gene_length) * 100, 2) if gene_length > 0 else 0

                if nonsynonymous_count == 0 and synonymous_count == 0:
                    ratio = "not applicable, single SNP"
                elif synonymous_count == 0:
                    ratio = "not applicable, single SNP"
                else:
                    ratio = nonsynonymous_count / synonymous_count

                # Append the extracted data as a new row
                summary_data.append([sample_name, gene_name, protein_name, start, end, gene_length, nonsynonymous_count, synonymous_count, n_of_nonsynonSNP_to_100bp, n_of_synonSNP_to_100bp, ratio])
                if isinstance(ratio, (int, float)):
                    ratio_data.append([sample_name, gene_name, protein_name, start, end, ratio, nonsynonymous_count, synonymous_count, gene_length, n_of_nonsynonSNP_to_100bp, n_of_synonSNP_to_100bp])

    # Create a DataFrame from the summary data
    summary_columns = ["Sample name", "Gene name", "Protein name", "Start", "End", "Gene length", "Number_of_nonsynonSNP", "Number_of_synonSNP", "N_of_nonsynonSNP_to_100bp", "N_of_synonSNP_to_100bp", "Ratio"]
    summary_df = pd.DataFrame(summary_data, columns=summary_columns)

    # Convert ratios to numeric where possible for sorting, handle non-numeric ratios separately
    summary_df['Numeric Ratio'] = pd.to_numeric(summary_df['Ratio'], errors='coerce')
    summary_df = summary_df.sort_values(by='Numeric Ratio', ascending=False, na_position='last').drop(columns=['Numeric Ratio'])

    # Save the summary DataFrame to a CSV file
    summary_df.to_csv(summary_output_file, index=False)
    print(f"Summary table has been saved to {summary_output_file}.")

    # Create a DataFrame from the ratio data
    ratio_columns = ["Sample name", "Gene name", "Protein name", "Start", "End", "Ratio", "Total nonsynonymous SNV", "Total synonymous SNV", "Gene length", "N_of_nonsynonSNP_to_100bp", "N_of_synonSNP_to_100bp"]
    ratio_df = pd.DataFrame(ratio_data, columns=ratio_columns)

    # Calculate statistics for each gene and protein combination
    stats_df = ratio_df.groupby(["Gene name", "Protein name", "Start", "End"]).agg(
        Gene_length=('Gene length', 'mean'),
        Sample_Count=('Ratio', 'count'),
        Total_nonsynonymous_SNV=('Total nonsynonymous SNV', 'sum'),
        Total_synonymous_SNV=('Total synonymous SNV', 'sum'),
        Mean_N_of_nonsynonSNP_to_100bp=('N_of_nonsynonSNP_to_100bp', lambda x: round(x.mean(), 2)),
        SD_of_nonsynonSNP_to_100bp=('N_of_nonsynonSNP_to_100bp', lambda x: round(x.std(), 2)),
        Mean_N_of_synonSNP_to_100bp=('N_of_synonSNP_to_100bp', lambda x: round(x.mean(), 2)),
        SD_of_synonSNP_to_100bp=('N_of_synonSNP_to_100bp', lambda x: round(x.std(), 2))
    ).reset_index()

    # Sort by Mean_N_of_nonsynonSNP_to_100bp column
    stats_df = stats_df.sort_values(by='Mean_N_of_nonsynonSNP_to_100bp', ascending=False)

    # Save the stats DataFrame to a CSV file
    stats_df.to_csv(stats_output_file, index=False)
    print(f"Statistics table has been saved to {stats_output_file}.")



# Example usage
input_dir = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/Italy_noBQSR'  # Replace with the actual path to your directory containing the processed files
summary_output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/Italy_noBQSR/summary_table.csv'  # Replace with the desired summary output file path
stats_output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/Italy_noBQSR/stats_table.csv'  # Replace with the desired statistics output file path

process_summary(input_dir, summary_output_file, stats_output_file)

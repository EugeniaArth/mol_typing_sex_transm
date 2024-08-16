import pandas as pd
import os

def extract_gene_name(info):
    match = re.search(r"Name=([^;]+)", info)
    return match.group(1) if match else None

def extract_protein_name(info):
    match = re.search(r"product=([^;]+)", info)
    return match.group(1) if match else None

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

                if nonsynonymous_count == 0 and synonymous_count == 0:
                    ratio = "not applicable, single SNP"
                elif synonymous_count == 0:
                    ratio = "not applicable, single SNP"
                else:
                    ratio = nonsynonymous_count / synonymous_count

                # Append the extracted data as a new row
                summary_data.append([sample_name, gene_name, protein_name, start, end, ratio, nonsynonymous_count, synonymous_count])
                if isinstance(ratio, (int, float)):
                    ratio_data.append([sample_name, gene_name, protein_name, start, end, ratio, nonsynonymous_count, synonymous_count])

    # Create a DataFrame from the summary data
    summary_columns = ["Sample name", "Gene name", "Protein name", "Start", "End", "Ratio", "nonsynSNP", "synSNP"]
    summary_df = pd.DataFrame(summary_data, columns=summary_columns)

    # Convert ratios to numeric where possible for sorting, handle non-numeric ratios separately
    summary_df['Numeric ratio'] = pd.to_numeric(summary_df['Ratio'], errors='coerce')
    summary_df = summary_df.sort_values(by='Numeric ratio', ascending=False, na_position='last').drop(columns=['Numeric ratio'])

    # Save the summary DataFrame to a CSV file
    summary_df.to_csv(summary_output_file, index=False)
    print(f"Summary table has been saved to {summary_output_file}.")

    # Create a DataFrame from the ratio data
    ratio_columns = ["Sample name", "Gene name", "Protein name", "Start", "End", "Ratio", "Total nonsynSNP", "Total synSNP"]
    ratio_df = pd.DataFrame(ratio_data, columns=ratio_columns)

    # Calculate statistics for each gene and protein combination
    stats_df = ratio_df.groupby(["Gene name", "Protein name", "Start", "End"]).agg(
        Mean_Ratio=('Ratio', lambda x: round(x.mean(), 2)),
        SD=('Ratio', lambda x: round(x.std(), 2)),
        Sample_Count=('Ratio', 'count'),
        Total_nonsynSNP=('Total nonsynSNP', 'sum'),
        Total_synSNP=('Total synSNP', 'sum')
    ).reset_index()

    # Calculate the percentage of samples where the ratio was calculated
    total_samples = len(ratio_df['Sample name'].unique())
    stats_df['% of Samples with calculated Ratio'] = round((stats_df['Sample_Count'] / total_samples) * 100, 2)

    # Calculate the percentage of samples with ratio >= mean
    ratio_df = ratio_df.merge(stats_df[['Gene name', 'Protein name', 'Start', 'End', 'Mean_Ratio']], on=['Gene name', 'Protein name', 'Start', 'End'])
    ratio_df['Ratio >= Mean_Ratio'] = ratio_df['Ratio'] >= ratio_df['Mean_Ratio']
    percentage_df = ratio_df.groupby(["Gene name", "Protein name", "Start", "End"])['Ratio >= Mean_Ratio'].mean().reset_index()
    percentage_df['Ratio >= Mean_Ratio'] = round(percentage_df['Ratio >= Mean_Ratio'] * 100, 2)

    stats_df = stats_df.merge(percentage_df[['Gene name', 'Protein name', 'Start', 'End', 'Ratio >= Mean_Ratio']], on=['Gene name', 'Protein name', 'Start', 'End'])

    # Sort by Mean_Ratio column
    stats_df = stats_df.sort_values(by='Mean_Ratio', ascending=False)

    # Save the stats DataFrame to a CSV file
    stats_df.to_csv(stats_output_file, index=False)
    print(f"Statistics table has been saved to {stats_output_file}.")

# Example usage
input_dir = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/test'  # Replace with the actual path to your directory containing the processed files
summary_output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/test/summary_table.csv'  # Replace with the desired summary output file path
stats_output_file = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/GATK4/analysis/test/stats_table.csv'  # Replace with the desired statistics output file path

process_summary(input_dir, summary_output_file, stats_output_file)

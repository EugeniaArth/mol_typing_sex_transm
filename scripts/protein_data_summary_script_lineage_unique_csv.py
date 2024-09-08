
import pandas as pd
import warnings
import os

# Suppress specific warnings
warnings.filterwarnings("ignore", category=UserWarning)

# Load the CSV file instead of the Excel file
input_csv_path = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/F_SW4_L2_434_Bu.csv'
df = pd.read_csv(input_csv_path)


# Drop duplicates based on 'Genbank.protein.name' and 'Sample.name' to ensure unique samples
df_unique = df.drop_duplicates(subset=['Genbank.protein.name', 'Sample.name']).copy()

# Use .loc to ensure we're modifying the DataFrame safely
df_unique.loc[:, 'Ratio'] = pd.to_numeric(df_unique['Ratio'], errors='coerce')


# Group by 'Genbank.protein.name' and calculate the required statistics
summary_df = df_unique.groupby('Genbank.protein.name').agg({
    'Orthogroup', # Copy other data
    'Gene.name',
    'Protein.name',
    'Protein.ID',
    'Locus.tag',
    'Start',
    'End',
    'Ratio': ['mean', 'std'],  # Mean and standard deviation of Ratio (ignoring NaNs)
    'nonsynSNP': 'sum',  # Sum of nonsynSNP
    'synSNP': 'sum',     # Sum of synSNP
    'Sample.name': 'count',  # Total number of samples
})

# Flatten the MultiIndex columns created by aggregation
summary_df.columns = ['_'.join(col).strip() for col in summary_df.columns.values]

# Save the processed data to a new CSV file
output_csv_path = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/processed_F_SW4_L2_434_Bu.csv'
summary_df.to_csv(output_csv_path, index=True)

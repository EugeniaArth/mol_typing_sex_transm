import pandas as pd

# Load the SNP tables and gene info files
df_F_SW4 = pd.read_csv('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_F_SW4.csv')
gene_info_F_SW4 = pd.read_csv('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/gene_info_F_SW4.csv', delimiter=';')

# Step 1: Create the "Genbank protein name" column in the summary tables
df_F_SW4['Genbank protein name'] = None

# Step 2: Update "Genbank protein name" based on specific strain names

# For F_SW4
for index, row in df_F_SW4.iterrows():
    if row['Gene name'].startswith('FSW4_'):
        match = gene_info_F_SW4[gene_info_F_SW4['locus_tag'] == row['Gene name']]
        if not match.empty:
            df_F_SW4.at[index, 'Genbank protein name'] = match['product'].values[0]

# Step 3: Create "Genbank protein name 2" based on common gene names
df_F_SW4['Genbank protein name 2'] = None


# For F_SW4
for index, row in df_F_SW4.iterrows():
    match = gene_info_F_SW4[gene_info_F_SW4['gene_name'] == row['Gene name']]
    if not match.empty:
        df_F_SW4.at[index, 'Genbank protein name 2'] = match['product'].values[0]

# Step 4: Copy data from "Genbank protein name 2" to "Genbank protein name" where the latter is empty
df_F_SW4['Genbank protein name'] = df_F_SW4['Genbank protein name'].combine_first(df_F_SW4['Genbank protein name 2'])

# Step 5: Drop the "Genbank protein name 2" column
df_F_SW4.drop(columns=['Genbank protein name 2'], inplace=True)


# Save the updated SNP tables to new CSV files
df_F_SW4.to_csv('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/F_SW4_final_snp_table.csv', index=False)

# Debug: Check the final result before saving
print("Final table head:\n", df_F_SW4.head())

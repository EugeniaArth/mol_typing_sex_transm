import pandas as pd

# Step 1: Load the Excel files
gene_info_F_SW4 = pd.read_excel('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/gene_info_F_SW4.xlsx')
gene_info_L2_434_Bu = pd.read_excel('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/gene_info_L2_434_Bu.xlsx')

# Step 2: Prepare the data (Convert relevant columns to strings and clean up any extra spaces)
gene_info_F_SW4['locus_tag'] = gene_info_F_SW4['locus_tag'].astype(str).str.strip()
gene_info_F_SW4['old_locus_tag'] = gene_info_F_SW4['old_locus_tag'].astype(str).str.strip()

gene_info_L2_434_Bu['locus_tag'] = gene_info_L2_434_Bu['locus_tag'].astype(str).str.strip()
gene_info_L2_434_Bu['old_locus_tag'] = gene_info_L2_434_Bu['old_locus_tag'].astype(str).str.strip()

# Load the SNP table
snp_table = pd.read_csv('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_L2_434_Bu.csv', delimiter=',')

# Convert 'Locus tag' column to string and clean up any extra spaces
snp_table['Locus tag'] = snp_table['Locus tag'].astype(str).str.strip()

# Ensure the "Genbank protein name" column exists in the SNP table
if 'Genbank protein name' not in snp_table.columns:
    snp_table['Genbank protein name'] = None

# Merge using both 'locus_tag' and 'old_locus_tag'
merged_table = snp_table.merge(gene_info_L2_434_Bu[['locus_tag', 'product']],
                               left_on='Locus tag',
                               right_on='locus_tag',
                               how='left')

# Fill missing values by merging on the 'old_locus_tag'
merged_table = merged_table.merge(gene_info_L2_434_Bu[['old_locus_tag', 'product']],
                                  left_on='Locus tag',
                                  right_on='old_locus_tag',
                                  how='left',
                                  suffixes=('', '_old'))

# Combine the 'Genbank protein name' with the product names from both merges
merged_table['Genbank protein name'] = merged_table['Genbank protein name'].combine_first(merged_table['product'])
merged_table['Genbank protein name'] = merged_table['Genbank protein name'].combine_first(merged_table['product_old'])

# Drop the extra columns used for merging
merged_table.drop(columns=['locus_tag', 'product', 'old_locus_tag', 'product_old'], inplace=True)

# Save the updated SNP table with Genbank protein names filled
merged_table.to_csv('/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/L2_434_Bu_final_snp_table.csv', index=False)

# Debug: Check the final result before saving
print("Final table head:\n", merged_table.head())

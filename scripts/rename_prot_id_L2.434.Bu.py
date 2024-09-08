import pandas as pd

# Load CSV file
new_summary_table_path = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_L2_434_Bu_with_Orthogroup.csv'
new_summary_df = pd.read_csv(new_summary_table_path)

# Load the gene info table (Excel file)
gene_info_path = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/gene_info_L2_434_Bu.xlsx'
gene_info_df = pd.read_excel(gene_info_path)

# Debugging: Print out the columns of both dataframes to verify they exist
print("New Summary Table Columns:", new_summary_df.columns)
print("Gene Info Table Columns:", gene_info_df.columns)

# Perform the merge between the new summary table and the gene info table
merged_df = pd.merge(new_summary_df, gene_info_df[['old_locus_tag', 'protein_id']], how='left', left_on='Locus tag', right_on='old_locus_tag')

# Debugging: Verify the merge result and the columns in the merged dataframe
print("Merged DataFrame Columns:", merged_df.columns)
print("Sample of Merged DataFrame:", merged_df.head())

# Check if 'protein_id' exists in the merged DataFrame
if 'protein_id' in merged_df.columns:
    # Replace the "Protein ID" in the new summary table with the "protein_id" from the gene info table
    merged_df['Protein ID'] = merged_df['protein_id']
    print("Protein ID column successfully updated.")
else:
    print("Error: 'protein_id' column not found in the merged DataFrame.")

# Drop the extra columns used for merging
final_df = merged_df.drop(columns=['old_locus_tag', 'protein_id'], errors='ignore')

# Save the updated dataframe to a new CSV file
output_path = '/Users/eugenianikonorova/Genomes_mol_typing/Chlamydiales/snp/summary_table_L2_434_Bu_with_Orthogroup.csv'
final_df.to_csv(output_path, index=False)

print(f"Updated summary table saved to {output_path}")

# Additional analysis
# For example, you can perform summary statistics or any other analysis on final_df here

# Example: Print basic summary statistics
print(final_df.describe())

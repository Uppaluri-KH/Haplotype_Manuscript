import sys
import pandas as pd
variants = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Variants_DP15/{}_final_DP.vcf".format(sys.argv[1], sys.argv[1]), comment= '#', sep = '\t', header=None, low_memory=False)
variants.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

df = pd.read_excel("~/Downloads/13genes_coordinates_haplotypes_07122023.xlsx", header = None)
df.columns = ['Gene', 'chromosome', 'Extended_Start_pos', 'Extended_End_pos']

# Step 1: Create a dictionary from the df DataFrame
chromosome_dict = {}
for _, row in df.iterrows():
    chromosome = row['chromosome']
    start_pos = row['Extended_Start_pos']
    end_pos = row['Extended_End_pos']
    gene = row['Gene']  # Assuming 'Gene' is the name of your gene column
    if chromosome not in chromosome_dict:
        chromosome_dict[chromosome] = []
    chromosome_dict[chromosome].append((start_pos, end_pos, gene))

# Step 2: Define a function to check coverage
def check_coverage(row):
    pos = row['POS']
    chromosome = row['CHROM']
    if chromosome in chromosome_dict:
        ranges = chromosome_dict[chromosome]
        for start, end, gene in ranges:
            if start <= pos <= end:
                return 'Covered', start, end, gene
    return 'Not_Covered', None, None, None  # Return None for start, end, and gene if not covered

# Step 3: Apply the function to create the new columns in data
variants['Covered/Not_Covered'], variants['Start_Pos_Covered'], variants['End_Pos_Covered'], variants['Gene'] = zip(*variants.apply(check_coverage, axis=1))

# Step 4: Create new columns for Covered rows
covered_rows = variants['Covered/Not_Covered'] == 'Covered'
variants.loc[covered_rows, 'Covered_Chromosome'] = variants.loc[covered_rows, 'CHROM']
variants.loc[covered_rows, 'Covered_Start_Pos'] = variants.loc[covered_rows, 'Start_Pos_Covered']
variants.loc[covered_rows, 'Covered_End_Pos'] = variants.loc[covered_rows, 'End_Pos_Covered']
variants.loc[covered_rows, 'Gene'] = variants.loc[covered_rows, 'Gene']

# Drop temporary columns
variants.drop(['Start_Pos_Covered', 'End_Pos_Covered'], axis=1, inplace=True)

# Display the DataFrame
variants = variants[variants['Covered/Not_Covered'] == 'Covered']

df1 = pd.read_excel("~/Downloads/main_new_data_updated.xlsx")
df_x = df1[['CHROM', 'POS', 'REF', 'ALT', 'Haplotype', 'Type']]

merged_df = pd.merge(variants, df_x, on=['CHROM', 'POS', 'REF', 'ALT'], how='left', sort = False)
merged_df['Type'].fillna('Snp', inplace=True)
merged_df['Haplotype'].fillna('Snp', inplace=True)


df2 = df1.copy()
df2['Gene'] = df2['Haplotype'].str.split('*').str.get(0)
df2 = df2[['CHROM', 'POS', 'Gene', 'Haplotype']]
df2 = df2.rename(columns={'Haplotype': 'Haplotype_mapped'})
df2 = df2.groupby(['CHROM', 'POS', 'Gene']).agg({'Haplotype_mapped': lambda x: ','.join(x.unique())}).reset_index()

# Merge datasets on the 'Gene' column
merged_df_new = pd.merge(merged_df, df2, on='Gene', suffixes=('_merged_df', '_df2'))

# Calculate the absolute difference between 'POS' values
merged_df_new['POS_Difference'] = abs(merged_df_new['POS_merged_df'] - merged_df_new['POS_df2'])

# Find the row with the minimum absolute difference for each gene
min_diff_df = merged_df_new.loc[merged_df_new.groupby('Gene')['POS_Difference'].idxmin()]

# Merge the minimum difference information back to df1
result_df = pd.merge(merged_df, min_diff_df[['Gene', 'POS_Difference', 'CHROM_df2', 'POS_df2']], on='Gene', how='left')

result_df['Haplotype_mapped'] = min_diff_df.set_index(['Gene', 'CHROM_df2', 'POS_df2']).loc[result_df.set_index(['Gene', 'CHROM_df2', 'POS_df2']).index, 'Haplotype_mapped'].values

# Display the result or save it to a new dataset
result_df['POS_Difference'] = result_df['POS'] - result_df['POS_df2']
mask = result_df['Haplotype'] != 'Snp'
result_df.loc[mask, 'POS_Difference'] = None
result_df.loc[mask, 'CHROM_df2'] = result_df.loc[mask, 'CHROM']
result_df.loc[mask, 'POS_df2'] = result_df.loc[mask, 'POS']

result_df = result_df.drop(['QUAL', 'FILTER', 'Type'], axis=1)
result_df['Haplotype_mapped'] = result_df['Haplotype_mapped'].apply(lambda x: ','.join(set(x.split(','))))
result_df.to_excel("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Key_Variants/Cond3_All_Star_Alleles/new_distance_files_mul_pos/{}_snp_pos_distance.xlsx".format(sys.argv[1], sys.argv[1]), index = False)
print("Condition3 is finished for the sample")
# print(result_df)

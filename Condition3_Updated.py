import sys
import pandas as pd
variants = pd.read_csv("/{}_final_DP.vcf".format(sys.argv[1], sys.argv[1]), comment= '#', sep = '\t', header=None, low_memory=False)
variants.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
df = pd.read_excel("/13genes_coordinates_haplotypes_07122023.xlsx", header = 0)
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
df1 = pd.read_excel("~/Downloads/Madhu_Updated/main_new_data.xlsx")
df_x = df1[['CHROM', 'POS', 'REF', 'ALT', 'Haplotype', 'Type']]
df_x = df_x.groupby(['CHROM', 'POS', 'REF', 'ALT']).agg({'Haplotype': lambda x: ','.join(x.unique()), 
                                                        'Type' : 'first'}).reset_index()

merged_df = pd.merge(variants, df_x, on=['CHROM', 'POS', 'REF', 'ALT'], how='left', sort = False)
merged_df['Type'].fillna('Snp', inplace=True)
merged_df['Haplotype'].fillna('Snp', inplace=True)
df2 = df1.copy()
df2['Gene'] = df2['Haplotype'].str.split('*').str.get(0)
df2 = df2[['CHROM', 'POS', 'Gene']]
df2 = df2.drop_duplicates(subset=['CHROM', 'POS', 'Gene'])
def find_nearest_pos(row):
    df2_subset = df2[(df2['Gene'] == row['Gene']) & (df2['CHROM'] == row['CHROM'])]
    if not df2_subset.empty:
        nearest_pos_index = (df2_subset['POS'] - row['POS']).abs().idxmin()
        nearest_pos = df2_subset.loc[nearest_pos_index, 'POS']
        return nearest_pos
    else:
        return None

# Apply the function to create 'POS_df2' column in data
merged_df['POS_df2'] = merged_df.apply(find_nearest_pos, axis=1)
merged_df['POS_diff'] = merged_df['POS_df2'] - merged_df['POS']
merged_df['CHROM_df2'] = merged_df['Covered_Chromosome']
df3 = pd.read_excel("~/Downloads/Madhu_Updated/main_new_data.xlsx")
df3 = df3.groupby(['CHROM', 'POS']).agg({'Haplotype': lambda x: ','.join(x.unique())}).reset_index()
df3.rename(columns={'POS':'POS_df2', 'Haplotype':'Haplotype_updated'}, inplace=True)
final = pd.merge(merged_df, df3, on = ['CHROM', 'POS_df2'], how = 'left', sort = False)
final['Haplotype_updated'] = final['Haplotype_updated'].apply(lambda x: ','.join(sorted(set(x.split(',')))))
final.to_excel("/{}_snp_pos_distance.xlsx".format(sys.argv[1], sys.argv[1]), index = False)
print("Condition3 is finished for the sample")

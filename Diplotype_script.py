import pandas as pd
import sys

df1 = pd.read_excel("~/Downloads/main_new_data.xlsx")
df1

import pandas as pd

vcf = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Key_Variants/Cond2_Key_Star_Alleles/Out_Vcfs/{}_final_DP_Cond2.vcf".format(sys.argv[1], sys.argv[1]), comment= '#', sep = '\t', header=0, usecols=[0,1,2,3,4,5,6,7,8,9],low_memory=False)
# vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)
vcf['HET'] = vcf['INFO'].str.extract(r'HET=(\d)')
vcf['HOM'] = vcf['INFO'].str.extract(r'HOM=(\d)')

# Create a new column 'Zygosity' based on conditions
vcf['Zygosity'] = ''

vcf.loc[vcf['HOM'] == '1', 'Zygosity'] = 'Homozygous'
vcf.loc[vcf['HET'] == '1', 'Zygosity'] = 'Heterozygous'
vcf = vcf[['CHROM', 'POS', 'rsID_x', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL','RDF', 'RDR', 'ADF', 'ADR', 'Zygosity']]
vcf

df = pd.read_csv("~/Downloads/KAPA_HyperExome_hg38_primary_targets_extended.bed", sep = '\t', header = None, error_bad_lines=False)
df.columns = ['chromosome', 'Start_pos', 'End_pos']
df['Extended_Start_pos'] = df['Start_pos'] - 20
df['Extended_End_pos'] = df['End_pos'] + 20
df = df[['chromosome', 'Extended_Start_pos', 'Extended_End_pos']]
# vcf['POS']=vcf['POS'].astype(int(float))
# df['Extended_Start_pos']=df['Extended_Start_pos'].astype(int(float))
# df['Extended_End_pos']=df['Extended_End_pos'].astype(int(float))

# Step 1: Create a dictionary from the df DataFrame
chromosome_dict = {}
for _, row in df.iterrows():
    chromosome = row['chromosome']
    start_pos = row['Extended_Start_pos']
    end_pos = row['Extended_End_pos']
    if chromosome not in chromosome_dict:
        chromosome_dict[chromosome] = []
    chromosome_dict[chromosome].append((start_pos, end_pos))

# Step 2: Define a function to check coverage
def check_coverage(row):
    pos = row['POS']
    chromosome = row['CHROM']
    if chromosome in chromosome_dict:
        ranges = chromosome_dict[chromosome]
        for start, end in ranges:
            if start <= pos <= end:
                return 'Covered'
    return 'Not_Covered'

# Step 3: Apply the function to create the new column in data
vcf['Covered/Not_Covered'] = vcf.apply(check_coverage, axis=1)
vcf

merged = pd.merge(vcf, df1, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'inner', sort=False)
merged

df_1 = vcf.copy()
df_2 = merged.copy()

df_2['Gene'] = df_2['Haplotype'].str.split('*').str[0]
df_2

df_2 = df_2.drop_duplicates(subset=['Gene', 'Haplotype'], keep='first').reset_index(drop=True)
df_2['Haplotype_new'] = "*" + df_2['Haplotype'].str.split('*').str[1]
df_2

import pandas as pd
from itertools import combinations_with_replacement
# Function to generate combinations
def get_combinations(group):
    return list(combinations_with_replacement(group, 2))

# Apply the function to each group of rsID
result = df_2.groupby('Gene')['Haplotype'].apply(get_combinations).explode().reset_index(name='Diplotype')

# Extract the relevant portions of the haplotypes and keep the asterisk (*)
result['Diplotype'] = result['Diplotype'].apply(lambda x: '/'.join([hap.replace(hap.split('*')[0], '') if len(hap.split('*')) > 1 else hap for hap in x]))

# Filter out reversed combinations
result = result[result['Diplotype'].apply(lambda x: x[0] != x[2])]

# Drop duplicates
result = result.drop_duplicates(subset=['Gene', 'Diplotype'], keep='first').reset_index(drop=True)

# Print the result
result

df2 = df_2.copy()
df2 = df2.rename(columns={'Haplotype_new': 'Diplotype', 'Covered/Not_Covered':'Coverage'})
df2 = df2[['Gene', 'Diplotype', 'Coverage', 'Zygosity']]
df2

# Create an empty 'Coverage' column in dataset1
result['Coverage'] = ''
result['Zygosity'] = ''

# Iterate over rows in dataset1
for index, row in result.iterrows():
    gene = row['Gene']
    diplotype_parts = row['Diplotype'].split('/')
    
    # Check if each part of the diplotype is covered in dataset2
    coverages = []
    zygosities = []
    for part in diplotype_parts:
        coverage = df2.loc[(df2['Gene'] == gene) & (df2['Diplotype'] == part), 'Coverage'].values
        zygosity = df2.loc[(df2['Gene'] == gene) & (df2['Diplotype'] == part), 'Zygosity'].values
        
        coverages.append(coverage[0] if coverage else 'Not Covered')
        zygosities.append(zygosity[0] if zygosity else 'Not Specified')

    # Set the 'Coverage' and 'Zygosity' column values in dataset1
    result.at[index, 'Coverage'] = '/'.join(coverages)
    result.at[index, 'Zygosity'] = '/'.join(zygosities)

result[['Gene', 'Diplotype', 'Coverage', 'Zygosity']]
result

data = pd.read_excel("~/Downloads/Diplotype_phenotype_Updated.xlsx")
data

result.Gene.value_counts()

df_merged = pd.merge(result, data, on = ['Gene', 'Diplotype'], how = 'left', sort = False)
df_merged

df_1 = pd.read_excel("~/Downloads/Star_allele_function_Updated.xlsx")
df_1

# Function to map Diplotype to Function
def map_diplotype_to_function(row):
    gene = row['Gene']
    diplotype_parts = row['Diplotype'].split('/')
    
    function_parts = [
        df_1.loc[(df_1['Gene'] == gene) & (df_1['Allele'] == part), 'Function'].values[0] 
        if not df_1.loc[(df_1['Gene'] == gene) & (df_1['Allele'] == part), 'Function'].empty 
        else '-' for part in diplotype_parts
    ]
    
    return '/'.join(function_parts)

# Apply the mapping function to create a new column
df_merged['Function'] = df_merged.apply(map_diplotype_to_function, axis=1)

# Display the result
result_df = df_merged[['Gene', 'Diplotype', 'Coverage','Function', 'Zygosity', 'Coded Diplotype/Phenotype Summary']]
result_df.to_excel("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Key_Variants/Cond2_Key_Star_Alleles/Out_Vcfs/Diplotypes/{}_Cond2_Diplotypes.xlsx".format(sys.argv[1], sys.argv[1]), index=False)
result_df
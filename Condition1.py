import pandas as pd
import sys

### Covered/Not_Covered
Bed = pd.read_csv("~/Downloads/KAPA_HyperExome_hg38_primary_targets_extended.bed", sep = '\t', header = None, error_bad_lines=False)
Bed.columns = ['chromosome', 'Start_pos', 'End_pos']
Bed['Extended_Start_pos'] = Bed['Start_pos'] - 20
Bed['Extended_End_pos'] = Bed['End_pos'] + 20
Bed = Bed[['chromosome', 'Extended_Start_pos', 'Extended_End_pos']]

# Step 1: Create a dictionary from the Bed DataFrame
chromosome_dict = {}
for _, row in Bed.iterrows():
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

# ### Condition1: 13 Genes Main file
# print("### Condition1: 13 Genes Main file")
# Genes_13 = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/check/13_Genes_Key_Variants_Coordinates.csv", sep="\t")
# vcf = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Variants_DP15/{}_final_DP.vcf".format(sys.argv[1], sys.argv[1]), sep="\t", comment="#")
# vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
# Genes_vcf = pd.merge(vcf, Genes_13, on=['CHROM', 'POS', 'REF', 'ALT'], how="inner", sort=False)

# Genes_vcf['Covered/Not_Covered'] = Genes_vcf.apply(check_coverage, axis=1)
# Genes_vcf.to_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Key_Variants/Cond1_All_Star_Alleles/Intermediate/{}_final_DP_Cond1.csv".format(sys.argv[1], sys.argv[1]), sep="\t", index=False)
# print("Condition1 is done and placed output")

### Condition2: Single shared file
print("### Condition2: Single shared file")
Shared_single = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/check/13_genes_single_shared_Updated.csv", sep="\t")
vcf = pd.read_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Variants_DP15/{}_final_DP.vcf".format(sys.argv[1], sys.argv[1]), sep="\t", comment="#")
vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
Shared_vcf = pd.merge(vcf, Shared_single, on=['CHROM', 'POS', 'REF', 'ALT'], how="inner", sort=False)

Shared_vcf['Covered/Not_Covered'] = Shared_vcf.apply(check_coverage, axis=1)

Shared_vcf.to_csv("/run/media/administrator/Expansion1/azure_backup/final_vcf/s3_uploaded_snps/1721_samples/PRX/Pharmacogenomics/PGx_Genes_Key_Variants/Cond2_Key_Star_Alleles/Intermediate/{}_final_DP_Cond2.csv".format(sys.argv[1], sys.argv[1]), sep="\t", index=False)
print("Condition2 is done and placed output")


import os
import pandas as pd

# Set the path to the folder containing your TSV files
#folder_path = r'C:/Users/GenepoweRx_Madhu/Downloads/haplo_20_samples/'
folder_path = '.'

# Initialize an empty DataFrame to store the concatenated data
concatenated_data = pd.DataFrame()

# Iterate over each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".tsv"):
        file_path = os.path.join(folder_path, filename)

        # Read the TSV file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')

        # Extract the "Sample" value from the file name
        sample_name = filename.split('_')[0]

        # Add a new column "Sample" with the extracted value
        df['Sample'] = sample_name

        # Extract 'DP' from the 'SAMPLE' column
        df['DP'] = df['SAMPLE'].str.split(':').str[3].astype(float)

        # Concatenate the data to the main DataFrame
        concatenated_data = pd.concat([concatenated_data, df], ignore_index=True)

# Group by specified columns and aggregate unique 'Sample' values and count
grouped_data = concatenated_data.groupby(['CHROM', 'POS', 'REF', 'ALT', 'Haplotype', 'rsID_y', 'Zygosity'], as_index=False).agg({
    'Sample': lambda x: ','.join(x.unique()),  # Comma-separated unique samples
    'DP': 'mean'
})

# Rename the 'Sample' column to 'Sample Count'
grouped_data.rename(columns={'Sample': 'Sample_List'}, inplace=True)
grouped_data['Sample_count'] = grouped_data['Sample_List'].apply(lambda x: len(x.split(',')))

# Save the grouped DataFrame to a new TSV file
#grouped_data.to_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/haplo_20_samples/grouped_data.tsv', sep='\t', index=False)
grouped_data.to_csv("Final_Counts_Condition1a.tsv", sep='\t', index=False)
grouped_data
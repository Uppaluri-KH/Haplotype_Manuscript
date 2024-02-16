import os
import pandas as pd

# Set the path to the folder containing your Excel files
folder_path = "." 

# Initialize an empty list to store individual DataFrames
dfs = []

# Iterate over each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".xlsx"):
        file_path = os.path.join(folder_path, filename)
        
        # Read the Excel file into a DataFrame
        df = pd.read_excel(file_path)
        
        # Extract the "Sample" value from the base name of the file
        sample_name = os.path.splitext(filename)[0]
        
        # Add a new column "Sample" with the extracted value
        df['Sample'] = sample_name
        #df['Sample'] = df['Sample'].str.replace('_snp', '')

        # Append the DataFrame to the list
        dfs.append(df)

# Concatenate all DataFrames in the list
concatenated_data = pd.concat(dfs, ignore_index=True)

# Save the final DataFrame to a new Excel file
output_file_path = "./All_1232_updated_haplo.xlsx"
concatenated_data.to_excel(output_file_path, index=False)

# Print the concatenated DataFrame
concatenated_data

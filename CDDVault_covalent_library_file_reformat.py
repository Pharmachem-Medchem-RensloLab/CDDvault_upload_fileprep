import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from datetime import datetime

# Function to convert SMILES to image
def smiles_to_image(smiles, title=""):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        return img
    else:
        return None

# Main function to process data
def process_data(project, chemist, input_file_1, input_file_2):
    try:
        # Check the extension of input_file_1
        file_ext_1 = os.path.splitext(input_file_1)[1].lower()
        if file_ext_1 != '.csv':
            raise ValueError(f"Error: input_file_1 must be a CSV file. Found: {file_ext_1}")

        # Load input files
        df1 = pd.read_csv(input_file_1)

        # Check the extension of input_file_2
        file_ext = os.path.splitext(input_file_2)[1].lower()
        if file_ext == '.csv':
            df2 = pd.read_csv(input_file_2)
        elif file_ext == '.xlsx':
            df2 = pd.read_excel(input_file_2)
        else:
            raise ValueError(f"Unsupported file format for input_file_2, needs to be either .csv or .xlsx: {file_ext}")

        # Display the first few rows of each DataFrame
        print("Input File 1:")
        display(df1.head())

        print("Input File 2:")
        display(df2.head())

        # Ensure the correct columns are present
        required_columns_df1 = [
            "Container Id", "Orientation Barcode", "Row", "Column", "Barcode",
            "Scan Time"
        ]
        required_columns_df2 = [
            "VIAL_QR_CODE", "SYNONYMS", "SMILES", "INITIAL_VOLUME_UL", "CONC_mM",
            "SALT", "RLA_Number"  # Added "RLA Number" to the required columns
        ]

        missing_columns_df1 = [col for col in required_columns_df1 if col not in df1.columns]
        missing_columns_df2 = [col for col in required_columns_df2 if col not in df2.columns]

        if missing_columns_df1:
            raise ValueError(
                f"Error: Input file 1 is missing the following required columns: "
                f"{missing_columns_df1}"
            )
        if missing_columns_df2:
            raise ValueError(
                f"Error: Input file 2 is missing the following required columns: "
                f"{missing_columns_df2}"
            )

        # Attempt to merge and check specifically for unmatched entries in 'VIAL_QR_CODES'
        merged_df = pd.merge(df1, df2, left_on="Barcode", right_on="VIAL_QR_CODE", how="outer", indicator=True)
        # Check for unmatched 'VIAL_QR_CODES'
        unmatched_df2 = merged_df[merged_df['_merge'] == 'right_only']
        if not unmatched_df2.empty:
            # Convert 'VIAL_QR_CODE' from float to integer for display
            unmatched_df2['VIAL_QR_CODE'] = unmatched_df2['VIAL_QR_CODE'].astype(int)
            # Generate error message listing all unmatched 'VIAL_QR_CODES'
            error_message = "Mismatch found: The following 'VIAL_QR_CODES' in File 2 did not match any 'Barcode' in File 1:\n"
            error_message += ', '.join(unmatched_df2['VIAL_QR_CODE'].astype(str).tolist())
            error_message += "\nPlease double check these 'VIAL_QR_CODE' values."
            raise ValueError(error_message)
        # Proceed with rows that matched in both
        merged_df = merged_df[merged_df['_merge'] == 'both'].drop(columns=['_merge'])

        # Create the 'PLATE_WELL' column
        merged_df["PLATE_WELL"] = merged_df["Row"].astype(str) + merged_df["Column"].astype(str)

        # Select and rename the desired columns for the output
        output_df = merged_df[[
            "Container Id", "PLATE_WELL", "INITIAL_VOLUME_UL", "CONC_mM", "VIAL_QR_CODE",
            "SMILES", "SYNONYMS", "RLA_Number"  # Added "RLA Number" to the output columns
        ]].copy()

        output_df.rename(columns={
            "Container Id": "PLATE_BARCODE",
            "SMILES": "CDD_SMILES"
        }, inplace=True)

        # Add the user-defined variables
        output_df["Project"] = project
        output_df["Chemist"] = chemist

        # Reorder columns
        output_df = output_df[[
            "Project", "Chemist", "PLATE_BARCODE", "PLATE_WELL", "INITIAL_VOLUME_UL",
            "CONC_mM", "VIAL_QR_CODE", "CDD_SMILES", "SYNONYMS", "RLA_Number"  # Added "RLA Number" to the reordered columns
        ]]

        # Display the first few rows of the output DataFrame
        display(output_df.head())

        # Handling the DataFrame for SMILES visualization
        # Convert 'VIAL_QR_CODE' from float to integer if necessary
        merged_df['VIAL_QR_CODE'] = merged_df['VIAL_QR_CODE'].astype(int)

        # Create a new DataFrame for displaying images
        smiles_data = merged_df[['VIAL_QR_CODE', 'SYNONYMS', 'SMILES']].copy()
        smiles_data['Image'] = smiles_data['SMILES'].apply(lambda x: smiles_to_image(x))

        # Display the table with images (as much as possible in a text-based environment)
        for index, row in smiles_data.iterrows():
            print(f"\nVIAL_QR_CODE: {row['VIAL_QR_CODE']}, SYNONYMS: {row['SYNONYMS']}, SMILES: {row['SMILES']}")
            plt.imshow(row['Image'])
            plt.xticks([])  # Disable x-axis ticks
            plt.yticks([])  # Disable y-axis ticks
            plt.title(f"VIAL_QR_CODE: {row['VIAL_QR_CODE']}")
            plt.show()

        # Save Output File Section
        # Get the current date in YYYYMMDD format
        current_date = datetime.now().strftime("%Y%m%d")

        # Define the output file path
        output_file_path = f"CDDupload_input_file_{current_date}.csv"
        # Saves file in the current working directory

        # Save the output DataFrame to a CSV file
        output_df.to_csv(output_file_path, index=False)
        print(f"Output file saved to: {output_file_path}")

    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")

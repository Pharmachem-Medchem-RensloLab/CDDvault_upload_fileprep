# Import libraries
import pandas as pd
from datetime import datetime

# Define functions

## General Check for File Availability
def check_file(path):
    try:
        return pd.read_excel(path, header=None)
    except FileNotFoundError:
        print(f"File at {path} is missing.")
        return None

## Reformat Protein Percent Bound Data
def percent_protein_bound_data_reformat(path, default_batch_value):
    df = check_file(path)
    if df is not None:
        df.columns = df.iloc[0]
        df = df.iloc[1:]
        df['Batch'] = default_batch_value
        required_columns = ['Compound', 'Concentration (uM)', 'Batch', '% Bound Average', '%Recovery', 'Species']
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            print(f"Warning: Missing expected columns {missing_cols}")
            return
        df = df[df['Compound'].notna() & df['% Bound Average'].notna()]
        float_columns = ['Concentration (uM)', '% Bound Average', '%Recovery']
        for col in float_columns:
            df[col] = df[col].apply(lambda x: round(float(x), 2) if isinstance(x, (int, float)) else x)
        df = df[required_columns]
        print("Preview of Percent Protein Bound Data:")
        print(df.head())
        output_filename = f'UCSF_Percent_Protein_Bound_Results_CDDformat_{datetime.now().strftime("%Y%m%d")}.csv'
        df.to_csv(output_filename, index=False)
        print(f"Output file created: {output_filename}")
    else:
        print("Processing skipped due to missing file.")

## Reformat Kinetic Solubility Data
def kinetic_solubility_data_reformat(path, default_batch_value):
    df = check_file(path)
    if df is not None:
        buffer_info = df.iloc[1, 2]
        df.columns = ['Compound', 'Final [DMSO]', f'KSOL (uM) {buffer_info} Value', f'KSOL (uM) {buffer_info} Mean']
        df = df.iloc[2:]
        df['Compound'] = df['Compound'].ffill()
        df['Final [DMSO]'] = df['Final [DMSO]'].ffill()
        df = df[df['Compound'].notna() & df[f'KSOL (uM) {buffer_info} Mean'].notna()]
        df[f'KSOL (uM) {buffer_info} Mean'] = df[f'KSOL (uM) {buffer_info} Mean'].apply(
            lambda x: round(x, 1) if isinstance(x, float) else x
        )
        df['Batch'] = default_batch_value
        selected_columns = ['Compound', 'Final [DMSO]', f'KSOL (uM) {buffer_info} Mean', 'Batch']
        df = df[selected_columns]
        print("Preview of the Kinetic Solubility Data:")
        print(df.head())
        output_filename = f'UCSF_KSOL_Results_CDDformat_{datetime.now().strftime("%Y%m%d")}.csv'
        df.to_csv(output_filename, index=False)
        print(f"Output file created: {output_filename}")
    else:
        print("Processing skipped due to missing file.")

## Reformat Liver Microsome Stability Data
def liver_microsome_stability_data_reformat(path, default_batch_value):
    df = check_file(path)
    if df is not None:
        special_header = pd.read_excel(path, header=None).iloc[0, 11]
        headers = df.iloc[0].ffill()
        subheaders = df.iloc[1].replace({pd.NA: ''})
        df.columns = [f'{h}' if str(sh).strip() == '' else f'{h}_{sh}' for h, sh in zip(headers, subheaders)]
        df = df.iloc[2:]
        print("Column names after processing:", df.columns)
        df['Batch'] = default_batch_value
        df[special_header] = pd.NA
        required_columns = ['Compound', 'Batch', 't1/2 (min)_Mean', 't1/2 (min)_SE',
                            'CLint (µL/min/mg protein)_Mean', 'CLint (µL/min/mg protein)_SE', special_header]
        def format_and_round(value, decimals):
            if isinstance(value, (int, float)):
                return round(float(value), decimals)
            return value
        df['t1/2 (min)_Mean'] = df['t1/2 (min)_Mean'].apply(lambda x: format_and_round(x, 1))
        df['t1/2 (min)_SE'] = df['t1/2 (min)_SE'].apply(lambda x: format_and_round(x, 1))
        df['CLint (µL/min/mg protein)_Mean'] = df['CLint (µL/min/mg protein)_Mean'].apply(lambda x: format_and_round(x, 1))
        df['CLint (µL/min/mg protein)_SE'] = df['CLint (µL/min/mg protein)_SE'].apply(lambda x: format_and_round(x, 3))
        missing_cols = [col for col in required_columns if col not in df.columns]
        if missing_cols:
            print(f"Warning: Missing expected columns {missing_cols}")
        df = df[required_columns]
        print("Preview of Liver Microsome Stability Data:")
        print(df.head())
        output_filename = f'UCSF_Liver_Microsome_Stability_Results_CDDformat_{datetime.now().strftime("%Y%m%d")}.csv'
        df.to_csv(output_filename, index=False)
        print(f"Output file created: {output_filename}")
    else:
        print("Processing skipped due to missing file.")

## Reformat Caco-2 Permeability Data
def caco2_permeability_data_reformat(path, default_batch_value):
    df = check_file(path)
    if df is not None:
        headers = df.iloc[0].ffill()
        subheaders = df.iloc[1].replace({pd.NA: ''})
        df.columns = [f'{h}' if str(sh).strip() == '' else f'{h} {sh}' for h, sh in zip(headers, subheaders)]
        df = df.iloc[2:]
        print("Column names after processing:", df.columns)
        df['Batch'] = default_batch_value
        required_columns = ['Compound', 'Batch', 'Papp, A-B (x10-6 cm/s) Mean',
                            'Papp, B-A (x10-6 cm/s) Mean', 'Ratio\nB-A/A-B', 'Recovery (%)']
        if not all(col in df.columns for col in required_columns):
            missing_cols = [col for col in required_columns if col not in df.columns]
            print(f"Warning: Missing expected columns {missing_cols}")
            return
        df = df[df['Compound'].notna() & df['Papp, A-B (x10-6 cm/s) Mean'].notna()]
        float_columns = ['Papp, A-B (x10-6 cm/s) Mean', 'Papp, B-A (x10-6 cm/s) Mean', 'Ratio\nB-A/A-B']
        for col in float_columns:
            df[col] = df[col].apply(lambda x: round(x, 2) if isinstance(x, float) else x)
        df['Recovery (%)'] = df['Recovery (%)'].apply(lambda x: round(x) if isinstance(x, float) else x)
        if not all(col in df.columns for col in required_columns):
            missing_cols = [col for col in required_columns if col not in df.columns]
            print(f"Warning: Missing expected columns {missing_cols}")
        else:
            df = df[required_columns]
            print("Preview of Caco-2 Permeability Data:")
            print(df.head())
            output_filename = f'UCSF_Caco-2_Results_CDDformat_{datetime.now().strftime("%Y%m%d")}.csv'
            df.to_csv(output_filename, index=False)
            print(f"Output file created: {output_filename}")
    else:
        print("Processing skipped due to missing file.")

## Reformat MDCK Permeability Data
def mdck_permeability_data_reformat(path, default_batch_value):
    df = check_file(path)
    if df is not None:
        cell_line = pd.read_excel(path, header=None).iloc[0, 0]
        headers = df.iloc[1].ffill()
        subheaders = df.iloc[2].replace({pd.NA: ''})
        df.columns = [f'{h}' if str(sh).strip() == '' else f'{h} {sh}' for h, sh in zip(headers, subheaders)]
        df = df.iloc[3:]
        print("Column names after processing:", df.columns)
        df['Batch'] = default_batch_value
        df['Cell Line'] = pd.NA
        df.loc[df.index[0], 'Cell Line'] = cell_line
        required_columns = ['Compound', 'Batch', 'Papp, A-B (x10-6 cm/s) Mean', 'Papp, B-A (x10-6 cm/s) Mean',
                            'Ratio\nB-A/A-B', 'Recover Rate (%)', 'Cell Line']
        if not all(col in df.columns for col in required_columns):
            missing_cols = [col for col in required_columns if col not in df.columns]
            print(f"Warning: Missing expected columns {missing_cols}")
            return
        df = df[df['Compound'].notna() & df['Papp, A-B (x10-6 cm/s) Mean'].notna() & df['Papp, B-A (x10-6 cm/s) Mean'].notna()]
        float_columns = ['Papp, A-B (x10-6 cm/s) Mean', 'Papp, B-A (x10-6 cm/s) Mean', 'Ratio\nB-A/A-B']
        for col in float_columns:
            df[col] = df[col].apply(lambda x: round(x, 2) if isinstance(x, float) else x)
        df['Recover Rate (%)'] = df['Recover Rate (%)'].apply(lambda x: round(x) if isinstance(x, float) else x)
        df = df[required_columns]
        print("Preview of MDCK Permeability Data:")
        print(df.head())
        output_filename = f'UCSF_MDCK_Permeability_Results_CDDformat_{datetime.now().strftime("%Y%m%d")}.csv'
        df.to_csv(output_filename, index=False)
        print(f"Output file created: {output_filename}")
    else:
        print("Processing skipped due to missing file.")

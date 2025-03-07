# CDDvault Upload File Reformatting
Contains two scripts:

(1) Reformat CRO small molecule ADME data for upload to CDDVault
        The script here is called by a Google Colab file to process ADME small molecule information from Quintara.
        The script can handle data for kinetic solubility, liver microsome stability, Caco2 permeability, MDCK permeability, and 
        plasma protein binding.

(2) Reformat user generated files to add covalent small molecules to bespoke library in CDDVault
        The script requires two input files. The first file is automatically generated by the barcode reader in the SMDC (CSV file).
        The second file is a user generated file containing information about the small molecules required for upload: vial QR code, 
        synonym(s), compound SMILE string, initial compound volume, compound concentration, and salt information. 

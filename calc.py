import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def process_smiles_file(input_csv, output_csv):
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv)

    # Initialize lists to store the calculated properties
    mol_wt = []
    log_p = []
    h_bond_donors = []
    h_bond_acceptors = []
    molar_refractivity = []
    lipinski = []

    # Process each SMILES string
    for smiles in df['SMILES']:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # If the SMILES string is invalid, append None for each property and "No" for Lipinski
            mol_wt.append(None)
            log_p.append(None)
            h_bond_donors.append(None)
            h_bond_acceptors.append(None)
            molar_refractivity.append(None)
            lipinski.append("No")
            continue

        # Calculate properties and round to 2 decimal places
        mw = round(Descriptors.MolWt(mol), 2)
        lp = round(Descriptors.MolLogP(mol), 2)
        donors = Descriptors.NumHDonors(mol)
        acceptors = Descriptors.NumHAcceptors(mol)
        mr = round(Descriptors.MolMR(mol), 2)

        # Append calculated properties
        mol_wt.append(mw)
        log_p.append(lp)
        h_bond_donors.append(donors)
        h_bond_acceptors.append(acceptors)
        molar_refractivity.append(mr)

        # Check Lipinski rule of five
        if (mw <= 500 and
            donors <= 5 and
            acceptors <= 10 and
            lp <= 5 and
            40 <= mr <= 130):
            lipinski.append("Yes")
        else:
            lipinski.append("No")

    # Add new columns to the DataFrame
    df['mol_wt'] = mol_wt
    df['log_p'] = log_p
    df['h_bond_donors'] = h_bond_donors
    df['h_bond_acceptors'] = h_bond_acceptors
    df['molar_refractivity'] = molar_refractivity
    df['Lipinski'] = lipinski

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Processed data saved to {output_csv}")

# Example usage
process_smiles_file('file_metergoline_mod.csv', 'output_drugs_with_properties.csv')

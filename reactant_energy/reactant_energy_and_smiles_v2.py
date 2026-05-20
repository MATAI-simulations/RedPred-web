import pickle

import psycopg
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator

import pandas as pd

from electrolyte_db.db import * 


# connect to mp_materials database 
connection = get_connection()
cursor = connection.cursor()

records = query_molecule(table="pubchem_molecules", result_cols=["iupac_name", "smiles"])


# separating the data
iupac_names = []
smiles_list = []

for r in records:
    iupac_name = r['iupac_name']
    smiles = r['smiles']

    # skipping null/empty SMILES
    if smiles:
        iupac_names.append(iupac_name)
        smiles_list.append(smiles)


# fingerprint generation
def get_ecfc(smiles_list, radius=2, nBits=2048):
    """
        Computes the Extended Connectivity Fingerprint (ECFP) for a given molecule, starting with iteration 0
        for just the atom identifiers, then expanding in further iterations. Additionally, collects all unique
        identifiers generated throughout the iterations.

        Parameters:
        - smiles_list: A list of smiles molecules 
        - radius: The radius of neighborhoods to consider, corresponding to the iteration count.
        - nBits: Int > 0; Length of the fingerprint bitarray.

        Returns:
        - A list of bitarray representing the ECFP of the molecule.

    """
    mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)

    ecfp_fingerprints = []
    valid_indices = []

    for idx, smiles in enumerate(smiles_list):

        mol = Chem.MolFromSmiles(smiles)

        

        if mol is not None:
            mol = Chem.AddHs(mol)

            # fingerprint changes the redox potentials 
            fp = mfgen.GetFingerprintAsNumPy(mol)

            ecfp_fingerprints.append(fp)
            valid_indices.append(idx)

        else:
            print(f"Invalid SMILES skipped: {smiles}")

    df_fp = pd.DataFrame(ecfp_fingerprints)

    return df_fp, valid_indices


# generating fingerprints
ecfc_encoder, valid_indices = get_ecfc(smiles_list)


# keeping only valid molecules
valid_iupac = [iupac_names[i] for i in valid_indices]
valid_smiles = [smiles_list[i] for i in valid_indices]


# loading ML model
rf_final_model = pickle.load(
    open('./final_models/rf_final_model.txt', "rb")
)

# predicting reaction energies
pred_rf = enumerate(rf_final_model.predict(ecfc_encoder))
pred_rf = sorted(pred_rf, key=lambda x: x[1])

# final data frame
results_df = pd.DataFrame({
    "iupac_name": [valid_iupac[i] for i, _ in pred_rf],
    "smiles": [valid_smiles[i] for i, _ in pred_rf],
    "reaction_energy": [x[1] for x in pred_rf]
})

# saving as a csv
results_df.to_csv("sorted_reaction_energy_predictions.csv", index=False)

print(results_df.head())
print("\nCSV saved as sorted_reaction_energy_predictions.csv")
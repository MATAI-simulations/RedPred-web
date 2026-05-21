import pickle

import psycopg
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator

import pandas as pd

from electrolyte_db.connection import * 
from electrolyte_db.query import * 


# connect to mp_materials database 
connection = get_connection()
cursor = connection.cursor()

records = query_molecule(table="pubchem_molecules", result_cols=["iupac_name", "n_atoms", "molecular_formula", "smiles"])


# retrieving only the records with a valid SMILES. Returns a list of dict where each dict corresponds to a row in the pubchem_molecules database table. The keys are the column names and values are entries associated with that particular column 
valid_records = []

for r in records:
    smiles = r['smiles']

    # skipping null/empty SMILES
    if smiles:
        valid_records.append(r)


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


# generating fingerprints for all smiles molecules 
smiles_list = [d['smiles'] for d in valid_records]
ecfc_encoder, valid_indices = get_ecfc(smiles_list)


# loading ML model
rf_final_model = pickle.load(
    open('./final_models/rf_final_model.txt', "rb")
)


# predicting reaction energies and storing them in valid_record (list of rows in database)
pred_rf = enumerate(rf_final_model.predict(ecfc_encoder))
pred_rf = sorted(pred_rf, key=lambda x: x[1])

for i in pred_rf:
    valid_records[i[0]]["reduction_potential"] = i[1] 


# final data frame
results_df = pd.DataFrame({
    "iupac_name": [valid_records[i]['iupac_name'] for i, _ in pred_rf],
    "n_atoms": [valid_records[i]['n_atoms'] for i, _ in pred_rf],
    "molecular_formula": [valid_records[i]['molecular_formula'] for i, _ in pred_rf],
    "smiles": [smiles_list[i] for i, _ in pred_rf],
    "reaction_energy": [x[1] for x in pred_rf]
})

# saving as a csv
results_df.to_csv("sorted_reaction_energy_predictions.csv", index=False)

print(results_df.head())
print("\nCSV saved as sorted_reaction_energy_predictions.csv")


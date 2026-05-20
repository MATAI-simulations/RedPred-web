import pickle

import psycopg
from rdkit import Chem
from rdkit.Chem import AllChem, rdFingerprintGenerator

import pandas as pd
import os 

# connection = psycopg2.connect(database="mp_materials", user="mp_materials_user", password="sEYXA8FeVc3vNnnOfXUpCFRX8E8dayn8", host="dpg-d74s49kr85hc73fu7620-a.oregon-postgres.render.com", port=5432)


db_url = os.environ["DATABASE_URL"]

connection = psycopg.connect(db_url, sslmode="require")

cursor = connection.cursor()

# getting the iupac_name and smiles
cursor.execute("""
    SELECT iupac_name, smiles
    FROM pubchem_molecules
""")

records = cursor.fetchall()

# separating the data
iupac_names = []
smiles_list = []

for r in records:
    iupac_name = r[0]
    smiles = r[1]

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

            # fp = list(
            #     AllChem.GetHashedMorganFingerprint(
            #         mol,
            #         radius,
            #         nBits
            #     )
            # )
    
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
pred_rf = rf_final_model.predict(ecfc_encoder)
# print(f'Predicted reaction energies are {pred_rf}')
# print(type(pred_rf))

# pred_rf = sorted(enumerate(pred_rf), key=lambda x: x[1])
# print(pred_rf)

# final data frame
results_df = pd.DataFrame({
    "iupac_name": valid_iupac,
    "smiles": valid_smiles,
    "reaction_energy": pred_rf
})

# saving as a csv
results_df.to_csv("reaction_energy_predictions.csv", index=False)

print(results_df.head())
print("\nCSV saved as reaction_energy_predictions.csv")
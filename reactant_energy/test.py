from electrolyte_db.connection import * 
from electrolyte_db.query import * 
from rdkit import Chem 

# connect to mp_materials database 
connection = get_connection()
cursor = connection.cursor()

records = query_molecule(table="pubchem_molecules", result_cols=["iupac_name", "n_atoms", "molecular_formula", "smiles"])

r1 = records[0]
smiles1 = r1["smiles"]

mol1 = Chem.MolFromSmiles(smiles1)

from rdkit.Chem import Draw
img = Draw.MolToImage(mol1)
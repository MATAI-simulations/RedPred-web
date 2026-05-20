from rdkit import Chem
from electrolyte_db.db import * 


query = {"n_atoms": 12}
results = query_molecule("pubchem_molecules", query_params=query, result_cols=['compound_cid', 'iupac_traditional_name', 'molecular_formula', 'smiles'])
print(results)


# # Molecules can be converted to and from text using Python’s pickling machinery:
# m = Chem.MolFromSmiles('c1ccncc1')

# import pickle

# pkl = pickle.dumps(m)
# m2=pickle.loads(pkl)

# smiles2 = Chem.MolToSmiles(m2)
# print(smiles2)




# from electrolyte_simulator.database.db import * 
# import pandas as pd
# from pathlib import Path




# def generate_csv(results, dir='.'):
#     data = {
#     "cid": [],
#     "name": [],
#     "molecular_formula": [],
#     "reactant_smiles": [],
# }

#     for res in results:
#         data['cid'].append(res['compound_cid'])
#         data['name'].append(res['iupac_traditional_name'])
#         data['molecular_formula'].append(res['molecular_formula'])
#         data['reactant_smiles'].append(res['smiles'])


#     df = pd.DataFrame(data)
#     df.to_csv(f'{dir}/reactant_smiles.csv')




# dir = Path(__file__).resolve().parent

# generate_csv(results, dir)
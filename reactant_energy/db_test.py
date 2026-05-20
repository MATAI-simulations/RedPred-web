from electrolyte_simulator.database.db import *

print(get_all_tables())
print(get_all_table_columns("pubchem_molecules"))

print(get_all_table_columns("calculations"))


# 1. Upsert all quinone derivatives from PubChem to the database 

# quinone = "C1=CC(=O)C=CC1=O"
# upsert_all_substructures(substructure=quinone)


# 2. Upsert a specific pubchem CID identifier 
conn = get_connection()
upsert_molecule_by_cid(conn, 8718)
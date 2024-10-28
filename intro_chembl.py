# import statements
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule  # initiates the query engine for molecules
# find molecule using the exact pref_name:
mols_pref = molecule.filter(pref_name__iexact='aspirin')
# or look in its synonyms:
mols_syn = molecule.filter(molecule_synonyms__molecule_synonym__icontains='Acetylsalicylic Acid')
print(mols_syn[0].keys())  # to see a list of the keys in the dictionary using the first entry
# now we can iterate through the results and print some details
for mol in mols_syn:
       print(f"{mol['molecule_chembl_id']} has {len(mol['molecule_synonyms'])} synonyms")

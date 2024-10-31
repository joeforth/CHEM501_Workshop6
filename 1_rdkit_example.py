# import statements
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors

# define a SMILES string for 3-aminobenzoic acid
smiles = 'Nc1cccc(c1)C(=O)O'
# convert SMILES to mol object
mol = Chem.MolFromSmiles(smiles)
# convert mol object to InChI
inchi = Chem.MolToInchi(mol)
# and InChI-Key
inchi_key = Chem.InchiToInchiKey(inchi)
# print the SMILES and InChI strings
print(f"SMILES: {smiles}\nInChI: {inchi}\nInChI-Key={inchi_key}")
# generate a Morgan fingerprint for the molecule
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
# print the fingerprint as a binary string
print(f'Fingerprint: {fp.ToBitString()}')
# calculate some properties for the molecule
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
hbd = Descriptors.NumHDonors(mol)
hba = Descriptors.NumHAcceptors(mol)

print(f"Molecular weight: {mw}\nLogP: {logp}", 
      f"\nNumber of hydrogen bond donors: {hbd}", 
      f"\nNumber of hydrogen bond acceptors: {hba}")
# generate a 2D depiction of the molecule
img = Draw.MolToImage(mol)
# display the image
img.show()

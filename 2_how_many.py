# first, we import the required modules
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, AllChem
from PIL import Image, ImageDraw, ImageFont

smiles_list = ["Nc1cccc(c1)C(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "CCO", "Nc1cccc(c1)C(=O)O",
               "CC(C)O", "CC(=O)O", "CCN", "COC", "CC(C)(C)O", "CC(C)C(=O)O", 
               "OC(=O)C(Br)(Cl)N", "ClC(Br)(N)C(=O)O"]

# Convert the SMILES strings into RDKit molecule objects:
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
# First, unique by SMILES, as written:
unique_smiles = set(smiles_list)
# Now we can use the RDKit algorithm to calculate canonical SMILES:
unique_canonical_smiles = set(Chem.MolToSmiles(mol, canonical=True) for mol in molecules if mol is not None)
# We can also calculate InChI 
unique_inchi = set(Chem.MolToInchi(mol) for mol in molecules if mol is not None)
# ... and InChIKeys
unique_inchi_key = set(Chem.MolToInchiKey(mol) for mol in molecules if mol is not None)

# What do you notice when you print the length of the sets?
print(f"Number of unique molecules by Non-Canonical SMILES: {len(unique_smiles)}")
print(f"Number of unique molecules by Canonical SMILES: {len(unique_canonical_smiles)}")
print(f"Number of unique molecules by InChI: {len(unique_inchi)}")
print(f"Number of unique molecules by InChIKey: {len(unique_inchi_key)}")

# Using fingerprints, we can also calculate the similarity of our molecules:
# Calculate fingerprints for all molecules
fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in molecules if mol is not None]
# Calculate similarity of one molecule (e.g., the first one in the list) against all others using Tanimoto coefficient
similarity_scores = [DataStructs.TanimotoSimilarity(fps[0], fp) for fp in fps]
print(similarity_scores)

# Finally, we can create the grid image with the molecules
img = Draw.MolsToGridImage(molecules, molsPerRow=4, subImgSize=(200,210))
# Create a drawing context
draw = ImageDraw.Draw(img)
# Specify the font to use for the text
font = ImageFont.truetype("arial.ttf", 12)
# Add scores as text to the image, centered horizontally and positioned below each molecule
for i, score in enumerate(similarity_scores):
    col = i % 4
    x = col * 200 + 100  # Center X coordinate
    row = i // 4
    y = (row + 1) * 200 + 10  # Position Y coordinate below the molecule
    draw.text((x, y), f"Score: {score:.2f}", font=font, fill=(0, 0, 0), anchor="mm")
# Display the image
img.show()

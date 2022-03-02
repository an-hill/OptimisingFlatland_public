import os, re, h5py
import pandas as pd

filename_and_smiles = []

with open('edges/mutated_mols.txt') as file:
    count = 1
    for line in file:
        smiles = line.strip()
        filename = "edge" + str(count)
        smiles_with_coo = re.sub(r'(\[Xe\])', r'(C(O)=O)', smiles)
        data = [filename, smiles_with_coo]
        filename_and_smiles.append(data)
        count += 1
        os.system(f'obabel -:"{smiles}" -O edges/edges_cml/{filename}.cml --gen3d')
        print(f"{filename} converted to .cml")

df = pd.DataFrame(filename_and_smiles, columns=['filename', 'smiles_with_coo'])
df.to_hdf('edge_smiles_database.h5', 'data')
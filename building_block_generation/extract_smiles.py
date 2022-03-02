import os, re, h5py
import pandas as pd
from openbabel import openbabel as ob

# THIS SCRIPT WORKS ONLY ON .CML FILES THAT DO NOT CONTAIN THE -CO(=O) IN THE SBU.
# TO USE THIS, ENSURE A SET OF SEPARATE .CMLS ARE USED FROM THE ONES BEING SENT TO 
# TOBACCO FOR CIF GENERATION.

path = 'nodes/smiles_extraction'

obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("cml", "smi")

filename_and_smiles = []

for f in os.scandir(path):    
    # Initialize an OBMol object
    mol = ob.OBMol()
    read_ok = obConversion.ReadFile(mol, f.path)
    if not read_ok:
        # There was an error reading the file
        raise Exception(f'Could not read file {f.path}')
    
    filename = f.name[:-4]
    smiles = obConversion.WriteString(mol).strip()
    
    data = [filename, smiles]
    filename_and_smiles.append(data)

df = pd.DataFrame(filename_and_smiles, columns=['filename', 'smiles'])
df.to_hdf('node_smiles_database.h5', 'data')
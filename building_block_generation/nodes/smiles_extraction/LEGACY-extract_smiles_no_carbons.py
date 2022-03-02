import os, re, h5py
import regex
import pandas as pd
from openbabel import openbabel as ob

# THIS SCRIPT WILL SCRUB ALL INSTANCES OF -CO(=O) FROM THE SMILES STRING GENERATED
# IT MIGHT BREAK IN SOME INSTANCES, IT WOULD BE BETTER TO CREATE A .CML WITH NO
# CARBONYL FUNCTIONAL GROUPS AND EXTRACT THE SMILES USING THE DEFAULT SCRIPT.

path = 'with_coo'

obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("cml", "smi")

for f in os.scandir(path):    
    # Initialize an OBMol object
    mol = ob.OBMol()
    read_ok = obConversion.ReadFile(mol, f.path)
    if not read_ok:
        # There was an error reading the file
        raise Exception(f'Could not read file {f.path}')

    smiles_string = obConversion.WriteString(mol)
    print(smiles_string)

    #(\(*\[O\]=CO\d*\)*)
    #(\(*OC=\[O\]\d*\)*)
    #(?<=^\[\w+\W*\])\d+ finds numbers at the start of the smiles string

    #smiles_withour_coo = re.sub(r'CO\(=O\)', r'', smiles_string)
    #regex = r'\[Xe\]'
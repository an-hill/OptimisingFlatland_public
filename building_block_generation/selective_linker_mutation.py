from rdkit import Chem
from rdkit.Chem import AllChem
import os

with open('edges/olinker_backbones_symmetry_limited.txt', 'r') as file:
    backbones=file.read().splitlines()

with open('edges/mutating_fgs.txt', 'r') as file:
    mutators=file.read().splitlines()

mutated_mols = []
rxn = AllChem.ReactionFromSmarts("[Rb][*:1].[Rb][*:2]>>[*:1][*:2]")

for mol in backbones:
    temp_mutated = []
    num_rb = mol.count("[Rb]")
    #print(num_rb)
    if "[Rb]" in mol:
        for mut in mutators:
            print('Mutating ' + mol + ' with ' + mut)
            to_mutate = []
            to_mutate.append(mol)
            for i in range(num_rb):
                for candidate in to_mutate:
                    backbone = Chem.MolFromSmiles(candidate)
                    mutator = Chem.MolFromSmiles(mut)
                    results = rxn.RunReactants( [backbone, mutator] )
                    for products in results:
                        for m in products:
                            smiles = Chem.MolToSmiles(m)
                            if "[Rb]" not in smiles:
                                temp_mutated.append(smiles)
                            elif smiles not in to_mutate:
                                to_mutate.append(smiles)
    else:
        print('Not mutating ' + mol + ', no [Rb] present.')
        temp_mutated.append(mol)
    for mutated in temp_mutated:
        if mutated not in mutated_mols:
            mutated_mols.append(mutated)

for mutated in mutated_mols:
    with open('mutated_mols.txt', 'a+') as myfile:
        myfile.write(mutated + '\n')

#print(mutated_mols)
print("appended " + str(len(mutated_mols)) + " molecules to file")
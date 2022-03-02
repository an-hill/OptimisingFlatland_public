from openbabel import openbabel
from openbabel import pybel
from openeye import oechem
import os
import pandas as pd
import re
import datetime

path = 'nodes/nodes_cml'

# Initialize the OpenBabel object that we will use later.
obConversion = openbabel.OBConversion()
obConversion.SetInFormat("cml")
ob = pybel.ob

get_unit_cell_from_file = False

for f in os.scandir(path):
    # Define containers
    bond_atom_0 = []
    bond_atom_1 = []
    bond_order = []
    bond_length = []
    atom_label = []
    x_coord = []
    y_coord = []
    z_coord = []
    element = []
    
    # Initialize an OBMol object
    mol = ob.OBMol()
    read_ok = obConversion.ReadFile(mol, f.path)
    if not read_ok:
        # There was an error reading the file
        raise Exception(f'Could not read file {f.path}')
    
    mol_name = f.name[:-4]
    num_atoms = mol.NumAtoms()
    
    # Assign unit cell data
    if get_unit_cell_from_file:
        molecule = pybel.Molecule(mol)
        unit_cell = molecule.unitcell
        a = unit_cell.GetA()
        b = unit_cell.GetB()
        c = unit_cell.GetC()
        alpha = unit_cell.GetAlpha()
        beta = unit_cell.GetBeta()
        gamma = unit_cell.GetGamma()
    else:
        a = 20
        b = 20
        c = 20
        alpha = 90
        beta = 90
        gamma = 90
    
    # Extracting atom labels and fractional coordinates
    for i in range(1, num_atoms + 1):
        atom = mol.GetAtom(i)
        # Divide by unit cell parameter to recover fractional coords
        # This is only doable simply because it is a cubic unit cell
        x_coord.append(round(atom.GetX() / a, 6))
        y_coord.append(round(atom.GetY() / b, 6))
        z_coord.append(round(atom.GetZ() / c, 6))
        atom_num = atom.GetAtomicNum()
        atom_symb = oechem.OEGetAtomicSymbol(atom_num)
        element.append(atom_symb)
        atom_id = atom_symb + str(i)
        atom_label.append(atom_id)
        
    coordinate_df = pd.DataFrame({'atom_label': atom_label,
                                  'element': element,
                                  'fract_x_coord': x_coord,
                                  'fract_y_coord': y_coord,
                                  'fract_z_coord': z_coord,
                                 })
    
        
    # Extract bonding information
    mol_bonds = openbabel.OBMolBondIter(mol) # iterate over all the bonds in the molecule
    for bond in mol_bonds:
        bond_atom_0.append(bond.GetBeginAtomIdx() - 1) # Must be 0-indexed
        bond_atom_1.append(bond.GetEndAtomIdx() - 1)
        bond_length.append(bond.GetLength())
        bond_order.append(bond.GetBondOrder())
        
    # Put bonding data into a dataframe
    bonding_df = pd.DataFrame({'atom_0': bond_atom_0,
                               'atom_1': bond_atom_1,
                               'length': bond_length,
                               'symmetry': '.',
                               'order': bond_order
                              })
    
    bonding_df = bonding_df.sort_values(['atom_0', 'atom_1']).reset_index(drop=True)    

    # Relabel bonding_df to work with ToBaCCo
    for i in range(0, num_atoms):
        bonding_df['atom_0'] = bonding_df['atom_0'].replace(i, atom_label[i])
        bonding_df['atom_1'] = bonding_df['atom_1'].replace(i, atom_label[i])
        bonding_df['order'] = bonding_df['order'].replace(1, 'S')
        bonding_df['order'] = bonding_df['order'].replace(2, 'D')
        bonding_df['order'] = bonding_df['order'].replace(3, 'T')

    # Checks to see if carbon is bonded to Xe, if true change to X
    num_bonds = len(bonding_df)
    y=[]
    z=[]
    for i in range(num_bonds):
        for j in range(2):
            if re.search(r"Xe\d", bonding_df.iloc[i]['atom_' + str(j)]):
                x = bonding_df.iloc[i]['atom_' + str(1-j)]
                y.append(x)
                z.append(x.replace('C','X'))

    for j in range(len(y)):
        for i in range(num_atoms):
            if re.search(y[j], coordinate_df.iloc[i]['atom_label']):
                coordinate_df.at[i, 'atom_label'] = z[j]
        for i in range(num_bonds):
            for k in range(2):
                if re.search(y[j], bonding_df.iloc[i]['atom_' + str(k)]):
                    bonding_df.at[i, 'atom_' + str(k)] = z[j]
    
    # Drop rows with Xe anywhere in them
    relabeled_coordinate_df = coordinate_df[~coordinate_df.element.str.contains('Xe')]
    temp_bonding_df = bonding_df[~bonding_df.atom_0.str.contains(r"Xe\d")]
    relabeled_bonding_df = temp_bonding_df[~temp_bonding_df.atom_1.str.contains(r"Xe\d")]

    with open('nodes/nodes_cif/%s.cif' % mol_name, 'w') as out:
        out.write('data_' + mol_name + '\n')
        out.write('_audit_creation_date              ' + datetime.datetime.today().strftime('%Y-%m-%d') + '\n')
        out.write("_audit_creation_method            'adam'" + '\n')
        out.write("_symmetry_space_group_name_H-M    'P1'" + '\n')
        out.write('_symmetry_Int_Tables_number       1' + '\n')
        out.write('_symmetry_cell_setting            triclinic' + '\n')
        out.write('loop_' + '\n')
        out.write('_symmetry_equiv_pos_as_xyz' + '\n')
        out.write('  x,y,z' + '\n')
        out.write('_cell_length_a                    ' + str(a) + '\n')
        out.write('_cell_length_b                    ' + str(b) + '\n')
        out.write('_cell_length_c                    ' + str(c) + '\n')
        out.write('_cell_angle_alpha                 ' + str(alpha) + '\n')
        out.write('_cell_angle_beta                  ' + str(beta) + '\n')
        out.write('_cell_angle_gamma                 ' + str(gamma) + '\n')
        out.write('loop_' + '\n')
        out.write('_atom_site_label' + '\n')
        out.write('_atom_site_type_symbol' + '\n')
        out.write('_atom_site_fract_x' + '\n')
        out.write('_atom_site_fract_y' + '\n')
        out.write('_atom_site_fract_z' + '\n')
        out.write(relabeled_coordinate_df.to_csv(sep=' ',index=False,header=False))
        out.write('loop_' + '\n')
        out.write('_geom_bond_atom_site_label_1' + '\n')
        out.write('_geom_bond_atom_site_label_2' + '\n')
        out.write('_geom_bond_distance' + '\n')
        out.write('_geom_bond_site_symmetry_2' + '\n')
        out.write('_ccdc_geom_bond_type' + '\n')
        out.write(relabeled_bonding_df.to_csv(sep=' ',index=False,header=False))

        print(f"{mol_name} converted to cif")
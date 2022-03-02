import os, random, re, h5py
import pandas as pd
#from mofid.run_mofid import cif2mofid

cif_with_data = []
count = 1

# For testing, ten mofs are selected
#random_selection_of_mofs = [random.choice(os.listdir('generated_mofs')) for i in range(10)]

node_df = pd.read_hdf('node_smiles_database.h5')
edge_df = pd.read_hdf('edge_smiles_database.h5')

for mof in os.listdir("tobacco-3.0/output_cifs"):
#for mof in random_selection_of_mofs:
    path_to_cif = os.path.join("generated_mofs", mof)
    print("Processing " + path_to_cif + f'({count})')
    count += 1

    ##### Manual MOFid generation #####
    topology = re.findall(r'^[a-z]{3,4}', mof)
    nodes = re.findall(r'(?<=v\d-)\d[a-zA-Z_]*(?=_v|_edge)', mof)
    edges = re.findall(r'edge_\d+', mof)

    user_comment = 'comment'

    nodes_as_smiles = []
    edges_as_smiles = []
    for i in range(len(node_df)):
        for j in range(len(nodes)):
            if node_df['filename'][i] == nodes[j]:
                nodes_as_smiles.append(node_df['smiles'][i])
    for i in range(len(edge_df)):
        for j in range(len(edges)):
            if edge_df['filename'][i] == edges[j]:
                edges_as_smiles.append(edge_df['smiles_with_coo'][i])

    mofid = '.'.join(nodes_as_smiles) + '.'
    mofid = mofid + '.'.join(edges_as_smiles) + ' '
    mofid = mofid + 'MOFid-v1' + '.'
    mofid = mofid + topology[0] + '.'  
    mofid = mofid + 'cat0'
    mofid = mofid + f';{user_comment}'

    data = [mof, mofid]
    # CURRENTLY NO MOFKEY IS GENERATED
    #data = [mof, mofid, mofkey]
    cif_with_data.append(data)

    ##### Automatic MOFid generation using MOFid library #####
    # generates a directory, can't just extract MOFid from the function
    #mofid_run = cif2mofid(path_to_cif, output_path=f'temp_mofids/{cif}')
    #with open(f'temp_mofids/{cif}/python_mofid.txt') as file:
    #    mofid = file.readline().strip()
    #with open(f'temp_mofids/{cif}/python_mofkey.txt') as file:
    #    mofkey = file.readline().strip()
    # Once the MOFid and MOFkey have been extracted, remove the directory
    #shutil.rmtree('temp_mofids/')
    #os.system(f'rm -r temp_mofids/{cif}')
    #print("MOFid = " + mofid)
    #print("MOFkey = " + mofkey)
    #data = [cif, mofid, mofkey]
   
df = pd.DataFrame(cif_with_data, columns=['cif_name', 'MOFid'])
print(df)
df.to_hdf('mofid_database.h5', 'data')

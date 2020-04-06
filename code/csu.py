import sys, os
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import *


import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)





def run_csu(input):
    """
    run csu
    """

    pdbID   = input[0]
    resID   = input[1]
    chainID = input[2]

    interactions = pd.DataFrame(columns=['Pos1', 'Pos2', 'AA', 'distance', 'surface', 'HB', 'Arom', 'Phob', 'DC'])

    proc = subprocess.Popen(['./resc', pdbID, resID, chainID], stdout=subprocess.PIPE)
    output = proc.stdout.read()

    output  = output.decode()

    counter = 0

    record = False
    for line in output.split('\n'):
        if line == "Table II":
            record = True

        elif line == "Table III":
            record = False

        if record:
            current_line = line.split()
            if len(current_line) == 8 and "Residues" not in current_line:
            
                pos    = current_line[0]
                AA     = str(current_line[1]) 
                dist   = float(current_line[2])
                surf   = float(current_line[3])
                HB     = str(current_line[4])
                Arom   = str(current_line[5])
                Phob   = str(current_line[6])
                DC     = str(current_line[7])

                interactions.loc[counter] = [resID, pos[:-1], AA, dist, surf, HB, Arom, Phob, DC ]
                counter += 1

    return interactions



def remove_hetatom(model):
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                residue_to_remove.append((chain.id, residue.id))
        if len(chain) == 0:
            chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return model



def nonneg_coords(model):

    rotation_matrix = rotaxis2m(0, Vector(1,0,0))
    translation_matrix = np.array((1000, 1000, 1000), 'f')

    for chain in model:
        for residue in chain:
            for atom in residue:
                atom.transform(rotation_matrix, translation_matrix)
            print(residue['CA'].get_vector() )

    return model


def get_network(pdb, outfile):
    """
    get hydrophobic contact network for pdb
    """

    G = nx.Graph()
    
    p = PDBParser()
    structure = p.get_structure('pdb', pdb)

    model = structure[0]    
    model = remove_hetatom(model)
    #model = nonneg_coords(model)

    chain_id = []
    for chain in model:
        chain_id.append(chain.get_id() )

    chain =  model[ chain_id[0] ] 

    lastres = None 
    for lastres in chain:
        pass

    L = len(chain) + lastres.get_id()[1]

    csu = np.zeros(( L, L ))

    for residue in chain:
        current_res = residue.get_id()[1] 

        inp = [ pdb, str(current_res), str( chain_id[0] ) ]
        I = run_csu( inp )

        for index, row in I.iterrows():
            if int(row['Pos1']) < int(row['Pos2']) and row['Phob'] == '+':
                G.add_edge( int(row['Pos1']), int(row['Pos2']) )
                csu[ int(row['Pos1']), int(row['Pos2']) ] = 1

    nx.write_gpickle(G, outfile+'.nx')


pdbIN = sys.argv[1]
pdbOUT = sys.argv[2]

get_network(pdbIN, pdbOUT)

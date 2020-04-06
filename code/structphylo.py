import sys, os
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import *
from Bio import SeqIO
import itertools as it


import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)


yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))

structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
current_structures = structures[structures['sf']=='c.37.1']
list_c371 = list(current_structures['ORF'])
    
    


#sf_counts = structures.groupby(['scop', 'substrate'] ).count()  

struct2 = structures[['sf', 'substrate', 'ORF']]

sf_counts = struct2.pivot_table(index='sf', columns='substrate', aggfunc=len).fillna(0).astype('int')
sf_counts.to_csv("counts.txt", sep='\t', header=False, index=True)




def compute_distance(pdb1, pdb2):
    """
    computer Z score between two superposed PDB structures
    """

    print(pdb1, pdb2)

    proc = subprocess.Popen(['superpose', pdb1, pdb2], stdout=subprocess.PIPE)

    output = proc.stdout.read()
    output  = output.decode()

    Q = 0
    rmsd = 0
    Nalign = 0 

    sel = False

    for line in output.split('\n'):

        if 'Scores' in line:
            sel = True

        if line[0:2] == "$$":
            sel = False

        if sel: 
            if 'quality' in line:
                Q = line.split()[2]
            if 'r.m.s.d' in line:
                rmsd = line.split()[1]
            if 'Nalign' in line: 
                Nalign = line.split()[1]

    return Q, rmsd, Nalign

#Q, rmsd, N = compute_distance(pdbA, pdbB)



def compute_distmat(inplist):
    """
    computes Q score for pairwise PDB superposition
    """
    
    L = len(inplist)

    Qmat    = np.zeros((L,L))
    RMSDmat = np.zeros((L,L))
    Nmat    = np.zeros((L,L))

    family = []
    sbstr = []

    for orf in inplist:
        current_fam = structures[structures['ORF']==orf]['fam'].item()
        family.append(current_fam)
        current_sbstr = structures[structures['ORF']==orf]['substrate'].item()
        sbstr.append(current_sbstr)

    for i,j in it.combinations(inplist,2):
 
        i_ind = inplist.index(i)
        j_ind = inplist.index(j)

        pdbA = '../data/pdb/minimized/' + i + '.pdb'
        pdbB = '../data/pdb/minimized/' + j + '.pdb'

        Q, rmsd, N = compute_distance(pdbA, pdbB)
    
        Qmat[i_ind, j_ind] = Qmat[j_ind, i_ind] = float(Q)
        RMSDmat[i_ind, j_ind] = RMSDmat[j_ind, i_ind] = float(rmsd)
        Nmat[i_ind, j_ind] = Nmat[j_ind, i_ind] = float(N) 


    np.savetxt('../data/processed/Qmat_dist.txt', Qmat)

    print(Qmat)

    Qmat = np.around(Qmat, 4)
    
    Qdf = pd.DataFrame(data=Qmat, columns=list(inplist) )
    Qdf.insert(loc=0, column='ORF', value=list(inplist) )
    Qdf.insert(loc=1, column='fam', value=family)
    Qdf.insert(loc=2, column='substrate', value=sbstr)

    Qdf.to_csv("../data/processed/Q_DF.txt", sep='\t', header=True, index=False)

    #np.savetxt('RMSD.txt', RMSDmat)
    #np.savetxt('N.txt', Nmat)

compute_distmat(list_c371)




def get_sequences(orflist, seqs):
    """
    pull out sequences from yeast fasta file
    """

    records = []
    for i in orflist:
        #print(i)
        #if i in yeast_seqs:
        records.append( seqs[i] )

        SeqIO.write(records, "tmp.fa", "fasta")

    return records


get_sequences(list_c371, yeast_seqs)



def load_dssp(fname):
    """
    reads pre-computed DSSP file into dataframe
    """

    #: max accessible surface area (square angstroms) for amino acids, from
    #: `Tien et al (2013) <https://doi.org/10.1371/journal.pone.0080635>`_.
    #: `MAX_ASA_TIEN[a]` is the max surface area for amino acid `a`.
    MAX_ASA_TIEN = {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0,
                'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0,
                'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0,
                'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}


    dssp = pd.DataFrame(columns=['Residue', 'Chain', 'AA', 'ASA', 'RSA', 'SS', 'X', 'Y', 'Z'])

    counter = 0
    df_counter = 0
    with open(fname) as f:
    #    next(f)
        for line in f:
            counter += 1

            if counter > 28:                    # remove header of the DSSP output
                res     = int(line[6:10])
                chain   = str(line[11])
                AA      = str(line[13])
                ASA     = int(line[35:38])
                RSA     = np.around(ASA / float(MAX_ASA_TIEN[AA]), 2)
                SS      = str(line[16])
 
                if SS.strip() in ['G', 'H', 'I']:
                    SStype = 'helix'
                elif SS.strip() in ['B', 'E']:
                    SStype = 'strand'
                else:
                    SStype = 'loop'

                X       = float(line[115:122])
                Y       = float(line[123:130])
                Z       = float(line[131:138])
            
                dssp.loc[df_counter] = [res, chain, AA, ASA, RSA, SStype, X, Y, Z]
                df_counter += 1


    return dssp



def load_data(structDF):
    """
    load DSSP and NX data into two dataframes
    """

    dssp = {}
    
    for i in structDF['ORF']:
        dssp[i]    = load_dssp('../data/pdb/dssp/'+i+'.dssp')
       
    return dssp




def dssp_mat(ORFs):
    """
    parse DSSP data into matrix
    """

    DSSP = load_data(structures)

    dssp_dict = {}

    maxL = 0

    for orf in ORFs:
        current_dssp = DSSP[orf]['SS']

        current_dssp = np.array(current_dssp)
        current_dssp[current_dssp=="loop"] = 0
        current_dssp[current_dssp=="helix"] = 1
        current_dssp[current_dssp=="strand"] = 2
        current_dssp = np.array(current_dssp, dtype=int)

        dssp_dict[orf] = current_dssp

        if len(current_dssp) > maxL:
            maxL = len(current_dssp)
    
    dsspMAT = np.ones((len(ORFs), maxL) ) * np.nan

    for ix, i in enumerate( dssp_dict.keys() ):
        dsspMAT[ix,0:len(dssp_dict[i])] = dssp_dict[i]

    dssp_df = pd.DataFrame(data=dsspMAT )
    dssp_df.insert(loc=0, column='ORF', value=list(ORFs) )



    return dssp_df


dssp_c371 = dssp_mat(list_c371)
dssp_c371.to_csv("../data/processed/dssp_c371.txt", sep='\t', header=True, index=False)



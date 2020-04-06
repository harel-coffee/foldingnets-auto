import sys, os
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import *
from Bio import SeqIO
from Bio import AlignIO
import itertools as it
import random

from sklearn.metrics import roc_auc_score
import multiprocessing as mp



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_network(fname):
    """
    read in pre-computed pickeled networkx graph
    """
    G = nx.read_gpickle(fname)

    return G



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_dssp(fname):
    """
    reads pre-computed DSSP file into dataframe
    #: max accessible surface area (square angstroms) for amino acids, from
    #: `Tien et al (2013) <https://doi.org/10.1371/journal.pone.0080635>`_.
    #: `MAX_ASA_TIEN[a]` is the max surface area for amino acid `a`.
    """
    MAX_ASA_TIEN = {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0,
                'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0,
                'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0,
                'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}

    dssp = pd.DataFrame(columns=['Residue', 'Chain', 'AA', 'ASA', 'RSA', 'SS', 'X', 'Y', 'Z'])

    counter = 0
    df_counter = 0
    with open(fname) as f:
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_data(structDF):
    """
    load DSSP and NX data into two dataframes
    """
    dssp = {}
    contmat = {}

    for i in structDF['ORF']:
        dssp[i]    = load_dssp('../data/pdb/dssp/'+i+'.dssp')
        contmat[i] = load_network('../data/pdb/contmat/'+i+'.nx')

    return dssp, contmat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map_nopt(orflist, theta_contact=2, theta_dist=10, randomized=False):
    """
    get clusters of nonoptimal codons associated with contacts
    """

    noptDF = pd.DataFrame(columns=['ORF', 'substrate', 'nopt', 'Cnopt', 'CC'])
    counter = 0

    for i in orflist:
        substr = structures[structures['ORF']==i]['substrate'].item()
        current_cm = CONTMAT[i]

        # create directed graph, direction in sequence order
        cm2 = nx.DiGraph()
        for e in current_cm.edges():
            e1 = e[0]
            e2 = e[1]
            if e1 < e2:
                cm2.add_edge(e1,e2)
            else:
                cm2.add_edge(e2,e1)

        ddist_cm2 = dict( cm2.in_degree() ) 
        degin = []
        for k in ddist_cm2.keys():
            if ddist_cm2[k] > theta_contact:
                degin.append(k-1)                           # this is a fix for python indexing!

        if substr == '-':
            current_substr = 0
        elif substr == '+':
            current_substr = 1

        if i in rna_seqs.keys():
            rna_seq = rna_seqs[i].seq
            L = int( np.floor( len(rna_seq) / 3. ) )
            profile = np.zeros(( L ))

            for pos in range(L):
                start = 3*int(pos)
                end   = start + 3
                current_codon = rna_seq[start:end]
                if randomized: 
                    current_codon = rand_syn(current_codon)
                if current_codon in nopt:
                    profile[pos] = 1

        # check for cluster,w=window_size, theta=threshold
        w = 4
        theta = 3
        cluster_profile = np.zeros(( L ))
        for pos in range(L-w):
            current_window = profile[pos:(pos+w)]
            if np.sum(current_window) >= theta:
                cluster_profile[pos] = 1

        CC = 0        
        current_contacts = np.array(degin)
        current_clusters = np.where(cluster_profile == 1)[0]
        for pos in current_clusters:
            current_dists = ( current_contacts - 35 ) - pos
            current_dists = current_dists[current_dists >= 0]
            if any( current_dists < theta_dist ):
                CC += 1
  
        noptDF.loc[len(noptDF)] = (i, current_substr,  np.sum(profile)/len(profile), np.sum(cluster_profile)/len(profile), CC )


    return  noptDF





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map_nopt_rand(packaged_inputs):
    """
    get clusters of nonoptimal codons associated with contacts
    duplicate function to handle packaged input for multiprocessing (careful with changes!)
    """

    structureDF = packaged_inputs[0]
    CONTMAT     = packaged_inputs[1]
    rna_seqs    = packaged_inputs[2]
    nopt        = packaged_inputs[3]
    theta_contact = packaged_inputs[4]
    theta_dist  = packaged_inputs[5]

    N = len(structures)
  
    nopt_mat  = np.zeros(( N ))
    Cnopt_mat = np.zeros(( N ))
    CC_mat    = np.zeros(( N ))

    for ix, i in enumerate(structureDF['ORF']):
        substr = structureDF[structureDF['ORF']==i]['substrate'].item()
        current_cm = CONTMAT[i]

        if substr == '-':
            current_substr = 0
        elif substr == '+':
            current_substr = 1

        # create directed graph, direction in sequence order
        cm2 = nx.DiGraph()
        for e in current_cm.edges():
            e1 = e[0]
            e2 = e[1]
            if e1 < e2:
                cm2.add_edge(e1,e2)
            else:
                cm2.add_edge(e2,e1)

        ddist_cm2 = dict( cm2.in_degree() ) 
        degin = []
        for k in ddist_cm2.keys():
            if ddist_cm2[k] > theta_contact:
                degin.append(k-1)                     # this is a fix for python indexing!


        if i in rna_seqs.keys():
            rna_seq = rna_seqs[i].seq
            L = int( np.floor( len(rna_seq) / 3. ) )
            profile = np.zeros(( L ))

            for pos in range(L):
                start = 3*int(pos)
                end   = start + 3
                current_codon = rna_seq[start:end]
                current_codon = rand_syn(current_codon)
                if current_codon in nopt:
                    profile[pos] = 1

        # check for cluster,w=window_size, theta=threshold
        w = 4
        theta = 3
        cluster_profile = np.zeros(( L ))
        for pos in range(L-w):
            current_window = profile[pos:(pos+w)]
            if np.sum(current_window) >= theta:
                cluster_profile[pos] = 1

        CC = 0        
        current_contacts = np.array(degin)
        current_clusters = np.where(cluster_profile == 1)[0]
        for pos in current_clusters:
            current_dists = ( current_contacts - 35 ) - pos             # 35 for ribosome tunnel
            current_dists = current_dists[current_dists >= 0]
            if any( current_dists < theta_dist ):
                CC += 1
  
        nopt_mat[ix] = np.sum(profile)/len(profile)
        Cnopt_mat[ix] = np.sum(cluster_profile)/len(cluster_profile)
        CC_mat[ix] = CC

    return  np.around(nopt_mat,3), np.around(Cnopt_mat,3), CC_mat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rand_syn(cdn):
    """
    get random synonymous codon
    """
    current_aa  = gencode[gencode['codon']==cdn]['aa'].item()
    current_syn = list(gencode[gencode['aa']==current_aa]['codon'])

    rand_syn = random.choice(current_syn)

    return rand_syn



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_rand_nopt(structures, CONTMAT, rna_seqs, nopt, theta_contact, theta_dist):
    """
    randomize codon sequences
    evaluate number of nopt, nopt-clusters, and nopt-clusters near high contacts
    """

    packaged_inputs = (structures, CONTMAT, rna_seqs, nopt, theta_contact, theta_dist)

    nrand = 10000


    pool = mp.Pool(processes=10)
    pool_inputs = [packaged_inputs for i in range(nrand)]
    rand_result = pool.map(map_nopt_rand, pool_inputs)
    pool.close()


    L = len(rand_result[0][0])

    rand_nopt = np.zeros(( L, nrand ))
    rand_Cnopt = np.zeros(( L, nrand ))
    rand_CC = np.zeros(( L, nrand ))

    for i in range(nrand):
        rand_nopt[:,i] = rand_result[i][0]
        rand_Cnopt[:,i] = rand_result[i][1]
        rand_CC[:,i] = rand_result[i][2]

    np.savetxt("../data/processed/rarecodon_rand_nopt_10000.txt", np.around(rand_nopt, 4), fmt='%.4f' )
    np.savetxt("../data/processed/rarecodon_rand_Cnopt_10000.txt", np.around(rand_Cnopt,4 ), fmt='%.4f')
    np.savetxt("../data/processed/rarecodon_rand_CC_10000.txt", rand_CC, fmt='%i')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_heatmap_CC():
    """
    scout parameters
    """
    resultDF = pd.DataFrame()
    for i in range(5):                          # theta contacts
        for j in range(5, 15):                  # theta distance
            colname = str("par_")+str(i)+str("_")+str(j)
            noptdf = map_nopt(structures['ORF'], i, j, randomized=False)
            resultDF[colname] = noptdf['CC']
    resultDF.to_csv("../data/processed/rarecodon_grid_CC.txt", sep='\t', header=True, index=False)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def expressionDF(orflist, yeast_pax):
    """
    parse substrate and protein abundance data into DF
    """

    abundanceDF = pd.DataFrame(columns=['ORF', 'substrate', 'abundance'])

    for i in orflist:
        substr = structures[structures['ORF']==i]['substrate'].item()
        if substr == '-':
            current_substr = 0
        elif substr == '+':
            current_substr = 1

        current_abundance = yeast_pax[yeast_pax['ORF']==i]['abundance'].item()
        abundanceDF.loc[len(abundanceDF)] = (i, current_substr, current_abundance)
    
    abundanceDF.to_csv("../data/processed/yeast_abundance.txt", sep='\t', header=True, index=False)






if __name__ == "__main__":

    # load data
    gencode = pd.DataFrame(
        {'codon':['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG', 
        'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG',
        'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 
        'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG'], 
        'aa':['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', '*', 'W', 
        'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 
        'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 
        'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'] })

    rna_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/orf_coding.fasta", "fasta"))
    structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)

    # bottom 20% of tAI scale
    nopt = ['CGA', 'CTT', 'AGT', 'CTG', 'ATA', 'TTT', 'AAT', 'CCT', 'TAT', 'GTA', 'TCA', 'CTC']

    pax = pd.read_csv("../data/yeast/yeast_paxdb.txt", header=11, sep='\t', index_col=False)
    pax['ORF'] = [entry.split('.')[1] for entry in list(pax['string_external_id'])]
    pax = pax[['ORF', 'abundance']]
    expressionDF(structures['ORF'], pax)


    DSSP, CONTMAT = load_data(structures)


    # run analyses
    get_rand_nopt(structures, CONTMAT, rna_seqs, nopt, 2, 10) 
    get_heatmap_CC() # lazy coding, should pass input args


    noptdf = map_nopt(structures['ORF'], 2, 10, randomized=False)
    noptdf.to_csv("../data/processed/nopt_clusters.txt", sep='\t', header=True, index=False)

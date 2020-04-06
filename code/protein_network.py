import sys, os
import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import *
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
import itertools as it
from sklearn import metrics




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

#coord = dssp[dssp['Residue'] == 3 ][['X', 'Y', 'Z']]
#print( dssp['SS'] )



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
def project_CM(alnfile, contmat):
    """
    project contact maps onto alignment
    """

    from Bio.Align import MultipleSeqAlignment

    alignment = AlignIO.read(open(alnfile), "clustal")

    # filter dssp entries from alignment file
    alignment_nr = MultipleSeqAlignment([])
    for ix, record in enumerate(alignment):
        if "dssp" not in record.id and "space" not in record.id:
            alignment_nr.append(record)
    alignment = alignment_nr

    nrow = len(alignment)
    ncol = alignment.get_alignment_length()

    mapping = {}
    list_id=[]

    CM = np.zeros(( ncol, ncol, nrow ))


    for ix, record in enumerate(alignment):
        seq_aln = np.array(record.seq)
        seq_ref = "".join(list(seq_aln[seq_aln!='-']))
        ident = record.id.split('_')[0]
        list_id.append(ident)
        current_substr = structures[structures['ORF'] == ident]['substrate'].item()

        current_mapping = {}
        # dict of aln position for seq position
        pos_aln = 0
        pos_ref = 1                                 # this is a FIX for python indexing vs. PDB
        while (pos_aln < ncol):
            if seq_aln[pos_aln] == '-':
                pos_aln += 1
            else:
                current_mapping[pos_ref] = pos_aln
                pos_aln += 1
                pos_ref += 1

        mapping[ident] = current_mapping

        current_cm = CONTMAT[ident]
        for n1, n2 in current_cm.edges():
            aln1 = current_mapping[n1]
            aln2 = current_mapping[n2]

            CM[aln1, aln2, ix] = CM[aln2, aln1, ix] = 1

    
    cons_cm = np.sum(CM, 2)
    np.savetxt("../data/processed/consensus_cm.txt", cons_cm, fmt='%i' )

    cm_df = pd.DataFrame(columns=['i', 'j', 'n'])
    counter = 0
    for i in range(ncol):
        for j in range(ncol):
            if cons_cm[i,j] > 0:
                cm_df.loc[counter] = [i, j, cons_cm[i,j]]
                counter += 1

    cm_df.to_csv("../data/processed/consensus_cm_df.txt", sep='\t', header=True, index=False)


    return CM, list_id





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_aligned_contacts():
    """
    buried cored from structure-based sequence alignment for c.37.1.8
    """

    
    def distances_buried_contacts(buriedcore, alnseq, cm):
        """
        subroutine to extract distances between buried contacts from aligned structural core
        buriedcore: alignment positions of buried structural core
        alnseq: aligned sequence with alignment gaps
        cm: contact map of protein
        """

        dists = []
        current_dist = 0
        contacts = np.where( np.sum(cm, 0) > 0 )[0]      # python indices of positions with contacts
        buried_contacts = np.intersect1d(buriedcore, contacts)

        for pos in range(len(alnseq)):
            if alnseq[pos] != "-":
                if pos in buried_contacts:
                    dists.append(current_dist)
                    current_dist = 0
                else:
                    current_dist += 1    
        
        dists = np.array(dists)

        return dists 



    asamat = pd.read_csv("../data/processed/c3718_aln_asa.txt", sep='\t', header=0)
    substr = list(asamat['substrate'])
    asamat.drop('substrate', 1, inplace=True)
    asamat = np.array(asamat)

    asa_mean = np.nanmean(asamat, 0)

    buried20  = np.where(asa_mean < 20)[0]
    buried30  = np.where(asa_mean < 30)[0]
    buried50  = np.where(asa_mean < 50)[0]

    exposed = np.where(asa_mean >= 50)[0]
    aligned = np.where( sum(np.isnan(asamat),0) < 12 )[0]           # at most 5 seqs not aligned at positions

    sel_buried_aligned20  = np.intersect1d(buried20, aligned)
    sel_buried_aligned30  = np.intersect1d(buried30, aligned)
    sel_buried_aligned50  = np.intersect1d(buried50, aligned)


    L = len(sel_buried_aligned20)
    N = len(asamat)

    histo = np.zeros((N, L))

    contacts = pd.DataFrame(columns=['ORF', 'substrate', 'mASA20', 'mASA30', 'mASA50', 'dists'])

    for i in range(16):
        current_cm = MAPCM[:,:,i]

        current_asa20 = current_cm[:,sel_buried_aligned20]
        current_asa30 = current_cm[:,sel_buried_aligned30]
        current_asa50 = current_cm[:,sel_buried_aligned50]

        ident = list_id[i]

        alignment = AlignIO.read(open('c3718.aln'), "clustal")
        current_seqrec = alignment[i]
        current_seq    = current_seqrec.seq

        histo[i,:] = np.sum(current_asa20,0)

        current_dists = distances_buried_contacts(sel_buried_aligned20, current_seq, current_cm)

        contacts.loc[len(contacts)] = ( ident, substr[i], np.mean(np.sum(current_asa20,0)), np.mean(np.sum(current_asa30,0)), np.mean(np.sum(current_asa50,0)), np.mean(current_dists) )



    contacts.to_csv("../data/processed/c3718_structaln_contacts.txt", sep='\t', header=True, index=False)

    histo = pd.DataFrame(histo)
    histo.insert(loc=0, column='substrate', value=substr)
    histo.to_csv("../data/processed/c3718_structaln_buried_degree.txt", sep='\t', header=True, index=False)

 
    return contacts





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def graph_to_matrix(G, nodelist, L):
    """
    nx graph to numpy matrix with full dimensions
    G: input graph
    nodelist: list of nodes to include
    L: length of protein
    THIS FIXES THE PYTHON VS. PDB INDEXING (SUBTRACTION BY 1)
    """

    cmatrix = np.zeros(( L, L ))
    #print(L, np.max(G.nodes() ))

    for edge in G.edges():
        if edge[0] in nodelist or edge[1] in nodelist:      # and
            if edge[0] < edge[1]:
                cmatrix[ edge[0]-1, edge[1]-1 ] = 1
            else:
                cmatrix[ edge[1]-1, edge[0]-1 ] = 1 

    return cmatrix





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_contacts(orflist):
    """
    get stats on buried/exposed hydrophobic contacts
    DSSP and PDB on same indexing! this function should be fine. 
    """

    contdens = pd.DataFrame(columns=['ORF', 'substrate', 'sf', 'fam', 'buriedCD', 'buried_alpha', 'buried_beta', 'exposedCD', 'buried_ldens', 'buried_wldens'])
    counter = 0

    for i in orflist:

        current_record = structures[structures['ORF'] == i]

        current_sf    = current_record['sf'].item()
        current_fam   = current_record['fam'].item()
        current_substr = current_record['substrate'].item()
        if current_substr == '-':
            current_substr = 0
        elif current_substr == '+':
            current_substr = 1

        current_cm    = CONTMAT[i]
        current_dssp  = DSSP[i]


        current_buried = current_dssp[current_dssp['ASA'] < 50]
        current_exposed = current_dssp[current_dssp['ASA'] >= 50]

        current_seq = "".join( list(current_dssp['AA']) )
        current_ddist  = dict( current_cm.degree() )

        current_buried_ix  = list(current_buried['Residue'])
        current_buried_alpha_ix = list( current_buried[current_buried['SS'] == 'helix']['Residue'])
        current_buried_beta_ix  = list( current_buried[current_buried['SS'] == 'strand']['Residue'])
        
        ddist_buried  = np.array([current_ddist.get(x, 0) for x in current_buried_ix] ) / 2.
        ddist_buried_alpha = np.array([current_ddist.get(x, 0) for x in current_buried_alpha_ix] ) / 2.
        ddist_buried_beta  = np.array([current_ddist.get(x, 0) for x in current_buried_beta_ix] ) / 2.

        current_buried_contact = np.intersect1d(current_buried_ix, np.array(list(current_cm.nodes() )))
        ldens_buried = np.array([ (current_buried_contact[pos+1]-current_buried_contact[pos]) for pos in np.arange( len(current_buried_contact)-1) ])
        ldens_buried_weighted = np.array([ (current_buried_ix[pos+1]-current_buried_ix[pos])*(current_ddist.get(pos,1) + current_ddist.get(pos+1,1))/2. for pos in np.arange( len(current_buried_ix)-1) ])
        
        current_cmat = graph_to_matrix(current_cm, current_buried_ix, len(current_seq) )
       
        current_exposed_ix = list(current_exposed['Residue'])
        ddist_exposed = np.array([current_ddist.get(x, 0) for x in current_exposed_ix] ) / 2.


        contact_density_buried = np.nansum(ddist_buried) / len(current_buried_ix)
        if len(current_buried_alpha_ix) > 0:
            contact_density_buried_alpha = np.nansum(ddist_buried_alpha) / len(current_buried_alpha_ix)
        else:
            contact_density_buried_alpha = 0

        if len(current_buried_beta_ix) > 0:
            contact_density_buried_beta  = np.nansum(ddist_buried_beta) / len(current_buried_beta_ix)
        else:
            contact_density_buried_beta = 0

        contact_density_exposed = np.nansum(ddist_exposed) / len(current_exposed_ix)
        

        if len(current_seq) == len(current_dssp):            
            contdens.loc[counter] = (i, current_substr, current_sf, current_fam, np.around(contact_density_buried, 3), np.around(contact_density_buried_alpha, 3), np.around(contact_density_buried_beta, 3), np.around(contact_density_exposed, 3), np.around(np.mean(ldens_buried),3), np.around(np.mean(ldens_buried_weighted),3) )
            counter += 1

    contdens.to_csv("../data/processed/all_contact_density.txt", sep='\t', header=True, index=False)

    return contdens






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def normalized_cd(contdensDF):
    """
    simple 'normalization' of contact density across protein families
    based on subtraction of best threshold
    """

    def evaluate_theta(substr, cd, theta):
        """
        subroutine to evaluate threshold
        substr and cd have to have same length
        """

        TP = 0; FP = 0; TN = 0; FN = 0

        pred = np.array( cd < theta, dtype=int)

        for i in range(len(substr)):
            if substr[i] == 1 and pred[i] == 1:
                TP += 1
            elif substr[i] == 0 and pred[i] == 1:
                FP += 1
            elif substr[i] == 0 and pred[i] == 0:
                TN += 1
            elif substr[i] == 1 and pred[i] == 0:
                FN += 1
        
        return TP, TN, FP, FN



    def compute_MCC(TP, TN, FP, FN):
        """
        compute Matthew's Correlation Coefficient
        """

        N = TN + TP + FN + FP
        if N > 0:
            S = (TP + FN) / float(N)
            P = (TP + FP) / float(N)
        else:
            S = np.nan
            P = np.nan

        if ( P*S*(1-S)*(1-P) ) > 0:
            MCC = ( (TP/float(N)) - (S*P) ) / np.sqrt( P*S*(1-S)*(1-P) ) 
        else:
            MCC = np.nan


        return MCC


    normCD = np.zeros(( len(contdensDF) ))
    normBD = np.zeros(( len(contdensDF) ))
    resultDF = pd.DataFrame(columns=['sf', 'theta', 'TP', 'TN', 'FP', 'FN'])

    for sf in list(set(contdensDF['sf'])):
        current_sf = contdensDF[contdensDF['sf']==sf]
        current_buriedCD = np.array( current_sf['buriedCD'] )
        current_substr = np.array( current_sf['substrate'] )
        current_BD = np.array(current_sf['buried_ldens'])

        current_result = np.zeros(( 200 ))
        for ix, i in enumerate( np.arange(200)/100. ) :
            TP, TN, FP, FN = evaluate_theta(current_substr, current_buriedCD, i)
            mcc = compute_MCC(TP, TN, FP, FN)
            current_result[ix] = mcc

        current_max = np.nanmax(current_result)
        current_best = np.where(current_result == current_max)[0]
        current_best_min = np.min(current_best)
        current_best_max = np.max(current_best)
        best_theta = np.mean([current_best_min, current_best_max]) / 100.
        TP, TN, FP, FN = evaluate_theta(current_substr, current_buriedCD, best_theta)       

        resultDF.loc[len(resultDF)] = (sf, best_theta, TP, TN, FP, FN)
        normCD[contdensDF['sf']==sf] = current_buriedCD - best_theta

        # dist
        current_result_dist = np.zeros(( 200 ) )
        for ix, i in enumerate( np.arange(200)/10. ) :
            TP2, TN2, FP2, FN2 = evaluate_theta(current_substr, current_BD, i)
            mcc2 = compute_MCC(TP2, TN2, FP2, FN2)
            current_result_dist[ix] = mcc2

        current_max2 = np.nanmax(current_result_dist)
        current_best2 = np.where(current_result_dist == current_max2)[0]
        current_best2_min = np.min(current_best2)
        current_best2_max = np.max(current_best2)
        best_theta2 = np.mean([current_best2_min, current_best2_max]) / 10.

        TP2, TN2, FP2, FN2 = evaluate_theta(current_substr, current_BD, best_theta2)       
        normBD[contdensDF['sf']==sf] = current_BD - best_theta2



    contdensDF['normCD'] = normCD
    contdensDF['normBD'] = normBD

    contdensDF.to_csv("../data/processed/all_contact_density_normalized.txt", sep='\t', header=True, index=False)
    resultDF.to_csv("../data/processed/superfamily_contact_density.txt", sep='\t', header=True, index=False)

    return contdensDF, resultDF





if __name__ == "__main__":


    # load data
    structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
    yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))

    # run analyses
    DSSP, CONTMAT = load_data(structures)
    MAPCM, list_id = project_CM('c3718.aln', CONTMAT)
    res_contacts = get_aligned_contacts()
    contdens = get_contacts(structures['ORF'])
    contdensN, sfstats = normalized_cd(contdens)


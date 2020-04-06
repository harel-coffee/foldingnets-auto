import os, sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from random import shuffle
import subprocess
from itertools import product
import weblogo
import shutil

import multiprocessing as mp

from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score

 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def peptide_matrix(fastaIN, w=7):
    """
    Generates matrix of counts of nr peptides
    FASTA record has to contain ORF and substrate key
    """

    N = 5000                    # preset matrix dimension, heuristic
    M_max = 1000000             # change to sparse matrix when with internet access to look up syntax

    pepmat = np.zeros(( N, M_max ), dtype=int)
    list_peps = []  
    sbstrmat = np.zeros(( N )) * np.nan 
    list_id = []
    id_idx = 0
    for record in fastaIN:
        current_seq = record.seq
        current_L = len(current_seq)
        current_id = record.id
        current_sbstr = int( record.description.split(' ')[1] ) 
        if current_id not in list_id:
            list_id.append(current_id)
            sbstrmat[id_idx] = current_sbstr
            id_idx += 1

        for j in range(current_L - w):
            current_pep = str( current_seq[j:(j+w)] )
            if current_pep not in list_peps:
                list_peps.append(current_pep)
            pepmat[list_id.index(current_id), list_peps.index(current_pep)] += 1

        if len(list_peps) > M_max - 1:
            print("exceeded preset matrix dimensions")
            break

    K = len(list_id)
    L = len(list_peps)
    pepmat = pepmat[0:K,0:L]
    sbstr = sbstrmat[0:K]

    return pepmat, list_peps, sbstr


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pairwise_score(pep1, pep2):
    """
    computes pairwise alignment score
    assumes peptides of same length
    default is BLOSUM62 for now
    """

    score = 0
    for i in range( min( len(pep1), len(pep2) ) ):
        score -= blosum62[pep1[i]][pep2[i]]

    return score


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pairwise_score_parallel( packaged_input ):
    """
    computes pairwise alignment score
    assumes peptides of same length
    default is BLOSUM62 for now
    """

    pep1 = packaged_input[0]
    pep2 = packaged_input[1]

    score = 0
    for i in range( min( len(pep1), len(pep2) ) ):
        score -= blosum62[pep1[i]][pep2[i]]

    return score


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pairwise_score_vectorized(pep, aln):
    """
    computes pairwise alignment score
    assumes peptide and aln of same length
    -> negative of BLOSUM score so that the 'distance' between seqs can be minimized
    """

    scores = np.zeros(( np.shape(aln) ))
    for i in range( len(pep) ):
        scores[:,i] = blosum62[pep[i]][aln[:,i]]

    scores = - np.sum(scores, 1)

    return scores


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def greedy_clustering(peps, min_size=5):
    """
    Performs greedy initial clustering of peptide seqs
    with threshold and BLOSUM62 similarity
    """
    print("Initial greedy clustering ...")

    index_table = np.zeros(( len(peps) ))
    theta = -11.9 # - len * 1.7  (default of Hammock paper)
    seeds = [0]
    seeds_peps = []
    seeds_peps.append(list(peps[0]))
    current_cluster = 0

    # check all seeds, join highest that is also above theta. else new seed
    for i in range( len(peps) ):
        if i not in seeds:       
            if len(seeds) <= 1:
                current_score = pairwise_score( peps[0], peps[i] )
                if current_score <= theta:
                    index_table[i] = np.copy(current_cluster)
                elif current_score > theta:
                    current_cluster += 1
                    index_table[i] = np.copy(current_cluster)
                    seeds.append(current_cluster)

            elif len(seeds) > 1:
                tmp_scores = pairwise_score_vectorized(list(peps[i]), np.array(seeds_peps))
                tmp_scores = np.array(tmp_scores)                
                current_best_idx = np.argmin(tmp_scores)
                current_best_score = tmp_scores[current_best_idx]
                if current_best_score <= theta:
                    index_table[i] = np.copy(seeds[current_best_idx])
                else:
                    current_cluster += 1
                    index_table[i] = np.copy(current_cluster)
                    seeds.append(current_cluster)
                    seeds_peps.append(list(peps[current_cluster]))

        if i%5000 == 0:
            print(i, "peptides processed")

    # move small clusters to sequence pool by setting to -1
    # re-number clusters after eliminating small clusters below threshold
    for i in list(set(index_table)):
        if np.sum(index_table == i) < min_size:
            index_table[index_table == i] = -1

    clusterID = 0
    for i in list(set(index_table)):
        if i > -1:
            current_idx = np.where(index_table == i)[0]
            index_table[current_idx] = clusterID
            clusterID += 1

    return index_table


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def write_aln(seqs, Fout, alnformat):
    """
    add sequences to record, write alignment file to disk
    'pseudo-alignment' as seqs have same length
    """

    records = []
    for idx, seq in enumerate(seqs):
        current_seq = SeqRecord( Seq(seq, IUPAC.protein), id="pep_"+str(idx), description="" )
        records.append(current_seq)
    aln = MultipleSeqAlignment(records)
    AlignIO.write(aln, Fout, alnformat)

    return aln




 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def peptide_hc(peplist, indextable, pepmatrix, substratematrix, RESULT_DIR, theta_hydro=0.5):
    """
    hierarchical clustering of input peptides
    peplist: list of peptides for clustering
    indextable: initial cluster assignment
    """

    def initialize_clustering(peplist, indextable, pepmatrix, write2disk=True):
        """
        prepares condensed peptide matrix
        compile dictionary of precomputed MSAs
        write alignments to tmp file to disk for use with hhsuit, hmmer, etc.
        """

        Nclust = len(set(list(indextable[indextable>=0])))

        pepmatrix_new = np.zeros(( np.shape(pepmatrix)[0], Nclust   ))
        clustermatrix_new = np.zeros(( Nclust )) * np.nan
        msadict = {}
        msaentropy = np.zeros(( Nclust ))
        msahydro   = np.zeros(( Nclust ))

        clusterID = 0
        for i in set(list(indextable)):
            if i > -1:              # omit sequence pool that has been set to -1
                idx_cluster = np.where(indextable == i)[0]
                pepmatrix_new[:,clusterID] = np.sum( pepmatrix[:,idx_cluster], axis=1)
                clustermatrix_new[clusterID] = clusterID
                current_msa = []
                records = []
                for j in idx_cluster:
                    current_msa.append( list(peplist[j]) )
                    records.append( SeqRecord( Seq(peplist[j], IUPAC.protein), id="pep_"+str(j), description="" ) )
                msadict[clusterID] =  np.array(current_msa)
                if write2disk:
                    aln = MultipleSeqAlignment(records)
                    AlignIO.write(aln, RESULT_DIR+"aln_"+str(clusterID)+".fa", "fasta")

                current_entro = compute_entropy_weblogo( np.array(current_msa) )
                msaentropy[clusterID] = current_entro
                current_hydro = compute_msa_hydrophobicity( np.array(current_msa) )
                msahydro[clusterID] = current_hydro 
            
                clusterID += 1

        return pepmatrix_new, clustermatrix_new, msadict, msaentropy, msahydro


    def write_msa2fasta(msa, id):
        """
        write msa in numpy array form to disk for use in post analysis
        msa: numpy array of aligned seqs
        id: identifier used for file name
        """
        records = []
        for i in range(np.shape(msa)[0]):
            current_seq = "".join( list(msa[i,:]))
            records.append( SeqRecord( Seq(current_seq, IUPAC.protein), id="pep_rel_"+str(i), description="" ) )
        aln = MultipleSeqAlignment(records)
        AlignIO.write(aln, RESULT_DIR+"aln_"+str(int(id))+".fa", "fasta")


    def update_indexed_array(arrayIN, newIDX, newVAL):
        """
        extend numpy array by value so that position corresponds to value ID
        """

        new_array = np.zeros(( int( np.max([np.max(arrayIN), newIDX]) + 1 )  ))
        new_array[0:len(arrayIN)] = np.copy(arrayIN)
        new_array[int(newIDX)] = newVAL

        return new_array


    def update_all(pepmatrix, clustermatrix, dmat, msa_dict, msa_entro, msa_hyd, dendro, cluster1, cluster2, prunelist):
        """
        update index table, cluster IDs as matrix indices
        intab: index table of cluster assignments
        cluster1, cluster2: clusters that are being merged
        """

        prunelist.append(int(cluster1))
        prunelist.append(int(cluster2))

        idx_cluster1 = np.where( clustermatrix == int(cluster1) )[0]
        idx_cluster2 = np.where( clustermatrix == int(cluster2) )[0]
        new_clusterID = int(np.nanmax(clustermatrix)) + 1

        # update msa dictionary
        new_msa = np.vstack( (msa_dict[cluster1], msa_dict[cluster2]) )
        msa_dict[new_clusterID] = new_msa
        write_msa2fasta(new_msa, new_clusterID)

        new_entropy = compute_entropy_weblogo(new_msa)
        msa_entro = update_indexed_array(msa_entro, new_clusterID, new_entropy)
        new_hydro = ( (np.shape(msa_dict[cluster1])[0] * msa_hyd[cluster1]) + (np.shape(msa_dict[cluster2])[0] * msa_hyd[cluster2]) ) / (np.shape(msa_dict[cluster1])[0] + np.shape(msa_dict[cluster2])[0])
        msa_hyd   = update_indexed_array(msa_hyd, new_clusterID, new_hydro)

        # update peptide matrix and index table
        pnr, pnc = np.shape(pepmatrix)
        newClust_counts = np.reshape( np.sum(pepmatrix[:,[int(idx_cluster1), int(idx_cluster2)]], 1), (pnr, 1) )
        pepmatrix = np.hstack( (pepmatrix, newClust_counts ) ) 
        clustermatrix = np.concatenate( (clustermatrix, [new_clusterID] ) )
        clustermatrix[int(idx_cluster1)] = np.nan
        clustermatrix[int(idx_cluster2)] = np.nan
        sel_cluster = np.isnan(clustermatrix)
        pepmatrix = pepmatrix[:, sel_cluster == False]        # not nan
        clustermatrix = clustermatrix[sel_cluster == False]
   
        # 3. update dendrogram
        dist = dmat[int(cluster1), int(cluster2)]
        size = np.shape(new_msa)[0]
        dendro.loc[len(dendro)] = (cluster1, cluster2, new_clusterID, dist, size)
    
        # update dmat, preserve dims so that index == clusterID
        nc,nr = dmat.shape
        new_dmat = np.zeros(( new_clusterID+1, new_clusterID+1 )) * np.nan  # copy to one added dimension
        new_dmat[0:nc,0:nc] = np.copy(dmat)

        pool = mp.Pool(processes=10)
        pool_set = np.where( [cluster not in prunelist for cluster in np.arange(new_clusterID+1)] )[0]
        pool_inputs = [(msa_dict, i, new_clusterID) for i in pool_set]
        new_dmat[pool_set, new_clusterID] = pool.map(average_linkage_parallel, pool_inputs)
        pool.close()

        # prune rows/columns of old clusters
        new_dmat[np.array(prunelist, dtype=int),:] = np.nan
        new_dmat[:,np.array(prunelist, dtype=int)] = np.nan
        new_dmat[new_clusterID, new_clusterID] = np.nan         # check that new diagonal is nan
  
        return pepmatrix, clustermatrix, new_dmat, dendro, msa_entro, msa_hyd, prunelist


    def evaluate_clusters(peptidematrix, clustermatrix, substratematrix, current_fts, bst_fts, bst, bst_mat, resultDF):
        """
        subroutine to evaluate clusters for discriminative power
        returns: best features, resultDF
        unpack/pack some vars for simpler code
        """
        bst_mat1, bst_mat3, bst_mat5 = bst_mat[0], bst_mat[1], bst_mat[2]
      
        if len(current_fts) >= 10:
            auc1   = roc_features(peptidematrix, substratematrix, current_fts[0:1])
            auc3   = roc_features(peptidematrix, substratematrix, current_fts[0:3])
            auc5   = roc_features(peptidematrix, substratematrix, current_fts[0:5])
            auc10  = roc_features(peptidematrix, substratematrix, current_fts[0:10])

            if auc1 > bst[0]:
                bst[0] = np.copy(auc1)
                bst_mat1 = np.copy( peptidematrix[:, current_fts[0:1]])
            if auc3 > bst[1]:
                bst[1] = np.copy(auc3)
                bst_mat3 = np.copy( peptidematrix[:, current_fts[0:3]])
                bst_fts = np.copy( clustermatrix[current_fts[0:3]] ) 
            if auc5 > bst[2]:
                bst[2] = np.copy(auc5)
                bst_mat5 = np.copy( peptidematrix[:, current_fts[0:5]])
            if auc10 > bst[3]:
                bst[3] = np.copy(auc10)
           
            bst_mat = (bst_mat1, bst_mat3, bst_mat5)

            resultDF.loc[len(resultDF)] = ( int(len(list(set(clustermatrix)))), int(clustermatrix[current_fts[0]]), int(clustermatrix[current_fts[1]]), int(clustermatrix[current_fts[2]]), int(clustermatrix[current_fts[3]]), int(clustermatrix[current_fts[4]]), auc1, auc3, auc5, auc10)

        return resultDF, bst_fts, bst, bst_mat



    ## INITIALIZE HIERARCHICAL CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print("Initializing hierarchical clustering ...")
    pepmat, clustmat, msa_dct, msa_entropy, msa_hydro = initialize_clustering(peplist, indextable, pepmatrix)
    np.savetxt(RESULT_DIR+"pepmat_clusters.txt", pepmat, fmt='%i')
    dendro = pd.DataFrame(columns=["cluster1", "cluster2", "clusterNew", "distance", "size"]) 
    dmat = compute_dmat(clustmat, msa_dct)     # initialize and compute initial dmat
    np.savetxt(RESULT_DIR+"dmat.txt", np.around(dmat,3), fmt='%.3f' )
    #dmat = np.loadtxt(RESULT_DIR+"dmat.txt")                           # more debugging


    ## INITIALIZED OUTPUT VARS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    predictDF = pd.DataFrame(columns=["N", "F1", "F2", "F3", "F4", "F5", "AUC1", "AUC3", "AUC5", "AUC10"])
    predict_hydDF = pd.DataFrame(columns=["N", "F1", "F2", "F3", "F4", "F5", "AUC1", "AUC3", "AUC5", "AUC10"])  
    best = np.zeros(( 4 ))
    best_hyd = np.zeros(( 4 ))
    best_mat = (None, None, None)
    best_mat_hyd = (None, None, None)
    best_fts = np.zeros(( 3 ))
    best_fts_hyd = np.zeros(( 3 ))
    list_pruned = []

    ## ITERATIVE CLUSTERING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print("Hierarchical clustering and feature selection ...")
    current_dist = 100                                            # smth high
    counter = 0
    while len(list(set(clustmat))) > 10: # and counter < 10:      # (FOR DEBUGGING)

        min_i, min_j = np.unravel_index(np.nanargmin(dmat, axis=None), dmat.shape) 
        current_dist = dmat[min_i, min_j]
        pepmat, clustmat, dmat, dendro, msa_entropy, msa_hydro, list_pruned = update_all(pepmat, clustmat, dmat, msa_dct, msa_entropy, msa_hydro, dendro, min_i, min_j, list_pruned)
        print(counter, len(list(set(list(clustmat)))), min_i, min_j,  current_dist, best[1])

        if np.sum(np.isnan(dmat)) == np.prod(np.shape(dmat)):
            break

        # feature selection and evalution
        fts = randomforest_features(pepmat, clustmat, substratematrix, msa_hydro, theta_hydro, hydrofilter=False)
        predictDF, best_fts, best, best_mat = evaluate_clusters(pepmat, clustmat, substratematrix, fts, best_fts, best, best_mat, predictDF)

        # only hydrophobic clusters
        fts_hyd = randomforest_features(pepmat, clustmat, substratematrix, msa_hydro, theta_hydro, hydrofilter=True)
        predict_hydDF, best_fts_hyd, best_hyd, best_mat_hyd = evaluate_clusters(pepmat, clustmat, substratematrix, fts_hyd, best_fts_hyd, best_hyd, best_mat_hyd, predict_hydDF)

        if counter % 100 == 0:  
            predictDF.to_csv(RESULT_DIR+"predict.txt", sep='\t', header=True, index=False)
            predict_hydDF.to_csv(RESULT_DIR+"predict_hyd.txt", sep='\t', header=True, index=False)
            dendro.to_csv(RESULT_DIR+"dendro.txt", sep='\t', header=True, index=False)
        counter += 1

    print(best_fts, "best1:", best[0], "best3:", best[1], "best5:", best[2])
    print(best_fts_hyd, "hyd_best1:", best_hyd[0], "hyd_best3:", best_hyd[1], "hyd_best5:", best_hyd[2])

    for i in best_fts:
        plot_weblogo(RESULT_DIR+'aln_'+str(int(i))+'.fa',   RESULT_DIR+'weblogo_aln_'+str(int(i))+'.png')
    for i in best_fts_hyd:
        plot_weblogo(RESULT_DIR+'aln_'+str(int(i))+'.fa',   RESULT_DIR+'weblogo_aln_'+str(int(i))+'.png')


    ## WRITE OUTPUT FILES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    predictDF.to_csv(RESULT_DIR+"predict.txt", sep='\t', header=True, index=False)
    predict_hydDF.to_csv(RESULT_DIR+"predict_hyd.txt", sep='\t', header=True, index=False)
    dendro.to_csv(RESULT_DIR+"dendro.txt", sep='\t', header=True, index=False)

    np.savetxt(RESULT_DIR+"roc1.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat[0]), 1), fmt='%i' )
    np.savetxt(RESULT_DIR+"roc3.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat[1]), 1), fmt='%i' )
    np.savetxt(RESULT_DIR+"roc5.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat[2]), 1), fmt='%i' )
    np.savetxt(RESULT_DIR+"roc1_hyd.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat_hyd[0]), 1), fmt='%i' )
    np.savetxt(RESULT_DIR+"roc3_hyd.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat_hyd[1]), 1), fmt='%i' )
    np.savetxt(RESULT_DIR+"roc5_hyd.txt", np.concatenate( (np.reshape( substratematrix,(len(substratematrix),1) ), best_mat_hyd[2]), 1), fmt='%i' )
  
    msa_entropy = pd.DataFrame({"clusterID":list(np.arange(len(msa_entropy))), "value": list(msa_entropy)})
    msa_entropy.to_csv(RESULT_DIR+"cluster_entropy.txt", sep='\t', header=True, index=False)
    msa_hydro = pd.DataFrame({"clusterID":list(np.arange(len(msa_hydro))), "value": list(msa_hydro)})
    msa_hydro.to_csv(RESULT_DIR+"cluster_hydrophobicity.txt", sep='\t', header=True, index=False)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def randomforest_features(ftsmat, indextable, classmat, msahydrophobicity, theta_hydrophobicity, hydrofilter=False):
    """
    random forest based feature selection
    only evaluate discriminative features that are enriched in class 1
    (rationale: sequence required for binding; if mixing, cannot simply sum feature counts in simple additive model)
    output:
        ftsmat: reduced dim feature matrix for ROC analysis
        features_idx: ranking of feature importance, indices corresponding to reduced dim feature matrix
    """

    class1 = classmat == 1
    class0 = classmat == 0
    counts_class1 = np.sum(ftsmat[class1,:], 0)
    counts_class0 = np.sum(ftsmat[class0,:], 0)
    potential_hits = counts_class1 > counts_class0

    if hydrofilter:
        current_hydro = msahydrophobicity[np.array(indextable, dtype=int)]
        hydro_sel = current_hydro > theta_hydrophobicity
        potential_hits = potential_hits * hydro_sel

    pothits_idx = np.where(potential_hits)[0] 
    new_ftsmat = ftsmat[:,potential_hits]
    new_indextable = indextable[potential_hits]

    if new_ftsmat.shape[1] > 0:            # sanity check that there are features left after filtering
        # Create a random forest classifier
        rfc = RandomForestClassifier(n_estimators=100, random_state = 0, n_jobs = -1)
        rfc.fit(new_ftsmat, classmat)

        feature_importances = pd.DataFrame(rfc.feature_importances_,
                                index=[n for n in range(np.shape(new_ftsmat)[1])],
                                columns=['importance']).sort_values('importance', ascending=False)
        features_idx = feature_importances.index
        features = np.array( pothits_idx[features_idx], dtype=int)

    else:
        features_idx = []
        features = []
        print("no features left ...")
    
    return features



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def roc_features(x, y, features):
    """
    get ROC AUC for select input features
    sum of feature counts as simple model of total number of putative binding sites
    """

    auc = roc_auc_score(y, np.sum(x[:,features], 1))
    
    return np.around(auc, 4)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def average_linkage(msa_dict, cluster1, cluster2):
    """
    compute average linkage of two multiple-sequence alignments
    msa_dict: dictionary with multiple-sequence alignments
    cluster1, cluster2: id's to retrieve alns from dictionary
    """
    
    if cluster1 > -1 and cluster2 > -1:
        msa1 = msa_dict.get( cluster1, np.nan )
        msa2 = msa_dict.get( cluster2, np.nan )
        nseqs1 = np.shape(msa1)[0]
        nseqs2 = np.shape(msa2)[0]
        scores = []
        distance = 0
        for i in range(nseqs1):
            for j in range(nseqs2):
                current_score = pairwise_score(msa1[i,:], msa2[j,:])
                scores.append(current_score)
        distance = -1 * np.mean( np.array(scores))
    else:
        distance = np.nan

    return distance


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def average_linkage_parallel( packaged_input ):
    """
    wrapper to run with multiprocessing library
    """
    
    msa_dict = packaged_input[0]
    cluster1 = packaged_input[1]
    cluster2 = packaged_input[2]

    aln1 = msa_dict[cluster1]
    aln2 = msa_dict[cluster2]
 
    if cluster1 < cluster2 :
        distance = pairwise_aln_vectorized(aln1, aln2)
    else:
        distance = np.nan

    return distance


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pairwise_aln_vectorized(aln1, aln2):
    """
    average linkage of alignments
    default is BLOSUM62 for now
    """

    n1, l1 = np.shape(aln1)
    n2, l2 = np.shape(aln2)
    L = np.min([l1, l2])
    scores = np.zeros(( n1, n2 ))
    
    for i in range( L ):
        scores += np.array( blosum62.loc[aln1[:,i]][aln2[:,i]] )
    avrg_score = - np.around(np.mean(scores),3) 

    return avrg_score


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_dmat(clustermatrix, msa_dict):
    """
    compute distance between clusters 
    """

    N = len(clustermatrix)

    pool = mp.Pool(processes=10)
    pool_inputs = [(msa_dict, i,j) for i in clustermatrix for j in clustermatrix]
    output = pool.map(average_linkage_parallel, pool_inputs)
    pool.close()

    dmat = np.reshape( np.array(output), (N,N) ) 
    print(dmat)

    return dmat


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_entropy_weblogo(MSA):
    """
    compute entropy of sequence alignment
    use weblogo implementation of relative entropy
    """

    aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    nrow, ncol = np.shape(MSA)
    counts  = np.zeros(( ncol, len(aminos) ))

    for ix, aa in enumerate(aminos):
        counts[:,ix] = np.sum(MSA == aa, 0)

    if np.sum(counts) > 0:        
        logodata = weblogo.LogoData.from_counts(weblogo.seq.protein_alphabet, counts )
        entropy = logodata.entropy
        S = np.sum(entropy)
    else:
        S = np.nan

    return S
      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_msa_hydrophobicity(MSA):
    """
    compute hydrophobicity of sequence alignment
    Kyte&Doolittle scale, normalized to unitary mean and sd (Pechmann et al. PNAS 2009)
    """

    hydroDF = pd.DataFrame(
        {"hyd": [4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.8, -0.9, -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5]},
        index = ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'E', 'Q', 'D', 'N', 'K', 'R']
        )

    N, L = np.shape(MSA)
    seqhyd = np.zeros(( N ))

    for i in range( N ):
        current_seq = MSA[i,:]
        current_hyd = np.array(hydroDF['hyd'][current_seq])
        seqhyd[i] = np.mean(current_hyd)
    
    msa_hyd = np.around(np.mean(seqhyd), 2)    

    return msa_hyd



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def save_intermediate_data(pepmatrix, substrmatrix, gclusters, peptidelist, RESULTDIR):
    """
    save intermediate data for quicker debugging
    """
    np.savetxt(RESULTDIR+"result_pepmat.txt", pepmatrix, fmt='%i')
    np.savetxt(RESULTDIR+"result_sbstrmat.txt", substrmatrix, fmt='%i')
    np.savetxt(RESULTDIR+"result_intab.txt", gclusters, fmt='%i')
    write_aln(peptidelist, RESULTDIR+"result_peplist.fa", "fasta")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_intermediate_data(RESULTDIR):
    """
    load intermediate data for quicker debugging
    """

    pepmat = np.loadtxt(RESULTDIR+"result_pepmat.txt")
    sbstrmat = np.loadtxt(RESULTDIR+"result_sbstrmat.txt")
    pc = np.loadtxt(RESULTDIR+"result_intab.txt")

    peplist = []
    tmp = SeqIO.parse(RESULTDIR+"result_peplist.fa", "fasta")
    for record in tmp:
        peplist.append(str(record.seq))

    return pepmat, sbstrmat, pc, peplist




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dendro_stats(RESULTDIR):
    """
    function to collect stats on clusters during hierarchical clustering
    size: average size of clusters 
    entropy: average entropy of clusters
    hydro: average hydrophobicity of clusters
    """

    dendrogram = pd.read_csv(RESULTDIR+"dendro.txt", header=0, sep='\t', index_col=False)
    predictions = pd.read_csv(RESULTDIR+"predict.txt", header=0, sep='\t', index_col=False)
    predictions_hyd = pd.read_csv(RESULTDIR+"predict_hyd.txt", header=0, sep='\t', index_col=False)
    indextable = np.loadtxt(RESULTDIR+"result_intab.txt")
    clent = pd.read_csv(RESULTDIR+"cluster_entropy.txt", header=0, sep='\t', index_col=False)
    clhyd = pd.read_csv(RESULTDIR+"cluster_hydrophobicity.txt", header=0, sep='\t', index_col=False)

    resultDF = pd.DataFrame(columns=['Iter', 'size_mean', 'size_sd', 'size_best3', 'size_best3_hyd', 'ent_mean', 'ent_sd', 'ent_best3', 'ent_best3_hyd', 'hyd_mean', 'hyd_sd', 'hyd_best3', 'hyd_best3_hyd'])

    maxgreedyClusterID = int( np.max(indextable) )
    maxClusterID = int( np.max(np.array(dendrogram['clusterNew'])) )

    sizearray = np.zeros(( maxClusterID + 1 )) * np.nan
    entarray = np.zeros(( maxClusterID + 1 )) * np.nan
    hydarray = np.zeros(( maxClusterID + 1 )) * np.nan

    # initial values at start of hierarchical clustering
    for i in range(maxgreedyClusterID):
        sizearray[i] = np.sum( indextable==i)
        entarray[i] = clent.iloc[i]['value'].item()
        hydarray[i] = clhyd.iloc[i]['value'].item()

    counter = 0 
    resultDF.loc[len(resultDF)] = (int(counter), np.around(np.nanmean(sizearray),2) , np.around(np.nanstd(sizearray),2), np.nan, np.nan,
                                                np.around(np.nanmean(entarray),2) , np.around(np.nanstd(entarray),2), np.nan, np.nan,
                                                np.around(np.nanmean(hydarray),2) , np.around(np.nanstd(hydarray),2), np.nan, np.nan  )

    # loop through dendrogram and update stats at each clustering iteration
    best3_size, best3_ent, best3_hyd = np.nan, np.nan, np.nan
    best3_size_hyd, best3_ent_hyd, best3_hyd_hyd = np.nan, np.nan, np.nan
    for i in range(len(dendrogram)):
        counter += 1

        current_cluster1 = int(dendrogram.iloc[i]['cluster1'])
        current_cluster2 = int(dendrogram.iloc[i]['cluster2'])
        current_clusterNew = int(dendrogram.iloc[i]['clusterNew'])
        current_size = int(dendrogram.iloc[i]['size'])

        sizearray[current_cluster1] = np.nan
        sizearray[current_cluster2] = np.nan
        sizearray[current_clusterNew] = current_size

        entarray[current_cluster1] = np.nan
        entarray[current_cluster2] = np.nan
        entarray[current_clusterNew] = float( clent.iloc[current_clusterNew]['value'].item() )

        hydarray[current_cluster1] = np.nan
        hydarray[current_cluster2] = np.nan
        hydarray[current_clusterNew] = float( clhyd.iloc[current_clusterNew]['value'].item() )

        N = len( sizearray[~np.isnan(sizearray)] ) + 1          # discount -1 'clusters' 
        current_pred = predictions[predictions['N']==N]
        if len(current_pred) > 0:
            current_f1 = int( np.array(current_pred['F1'])[0] )
            current_f2 = int( np.array(current_pred['F2'])[0] )
            current_f3 = int( np.array(current_pred['F3'])[0] )
            best3_size = np.around(np.mean(np.array([sizearray[current_f1], sizearray[current_f2], sizearray[current_f3] ] )), 2 )
            best3_ent = np.around(np.mean(np.array([entarray[current_f1], entarray[current_f2], entarray[current_f3] ] )), 2 )
            best3_hyd = np.around(np.mean(np.array([hydarray[current_f1], hydarray[current_f2], hydarray[current_f3] ] )), 2 )

        current_pred_hyd = predictions_hyd[predictions_hyd['N']==N]
        if len(current_pred_hyd) > 0:
            current_f1_hyd = int( np.array(current_pred_hyd['F1'])[0] )
            current_f2_hyd = int( np.array(current_pred_hyd['F2'])[0] )
            current_f3_hyd = int( np.array(current_pred_hyd['F3'])[0] )
            best3_size_hyd = np.around(np.mean(np.array([sizearray[current_f1_hyd], sizearray[current_f2_hyd], sizearray[current_f3_hyd] ] )), 2 )
            best3_ent_hyd = np.around(np.mean(np.array([entarray[current_f1_hyd], entarray[current_f2_hyd], entarray[current_f3_hyd] ] )), 2 )
            best3_hyd_hyd = np.around(np.mean(np.array([hydarray[current_f1_hyd], hydarray[current_f2_hyd], hydarray[current_f3_hyd] ] )), 2 )

        resultDF.loc[len(resultDF)] = (int(counter), np.around(np.nanmean(sizearray),2) , np.around(np.nanstd(sizearray),2), best3_size, best3_size_hyd,
                                                np.around(np.nanmean(entarray),2) , np.around(np.nanstd(entarray),2), best3_ent, best3_ent_hyd,
                                                np.around(np.nanmean(hydarray),2) , np.around(np.nanstd(hydarray),2), best3_hyd, best3_hyd_hyd  )

    resultDF.to_csv(RESULTDIR+"clustering_summary.txt", header=True, index=False, sep='\t', na_rep="NA")

    return resultDF



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_weblogo(file_aln, file_out):
    """
    plot sequence logo with Weblogo CLI
    """

    weblogo = '/home/sebastian/sw/anaconda3/bin/weblogo'
    formatting_args = " -F png_print -A protein -s large -c monochrome --annotate '1,2,3,4,5,6,7' --errorbars NO --fineprint '' "
    cmd_weblogo = weblogo + ' -f ' + file_aln + ' -o ' + file_out + formatting_args

    subprocess.call(cmd_weblogo, shell=True)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hcluster_motif(fastaIN, RESULTDIR, windowsize=7, theta_hydro=0.5):
    """
    wrapper function for discriminative motif finding through hierarchical clustering
    fastaIN: input fasta sequence, class has to be in seq record description
    RESULTDIR: directory for output files
    """

    fasta = SeqIO.parse(fastaIN, "fasta")

    if not os.path.exists(RESULTDIR):
        os.mkdir(RESULTDIR)
    shutil.copy(fastaIN, RESULTDIR)

    pepmat, peplist, sbstrmat = peptide_matrix(fasta, windowsize)
    greedyclusters = greedy_clustering(peplist)
    save_intermediate_data(pepmat, sbstrmat, greedyclusters, peplist, RESULTDIR) 
    print(np.max(greedyclusters), len(list(set(greedyclusters))), len(greedyclusters), pepmat.shape )

    peptide_hc(peplist, greedyclusters, pepmat, sbstrmat, RESULTDIR, theta_hydro) 
    summary = dendro_stats(RESULTDIR)



if __name__ == "__main__":


    ## LOAD AUX DATA
    blosum62 = pd.read_csv("../data/yeast/blosum62.txt", header=0, index_col=0, sep='\t')

    ## RUN 
    hcluster_motif('../data/processed/ssb_structsubstr.fasta', '../data/processed/result_structures_hyd05/', windowsize=7, theta_hydro=0.5)                          

    hcluster_motif('../data/processed/ssb_peaks_R1.fasta', '../data/processed/result_peaks_hyd05_R1/', windowsize=7, theta_hydro=0.5)       
    hcluster_motif('../data/processed/ssb_peaks_R2.fasta', '../data/processed/result_peaks_hyd05_R2/', windowsize=7, theta_hydro=0.5)      
    hcluster_motif('../data/processed/ssb_peaks_R3.fasta', '../data/processed/result_peaks_hyd05_R3/', windowsize=7, theta_hydro=0.5)      

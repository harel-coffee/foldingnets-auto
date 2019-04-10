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
 
 
yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))
orflist = pd.read_csv("../data/processed/orflist.txt", header=0, sep='\t', index_col=None)


blosum62 = pd.read_csv("../data/yeast/blosum62.txt", header=0, index_col=0, sep='\t')

aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def yeast_bgfreqs(seqs):
    """
    extracts yeast AA background frequencies
    """

    aafreq = np.zeros(( len(aminos) ))

    for seq in seqs.values():
        current_seq = np.array(seq)

        for ax, aa in enumerate(aminos):
            aafreq[ax] += sum(current_seq == aa)


    aafreq /= np.sum(aafreq)

    return aafreq

#f = yeast_bgfreqs(yeast_seqs)

#print(f)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def generate_peps(seq):
    """
    chopes sequence into peptides with sliding window
    """
 
    w = 7
    L = len(seq)
 
    seq_bag = []
    sequence = str(seq.seq)
 
    for i in range(L-w):
        pep = sequence[i:i+w]
        seq_bag.append(pep) 
 
    return seq_bag
 
#s = generate_peps(yeast_seqs['YAL001C'])
 




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_nrseqs(ORFlist):
    """
    extracts all peptide sequences from ORFs in ORFlist
    filters for non-redundant sequences
    """

    peptides = []

    for i in ORFlist:
        current_seq = yeast_seqs[i]
        current_peps = generate_peps(current_seq)
        peptides += current_peps

    peptides = sorted( set(peptides) )

    return peptides



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pairwise_score(pep1, pep2):
    """
    computes pairwise alignment score
    assumes peptides of same length
    default is BLOSUM62 for now
    """

    score = 0

    for i in range( min( len(pep1), len(pep2) ) ):
        score += blosum62[pep1[i]][pep2[i]]

    return score







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def greedy_clustering(peps):
    """
    performs greedy initial clustering of peptide seqs
    with threshold and BLOSUM62 similarity
    """

    min_size = 10

    index_table = np.zeros(( len(peps) ))
    theta = int( round( len(peps[0]) * 1.7, 0) ) #default of Hammock paper
    seeds = [0]
    current_cluster = 0

    # check all seeds, join highest that is also above theta. else new seed
    for i in range(1, len(peps) ):

        if i not in seeds:
            
            if len(seeds) <= 1:
                current_score = pairwise_score(peps[seeds[0]], peps[i] )

                if current_score >= theta:
                    index_table[i] = current_cluster

                elif current_score < theta:
                    current_cluster += 1
                    index_table[i] = current_cluster
                    seeds.append(i)

            elif len(seeds) > 1:
                tmp_scores = np.zeros(( len(seeds) ))
                for jx, j in enumerate(seeds):
                    tmp_scores[jx] = pairwise_score( peps[j], peps[i] )
                current_best = np.argmax(tmp_scores)

                if tmp_scores[current_best] >= theta:
                    index_table[i] = current_best
                else:
                    seeds.append(i)
                    current_cluster += 1
                    index_table[i] = current_cluster

        if i%1000 == 0:
            print(i, "processed")


    # move small clusters to sequence pool
    # set to -1 to indicate sequence pool
    for i in list(set(index_table)):
        if np.sum(index_table == i) < min_size:
            index_table[index_table == i] = -1


    return index_table




# introduce a step that only considers clusters of min size, and
# puts the rest of the seqs into the sequence pool





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def write_aln(seqs, Fout, alnformat):
    """
    add sequences to record
    write alignment file to disk
    'pseudo-alignment' as seqs have same length
    """

    records = []
    for idx, seq in enumerate(seqs):
        current_seq = SeqRecord( Seq(seq, IUPAC.protein), id="pep_"+str(idx) )
        records.append(current_seq)

    aln = MultipleSeqAlignment(records)

    # alnformat: clustal, fasta
    AlignIO.write(aln, Fout, alnformat)

    return aln







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hmm_align(msa1, msa2):
    """
    use HHsuite3 to compute similarity score of two MSAs via hhalign
    msa1, msa2 are files on disk
    """

    hhalign = '/home/user/sw/hhsuite/bin/hhalign'

    cmd_hhalign = hhalign + ' -i ' + msa1 + ' -t ' + msa2 
    output = subprocess.getoutput(cmd_hhalign).split('\n')


    result = []
    parse = False
    counter = 0
    for line in output:
        #print(line)
        if "No Hit" in line:
        	colnames = line.split()
        	result = pd.DataFrame(columns=colnames[0:8])
        	parse = True
        elif parse == True:
            scores = line.split()[0:8]
            if len(scores) == 8:
                result.loc[counter] = scores
                counter += 1


    if len(result) > 0:
        score = np.max(np.array( result['Score'] ) )
    else:
        score = 0

    #print(score)
    #print(result)

    return score


#hmm_align("tmp.pos.fa", "tmp.neg.fa")







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def pairwise_similarity(pep1, intab1, pep2, intab2):
    """
    compute pairwise similarity between clusters in 2 sets
    pep1/pep2: input sequences
    intab1/intab2: index table of initial clustering
    """

    min_size = 10

    clust_pos = []
    clust_neg = []

    for i in list(set(intab1)):
        if np.sum(intab1 == i) >= min_size and i != np.nan:
            clust_pos.append(i)

    for i in list(set(intab2)):
        if np.sum(intab2 == i) >= min_size and i != np.nan:
            clust_neg.append(i)

    scoremat = np.zeros(( len(intab1), len(intab2) )) * np.nan

    for ix,i in enumerate( clust_pos ):
        id_seqs_pos_cl1 = np.where(intab1 == i)[0]
        seqs_pos_cl1 = []
        for i in id_seqs_pos_cl1:
            seqs_pos_cl1.append(pep1[i])

        aln_pos = write_aln(seqs_pos_cl1, "tmp.pos.fa", "fasta")

        for jx, j in enumerate( clust_neg ):
            id_seqs_neg_cl1 = np.where(intab2 == j)[0]
            seqs_neg_cl1 = []
            for i in id_seqs_neg_cl1:
                seqs_neg_cl1.append(pep2[i])

            aln_neg = write_aln(seqs_neg_cl1, "tmp.neg.fa", "fasta")

            score = hmm_align("tmp.pos.fa", "tmp.neg.fa")
            #print(i, j, score)

            scoremat[ix,jx] = score

    #print(scoremat)
    #print( np.max(scoremat) )

    return scoremat


# NEXT

# - cluster above min size: HMM
# - search sequence pool with HMMs

# - cluster similarity with HMM-HMM alignment (HH-suite)

# - pos and neg clusters of highest similarity are removed (no discrim seqs)
# - cluster merging within sets


# assym: use clusters from pos set: keep if no simil to any cluster in neg set
# crit for scanning: pos cluster has min size and min coverage! 
# simil: either HMM-HMM dist, or HMM-scan and count of occs




def complete_linkage_score(cluster1, cluster2):
    """
    compute complete linkage similarity score
    """
    print("bla")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def exterminate_cluster(indextable, clusterid):
    """
    remove clusters from both input sets
    - set indextable to np.nan to indicate sequence is no longer considered
    """
    
    indextable2 = np.copy(indextable)

    indextable2[ indextable2 == clusterid ] = np.nan


    return indextable2





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def merge_clusters(indextable, clust_id1, clust_id2):
    """
    merge highly similar clusters in same input set
    set cluster id to that of cluster1 (arbitrary)
    """

    # re-assign labels
    indextable[indextable == clust_id2] = clust_id1

    return indextable





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def expand_cluster(indextable, clusterid, seqid):
    """
    expand cluster from sequence pool
    """
    
    indextable[indextable == seqid] = clusterid


    return indextable





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def seqs_from_df(orfdf):
    """
    get substrates from orflist DF
    use "get_nrseqs" fct for nr peptides
    """

    sbstr_pos = orflist[orfdf['substrate']=="+"]['orf']
    sbstr_neg = orflist[orfdf['substrate']=="-"]['orf']

    peptides_pos = get_nrseqs(sbstr_pos)
    peptides_neg = get_nrseqs(sbstr_neg)


    return peptides_pos, peptides_neg

# TO FINISH

# MERGE CLUSTER FUNCTION
# full pairwise simil within, merge top clusters within


# EXPAND CLUSTER FUNCTION
# full loop through clusteres + seqs, expand

# PRUNE CLUSTER FUNCTION
# full pairwise simil between, prune







if __name__ == "__main__":


    # generate peptide sequences
    pep_pos, pep_neg = seqs_from_df(orflist)

    pep_pos = pep_pos[0:1000]
    pep_neg = pep_neg[0:1000]

    # perform initial greedy clustering
    intab_pos = greedy_clustering(pep_pos)
    intab_neg = greedy_clustering(pep_neg)


    sm = pairwise_similarity(pep_pos, intab_pos, pep_neg, intab_neg)

    #print(sm[5:13,5:13])
    print( np.nansum(sm) )
    
    max_x, max_y = np.unravel_index(np.nanargmax(sm, axis=None), sm.shape)

    print(max_x, max_y)

    print( list(set(intab_pos)) )


    intab_pos2 = exterminate_cluster(intab_pos, 0)
    intab_neg2 = exterminate_cluster(intab_neg, 1)



    print( np.sum(intab_pos == intab_pos2) / len(intab_pos) )
    print( np.sum(intab_neg == intab_neg2) / len(intab_neg) )

    

    """
    sm1 = pairwise_similarity(pep_pos, intab_pos, pep_neg, intab_neg)

    #print(sm1[5:13,5:13])
    print( np.nansum(sm1) )
    
    max_ind = np.unravel_index(np.nanargmax(sm1, axis=None), sm1.shape)

    print(max_ind)

    print(sm1[ max_ind ])


    """

    #xy = np.unravel_index(np.argmin(simimat[simimat!=np.nan], axis=None), np.shape(simimat[simimat!=np.nan]) ) 

    #xy = max(simidict.iterkeys(), key=(lambda key: simidict[key]))
    #xy 	= max(simidict, key=lambda key: simidict[key])
    #xy2 = min(simidict, key=lambda key: simidict[key])

    #print(xy, xy2)
    #print(simidict[xy], simidict[xy2])

    

    #hmmdict_pos, intab_pos = exterminate_cluster(hmmdict_pos, intab_pos, xy2[0])
    #hmmdict_neg, intab_neg = exterminate_cluster(hmmdict_neg, intab_neg, xy2[1])

    
 

    #print(simimat[ xy[0], xy[1] ])

    
  
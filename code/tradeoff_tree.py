import sys,os
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import re
import subprocess

from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_motif(aln, seq):
    """
    test a motif/peptide in a sequence, return position of occurence using regex
    aln: seq alignment of possible motifs/peptides to test
    seq: seq to search
    """

    positions = []
    for rec in aln:
        current_motif = str(rec.seq)
        m = re.search(current_motif, seq[:])
        if m:
            m_start, m_end = m.span()
            positions.append(m_start)

    return positions


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def peptide_counts(structureDF, yeastseqs, aln_list):
    """
    compile dataframe of peptide counts given list of peptide alignments
    """

    resultDF = pd.DataFrame(columns=['ORF', 'substrate', 'counts'])

    for i in range(len(structureDF)):
        current_orf = structureDF.iloc[i]['ORF']
        current_seq = str( yeastseqs[current_orf].seq )
        current_sbstr = structureDF.iloc[i]['substrate']
        if current_sbstr == "+":
            current_sbstr = 1
        elif current_sbstr == "-":
            current_sbstr = 0

        current_counts = 0
        for j in range(len(aln_list)):
            current_positions = test_motif(aln_list[j], current_seq)
            current_counts += len(current_positions)

        resultDF.loc[len(resultDF)] = (current_orf, current_sbstr, current_counts)

    return resultDF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def crossref_data(contdensDF, pepcountDF, noptDF, fOUT):
    """
    compile DF with contact density, peptide counts, and nopt hd cluster counts
    """

    resultDF = pd.DataFrame(columns=['ORF', 'substrate', 'cd', 'counts', 'nopt'])

    for i in range(len(contdensDF)):
        current_orf = contdensDF.iloc[i]['ORF']
        current_sbstr = contdensDF.iloc[i]['substrate']
        current_CD = np.around( contdensDF.iloc[i]['normCD'], 3 )
        if current_sbstr == "+":
            current_sbstr = 1
        elif current_sbstr == "-":
            current_sbstr = 0

        current_pepcount = pepcountDF[pepcountDF['ORF']==current_orf]['counts'].item()
        current_noptCC   = noptDF[noptDF['ORF']==current_orf]['CC'].item()

        resultDF.loc[len(resultDF)] = (current_orf, current_sbstr, current_CD, current_pepcount, current_noptCC)

    resultDF.to_csv(fOUT, header=True, index=False, sep='\t')

    return resultDF


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def randomforest_prediction(inputDF):
    """
    leave-one-out cross-validation of predictive power of features
    """

    def cross_validation(substrate, features):
        """
        leave-one-out cross-validation
        """

        result = np.zeros(( len(substrate) )) * np.nan
        for i in range(len(substrate)):
            selection = np.ones(( len(substrate) ), dtype=bool)
            selection[i] = False

            substrate_train = np.copy( substrate[selection] )
            substrate_test  = np.copy( substrate[~selection] )

            features_train = np.copy( features[selection,:] )
            features_test = np.copy( features[~selection,:] )

            rfc = RandomForestClassifier(n_estimators=100, random_state = 0, n_jobs = -1)
            rfc.fit(features_train, substrate_train)
            y_predict = rfc.predict(features_test)
            acc = accuracy_score(substrate_test, y_predict)
            result[i] = acc

        accuracy = np.around( np.mean(result),3 )
        
        return accuracy


    labels   = np.array( inputDF['substrate'], dtype=int )

    features_all = np.array( inputDF[['cd', 'counts', 'nopt']] )
    features_cd_nopt = np.array( inputDF[['cd', 'nopt']] )
    features_cd = np.array( inputDF[['cd']] )
    features_nopt = np.array( inputDF[['nopt']] )
    features_counts = np.array( inputDF[['counts']] )

    score_all = cross_validation(labels, features_all)
    score_cd_nopt = cross_validation(labels, features_cd_nopt)
    score_cd = cross_validation(labels, features_cd)
    score_nopt = cross_validation(labels, features_nopt)
    score_counts = cross_validation(labels, features_counts)

    output = [score_all, score_cd_nopt, score_cd, score_nopt, score_counts]

    return output



if __name__ == "__main__":


    # load data
    structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
    yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))
    contdens = pd.read_csv("../data/processed/all_contact_density_normalized.txt", header=0, sep='\t', index_col=None)
    nopt = pd.read_csv("../data/processed/nopt_clusters.txt", sep='\t', header=0, index_col=None)

    # best features based on highest AUC3 score in '../data/processed/result_structures_hyd05/predict.txt'
    aln1 = AlignIO.read("../data/processed/result_structures_hyd05/aln_204.fa", "fasta")
    aln2 = AlignIO.read("../data/processed/result_structures_hyd05/aln_691.fa", "fasta")
    aln3 = AlignIO.read("../data/processed/result_structures_hyd05/aln_742.fa", "fasta")

    # best hydrophobic features based on highest AUC3 score in '../data/processed/result_structures_hyd05/predict_hyd.txt'
    aln_hyd1 = AlignIO.read("../data/processed/result_structures_hyd05/aln_1734.fa", "fasta")
    aln_hyd2 = AlignIO.read("../data/processed/result_structures_hyd05/aln_2064.fa", "fasta")
    aln_hyd3 = AlignIO.read("../data/processed/result_structures_hyd05/aln_1987.fa", "fasta")

    pepcounts = peptide_counts( structures, yeast_seqs, [aln1, aln2, aln3] )
    pepcounts_hyd = peptide_counts( structures, yeast_seqs, [aln_hyd1, aln_hyd2, aln_hyd3] )


    # run stuff
    data_all = crossref_data(contdens, pepcounts, nopt, "../data/processed/data_tradeoffs.txt")
    data_all_hyd = crossref_data(contdens, pepcounts_hyd, nopt, "../data/processed/data_tradeoffs_hyd.txt")

    res_all = randomforest_prediction(data_all)
    res_hyd = randomforest_prediction(data_all_hyd)


    resDF = pd.DataFrame({"feature":['all', 'cd+nopt', 'cd', 'nopt', 'counts'], 'all':list(res_all), 'hyd':list(res_hyd)})
    resDF.to_csv("../data/processed/pred_tradeoffs.txt", header=True, index=False, sep='\t')

    print(resDF)
import sys, os
import numpy as np
import pandas as pd
import glob
import shutil
import subprocess
import pickle
import glob

from Bio import AlignIO
from Bio import Align
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa



aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_hydrophobicity = {"I": 4.5,  "V":  4.2, "L":  3.8, "F":  2.8, "C":  2.5, "M":  1.9, "A":  1.8, "G": -0.4, "T": -0.7, "S": -0.8,
                     "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5 }




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def write_domain_file(family):
    """
    write STAMP domain file and seed file
    """

    fname = 'db.' + family + '.domains'
    fname2 = 'tmp.seed.domain'
    current_fam = structures[structures['fam']==family]

    file = open(fname, 'w')
    file2 = open(fname2, 'w')

    for ix, i in enumerate(current_fam['ORF']):
        pdb = '../data/pdb/minimized/' + i + '.pdb'
        outline = pdb + ' ' + str(i) + ' ' + "{ ALL } " + '\n'
        file.write(outline)
        if ix < 1:
            file2.write(outline)
    file.close() 
    file2.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def run_stamp(family):
    """
    wrapper function to generate a structure-based multiple-sequence alignment with STAMP
    """
    import subprocess
    from subprocess import check_output


    # get rough initial alignment with scan
    seed = 'tmp.seed.domain'
    db = 'db.' + family + '.domains'
    prefix ="".join(family.split("."))
    scan = prefix + '.scan'
    sortd = prefix + '.sorted'
 
    output1 = subprocess.run([r"stamp", "-l", seed, "-n", "2", "-s", "-slide", "5", "-d", db, "-prefix", prefix])
    
    output2 = check_output([r"sorttrans", "-f", scan, "-s", "Sc", "2.5"], encoding='UTF-8')
    fileOUT = open(sortd, 'w')
    fileOUT.write(output2)
    fileOUT.close()
        
    prefix2 = 'tmp.stamp' 
    output3 = check_output([r"stamp", "-l", sortd, "-prefix", prefix2], encoding='UTF-8')
   
    fileOUT2 = open("stamp.log", 'w')
    fileOUT2.write(output3)
    fileOUT2.close()

    structaln = pd.DataFrame(columns=['Orf1', 'Orf2', 'Sc', 'RMS', 'Len1', 'Len2', 'Aling', 'Nfit', 'SeqID'])
    df_counter = 0

    for line in output3.split('\n'):
        if "Pair" in line:
            current_line = line.split()
            ORF1 = current_line[2][0:7]
            ORF2 = current_line[3][0:7]
            Sc   = current_line[4]
            RMS  = current_line[5]
            L1   = current_line[6]
            L2   = current_line[7]
            ALen = current_line[8]
            Nfit = current_line[9]
            Sid  = current_line[12]

            structaln.loc[df_counter] = [ORF1, ORF2, Sc, RMS, L1, L2, ALen, Nfit, Sid]
            df_counter += 1

    structaln.to_csv("../data/processed/structaln.txt", sep='\t', header=True, index=False)


    aln_files = glob.glob('tmp.stamp.*')
    current_max = 0
    for entry in aln_files:
        current_entry = entry.split('.')[2]
        if 'mat' not in entry:
            if int(current_entry) > int(current_max):
                current_max = current_entry
    aln_file = 'tmp.stamp.' + str(current_max)
     
    cmd =  "aconvert -in b -out c < " + aln_file + " > " + prefix + ".aln"
    output4 = subprocess.call(cmd, shell=True)

    tmp = glob.glob("tmp.*")
    if len(tmp) > 0:
        for f in tmp:
            os.remove(f)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def convert_aln(alnfile):
    """
    simplified alignment for plotting in R
    """

    alignment = AlignIO.read(open(alnfile), "clustal")

    list_id = []
    dict_id = {}
    list_substr = []
    list_ix = []   

    for ix, record in enumerate(alignment):
        current_record = record.id
        current_orf = current_record.split("_")[0]
        if current_orf not in list_id and current_orf != "space":
            current_record = current_record.split('_')[0]
            list_id.append(current_record)
            list_ix.append(ix)
            dict_id[ix] = current_record

            substr = structures[structures['ORF'] == current_record]['substrate'].item()
            if substr == '+':
                list_substr.append(1)
            elif substr == '-':
                list_substr.append(0)

    alignment = np.array(alignment)
    alignment = alignment[np.array(list_ix),:]

    aln_df = pd.DataFrame(alignment)
    aln_df.insert(loc=0, column='idx', value=list(list_id) )
    aln_df['substrate'] = list_substr
    aln_df = aln_df.sort_values('substrate', ascending=True)
    aln_df.set_index('idx', inplace=True, drop=True) 
    aln_df.drop('substrate', 1, inplace=True)

    aln_df.to_csv("../data/processed/c3718_aln.txt", header=False, index=True, sep='\t')

    return list_id







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_alignment(fname):
    """
    read in multiple-sequence alignment
    """

    alignment = AlignIO.read(open(fname), "clustal")
    #print("Alignment length %i" % alignment.get_alignment_length())
    for record in alignment :
        seq = np.array(record.seq)
        seq = seq[seq!='-']
        seq = list(seq)
        seq = "".join(seq)

        ident = record.id.split('_')[0]
        print(seq, ident, len(seq), len(yeast_seqs[ident]))

#read_alignment("test.aln")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def project_hyd(alnfile):
    """
    compute and project hydrophobicity profiles onto alignment
    """

    alignment = AlignIO.read(open(alnfile), "clustal")

    nrow = len(alignment)
    ncol = alignment.get_alignment_length()
    alnmat = np.zeros(( nrow, ncol ))
    list_substr = []
    list_id = []
    list_ix = []

    for ix, record in enumerate(alignment):
        seq_aln = np.array(record.seq)
        seq_ref = "".join(list(seq_aln[seq_aln!='-']))
        current_record = record.id.split('_')[0]

        if "dssp" not in record.id and current_record not in list_id and current_record != "space":
            list_id.append(current_record)
            list_ix.append(ix)

            substr = structures[structures['ORF'] == current_record]['substrate'].item()
            if substr == '+':
                list_substr.append(1)
            elif substr == '-':
                list_substr.append(0)

            # get profile
            hyd = np.zeros(( len(seq_ref) ))
    
            for ax, aa in enumerate(seq_ref):
                hyd[ax] = aa_hydrophobicity[aa]

            hyd_aln = np.zeros(( len(seq_aln) )) * np.nan

            pos_aln = 0
            pos_ref = 0

            while (pos_aln < ncol):
                if seq_aln[pos_aln] == '-':
                    pos_aln += 1
                else:
                    hyd_aln[pos_aln] = hyd[pos_ref]
                    pos_aln += 1
                    pos_ref += 1

            alnmat[ix,:] = hyd_aln

    alnmat = alnmat[np.array(list_ix),:]

    alndf = pd.DataFrame(alnmat)
    alndf.insert(loc=0, column='substrate', value=list(list_substr) )
    alndf.sort_values('substrate', ascending=True, inplace=True)
    #alndf.drop('substrate', 1, inplace=True)
    alndf.to_csv("../data/processed/c3718_aln_hyd.txt", sep='\t', header=True, index=False)

    return alndf, alnmat




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def agg_scores(D, inplist):
    """
    global aggregation score from Tango profile
    """
    AP = pd.DataFrame(columns=['ORF','agg','substr'])
    counter = 0  

    for i in inplist:
        score = np.sum( agg_dict[i] )
        substr = structures[structures['ORF'] == i]['substrate'].item() 
        AP.loc[counter] = [i, score, substr]
        counter += 1

    AP.set_index('ORF', inplace=True, drop=True) 
    AP.to_csv("../data/processed/c3718_aggregation.txt", sep='\t', header=True, index=True)

    return AP




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def project_asa(alnfile, D, fileOUT):
    """
    compute and project profiles from dict onto alignment
    """

    alignment = AlignIO.read(open(alnfile), "clustal")

    nrow = len(alignment)
    ncol = alignment.get_alignment_length()
    alnmat = np.zeros(( nrow, ncol ))
    list_substr = []
    list_id = []
    list_ix = []

    for ix, record in enumerate(alignment):
        seq_aln = np.array(record.seq)
        seq_ref = "".join(list(seq_aln[seq_aln!='-']))
        current_record = record.id.split('_')[0]


        if "dssp" not in record.id and current_record not in list_id and current_record != "space":
            list_id.append(current_record)
            list_ix.append(ix)

            substr = structures[structures['ORF'] == current_record]['substrate'].item()
            if substr == '+':
                list_substr.append(1)
            elif substr == '-':
                list_substr.append(0)


            # get profile from dictionary D

            ident = record.id.split('_')[0]
            asa = D[ident]

            asa_aln = np.zeros(( len(seq_aln) )) * np.nan

            pos_aln = 0
            pos_ref = 0

            while (pos_aln < ncol):
                if seq_aln[pos_aln] == '-':
                    pos_aln += 1
                else:
                    asa_aln[pos_aln] = asa[pos_ref]
                    pos_aln += 1
                    pos_ref += 1

            alnmat[ix,:] = asa_aln
 
    alnmat = alnmat[np.array(list_ix), :]

    alndf = pd.DataFrame(alnmat)
    alndf.insert(loc=0, column='substrate', value=list(list_substr) )
    alndf.sort_values('substrate', ascending=True, inplace=True)
    #alndf.drop('substrate', 1, inplace=True)
    alndf.to_csv(fileOUT, sep='\t', header=True, index=False)


    return alndf, alnmat



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def save_df(df, sel, fileOUT):
    """
    convert to DF for easy handling of missing values / mixed types
    """
    seldf = df[sel]
    sbstr = list ( df['substrate'] )
    seldf.insert(loc=0, column='substrate', value=sbstr)
    seldf.to_csv(fileOUT, sep='\t', header=False, index=False)






if __name__ == '__main__': 


    structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
    yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))


    write_domain_file('c.37.1.8')
    run_stamp('c.37.1.8')
    if os.path.exists("c3718.aln"):
        list_c3718 = convert_aln("c3718.aln")

        # write ASA and RSA into dicts
        asa_c3718 = {}
        rsa_c3718 = {}
        for i in list_c3718:
            current_dssp = load_dssp("../data/pdb/dssp/" + i + ".dssp")
            asa_c3718[i] = np.array(current_dssp['ASA'])
            rsa_c3718[i] = np.array(current_dssp['RSA'])


        hyd_df, hydmat = project_hyd('c3718.aln')

        # load agg preds
        with open('../data/processed/allseqs.agg.pkl', 'rb') as f:
            agg_dict = pickle.load(f)

        aggscores = agg_scores(agg_dict, list_c3718)


        asa_df, asamat = project_asa('c3718.aln', asa_c3718, "../data/processed/c3718_aln_asa.txt")
        asamat = pd.DataFrame(asamat)
        asamat.to_csv("../data/processed/c3718.asa", sep='\t', header=False, index=False)
        asa_mean = np.nanmean(asamat, 0)

        agg_df, aggmat = project_asa('c3718.aln', agg_dict, "../data/processed/c3718_aln_agg.txt")

        buried  = np.where(asa_mean < 50)[0]
        exposed = np.where(asa_mean >= 50)[0]
        aligned = np.where( sum(np.isnan(hydmat),0) < 15 )[0]
        notaln  = np.where( sum(np.isnan(hydmat),0) >= 15)[0]

        sel_buried_aligned  = np.intersect1d(buried, aligned)
        sel_buried_notaln   = np.intersect1d(buried, notaln)
        sel_exposed_aligned = np.intersect1d(exposed, aligned)
        sel_exposed_notaln  = np.intersect1d(exposed, notaln)

        save_df(hyd_df, sel_buried_aligned, "../data/processed/fig2_hyd_buried_alinged.txt")
        save_df(hyd_df, sel_buried_notaln, "../data/processed/fig2_hyd_buried_notalinged.txt")
        save_df(hyd_df, sel_exposed_aligned, "../data/processed/fig2_hyd_exposed_alinged.txt")
        save_df(hyd_df, sel_exposed_notaln, "../data/processed/fig2_hyd_exposed_notalinged.txt")

        save_df(agg_df, sel_buried_aligned, "../data/processed/fig2_agg_buried_alinged.txt")
        save_df(agg_df, sel_buried_notaln, "../data/processed/fig2_agg_buried_notalinged.txt")
        save_df(agg_df, sel_exposed_aligned, "../data/processed/fig2_agg_exposed_alinged.txt")
        save_df(agg_df, sel_exposed_notaln, "../data/processed/fig2_agg_exposed_notalinged.txt")

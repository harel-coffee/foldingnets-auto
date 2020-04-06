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
import multiprocessing as mp
import subprocess

 

aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_substrate_seqs(listpos, listneg, yeast_fasta, prefixOUT):
    """
    parse sequences of orfdf SSB substrates  
    """

    F = open(prefixOUT+'.fasta', 'w')
    Fpos = open(prefixOUT+'_pos.fasta', 'w')
    Fneg = open(prefixOUT+'_neg.fasta', 'w')

    for orf in listpos:
        if orf in yeast_fasta.keys():
            current_record = yeast_fasta[orf]
            current_seq = current_record.seq
            F.write( ">" + orf + " " + str(1) + '\n' )
            F.write( str(current_seq) + '\n')
            Fpos.write( ">" + orf + '\n' )
            Fpos.write( str(current_seq) + '\n')
    Fpos.close()

    for orf in listneg:
        if orf in yeast_fasta.keys():
            current_record = yeast_fasta[orf]
            current_seq = current_record.seq
            F.write( ">" + orf + " " + str(0) + '\n' )
            F.write( str(current_seq) + '\n')
            Fneg.write( ">" + orf + '\n' )
            Fneg.write( str(current_seq) + '\n')

    Fneg.close()
    F.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_peak_seqs(peaks, listpos, listneg, yeast_fasta, prefixOUT, cntrl="neg"):
    """
    parse sequences of SSB substrates and SSB binding regions
    based on mRNA site translated while SSB binding
    binding site inferred through offset of allowing for min/max ribosome tunnel length + 20aa window
    CONTROL SEQS HAVE RANDOM ELEMENT to it to balance number of pos/neg peptides
    cntrl:  "pos" to get FASTA including non binding regions from same protein seqs
            "neg" to get FASTA with about same number of peps from non-binding protein seqs
    """

    w = 20                      # window to search for motif
    theta_length = 35

    F = open(prefixOUT+'.fasta', 'w')
  
    # putative binding regions
    Fpos = open(prefixOUT+'_pos.fasta', 'w')
    for i in range(len(peaks)):
        current_orf = peaks.iloc[i]['ORF']
        if current_orf in listpos:
            current_start = peaks.iloc[i]['Start'] - 50 - w # maximal tunnel length (published is ~50)
            current_end   = peaks.iloc[i]['End'] - 25       # minimal tunnel length (published is 23, but unlikely that chais is fully, maximally extended)
            if current_orf in yeast_fasta.keys():
                current_seq = yeast_fasta[current_orf].seq
                if current_start > 1:
                    current_peak = current_seq[current_start:current_end]
                    F.write( ">" + current_orf + "_" + str(i+1) + " " + "1" + '\n' )
                    F.write( str(current_peak) + '\n')
                    Fpos.write( ">" + current_orf + "_" + str(i+1) + '\n' )
                    Fpos.write( str(current_peak) + '\n')
    Fpos.close()
    
    # non-binding regions of positive sequences as control, random fragments at least 100 away from realy binding site
    if cntrl == "pos":
        Fneg = open(prefixOUT+'_neg.fasta', 'w')
        for orf in list(set(peaks['ORF'])):
            if orf in listpos:
                current_motifs = peaks[peaks['ORF']==orf]
                current_seq = str(yeast_fasta[orf].seq)

                for i in range(len(current_motifs)):
                    current_start = int(current_motifs.iloc[i]['Start'])
                    current_end = int(current_motifs.iloc[i]['End'])
                    current_fragment = current_seq[100:(current_start-100)]
                    if len(current_fragment) > theta_length:
                        current_offset = int( np.random.choice( len(current_fragment)-theta_length, 1) )
                        F.write( ">" + orf + "_" + str(i+1) + " " + "0" + '\n' )
                        F.write( str(current_fragment[current_offset:(current_offset+theta_length)]) + '\n')
                        Fneg.write( ">" + orf + "_" + str(i+1) + '\n' )
                        Fneg.write( str(current_fragment[current_offset:(current_offset+theta_length)]) + '\n')
                    current_seq = current_seq[current_end:len(current_seq)]

                if len(current_seq) > theta_length:
                    current_offset = int( np.random.choice( len(current_seq)-theta_length, 1) )
                    F.write( ">" + orf + "_" + str(i+2) + " " + "0" + '\n' )
                    F.write( str(current_seq[current_offset:(current_offset+theta_length)]) + '\n')
                    Fneg.write( ">" + orf + "_" + str(i+2) + '\n' )
                    Fneg.write( str(current_seq[current_offset:(current_offset+theta_length)]) + '\n')
        Fneg.close()

    # negative ssb seqs, but only limited to cyto/nuclear. get same amount of random fragments from seqs
    if cntrl == "neg":
        Fneg = open(prefixOUT+'_neg.fasta', 'w')
        for orf in list(listneg):
            if orf in yeast_fasta.keys():
                current_description = str( loqate[loqate['ORF']==orf]['Localization'] )
                if 'cyto' in current_description or 'nucl' in current_description or 'below' in current_description:
                    current_seq = str(yeast_fasta[orf].seq)
                    if len(current_seq) > 100:
                        current_offset = int( np.random.choice( len(current_seq)-100, 1) )
                        F.write( ">" + orf + " " + "0" + '\n' )
                        F.write( str(current_seq[current_offset:(current_offset+30)]) + '\n')
                        Fneg.write( ">" + orf + '\n' )
                        Fneg.write( str(current_seq[current_offset:(current_offset+30)]) + '\n')
        Fneg.close()

    if cntrl == "full":
        Fneg = open(prefixOUT+'_neg.fasta', 'w')
        for orf in list(listneg):
            if orf in yeast_fasta.keys():
                current_description = str( loqate[loqate['ORF']==orf]['Localization'] )
                if 'cyto' in current_description or 'nucl' in current_description or 'below' in current_description:
                    current_seq = str(yeast_fasta[orf].seq)
                    if len(current_seq) > 100:
                        F.write( ">" + orf + " " + "0" + '\n' )
                        F.write( str(current_seq) + '\n')
                        Fneg.write( ">" + orf + '\n' )
                        Fneg.write( str(current_seq) + '\n')
        Fneg.close()

    F.close()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def run_meme(FposIN, FnegIN, DIRout):
    """
    wrapper to run conventional MEME motif search
    requires local install of MEME motif finding software
    """

    meme = '/home/sebastian/sw/meme-5.1.0/bin/meme'
    flags = '-objfun de -minw 7 -maxw 9 -protein -nmotifs 5 -minsites 10 -maxiter 500 -nostatus'
    # -np (run in parallel once cpu's are freed upgre)

    if os.path.exists(DIRout):
        print("overwriting existing results directory!")

    cmd_meme = meme + ' ' + FposIN + ' -neg ' + FnegIN + ' -oc ' + DIRout + ' ' + flags
    subprocess.call(cmd_meme, shell=True)







if __name__ == "__main__":


    ## LOAD data
    yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))
    structures = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
    loqate = pd.read_csv("../data/yeast/yeast_loqate.txt", header=0, sep='\t', index_col=None)
    ssb_strong = pd.read_csv("../data/yeast/ssb_consensus.txt", header=0, index_col=None)
    ssb_not = pd.read_csv("../data/yeast/ssb_cons_not.txt", header=0, index_col=None)
    ssb_peaks = pd.read_csv("../data/yeast/ssb_peaks.txt", header=0, index_col=None, sep='\t')



    ## PARSE SEQUENCES INTO FASTA FILES
    # 3 'runs' or 'replicas' as negative control sequences are random
    parse_peak_seqs(ssb_peaks, list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_peaks_R1', cntrl="neg")
    parse_peak_seqs(ssb_peaks, list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_peaks_R2', cntrl="neg")
    parse_peak_seqs(ssb_peaks, list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_peaks_R3', cntrl="neg")

    parse_peak_seqs(ssb_peaks, list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_peaks_poscntrl', cntrl="pos")
    parse_peak_seqs(ssb_peaks, list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_peaks_negfull', cntrl="full")
    parse_peak_seqs(ssb_peaks, list(structures[structures['substrate']=='+']['ORF']), list(structures[structures['substrate']=='-']['ORF']), yeast_seqs, '../data/processed/ssb_structures_peaks_negfull', cntrl="full")

    parse_substrate_seqs(list(ssb_strong['ssb']), list(ssb_not['nssb']), yeast_seqs, '../data/processed/ssb_substrates')
    parse_substrate_seqs(list(structures[structures['substrate']=='+']['ORF']), list(structures[structures['substrate']=='-']['ORF']), yeast_seqs, '../data/processed/ssb_structsubstr')


    ## RUN ANALYSES
    run_meme('../data/processed/ssb_substrates_pos.fasta', '../data/processed/ssb_substrates_neg.fasta', '../data/processed/meme_ssb_substrates')
    run_meme('../data/processed/ssb_structsubstr_pos.fasta', '../data/processed/ssb_structsubstr_neg.fasta', '../data/processed/meme_ssb_structsubstr')
    run_meme('../data/processed/ssb_peaks_negfull_pos.fasta', '../data/processed/ssb_peaks_negfull_neg.fasta', '../data/processed/meme_ssb_peaks_full')
    run_meme('../data/processed/ssb_peaks_pos.fasta', '../data/processed/ssb_peaks_neg.fasta', '../data/processed/meme_ssb_peaks')
    run_meme('../data/processed/ssb_peaks_poscntrl_pos.fasta', '../data/processed/ssb_peaks_poscntrl_neg.fasta', '../data/processed/meme_ssb_peaks_poscntrl')

import sys, os, io
import numpy as np
import pandas as pd
from Bio import SeqIO
import re
import random




#----------------------------------------------------------------------------------------------------
def read_genelist(fname):
    """
    reads list of substrate ORF names into list
    """
    ORFs = [line.rstrip('\n') for line in open(fname)]

    return ORFs 


#----------------------------------------------------------------------------------------------------
def load_matchid(fname):
    """
    match ORF name to geneID
    """  
    matchid = {}

    fileIN = io.open(fname, 'r')
    line = fileIN.readline()
    while line:
        current_line = line.split()
        ID = int(current_line[0])
        orf = str(current_line[1])
        matchid[ID] = orf
        line = fileIN.readline()

    return matchid



# 1. LOAD SUBSTRATES  AND DATA

# yeast genome
yeast = read_genelist("../data/yeast/yeast_genome.txt")

# yeast proteome fasta
yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))



# substrate lists (SSB1: Doering et al. Cell 2017; SSB2: Willmund et al. Cell 2013)
ssb1_strong = read_genelist("../data/substrates/ssb1_substrates_strong.txt")
ssb1_all = read_genelist("../data/substrates/ssb1_substrates_all.txt")
ssb2_strong = read_genelist("../data/substrates/ssb2_substrates_strong.txt")
ssb2_all = read_genelist("../data/substrates/ssb2_substrates_all.txt")


# all genes of proteins that interaction somehow with SSB1/2
#ssb_all = list(set( ssb1_all + ssb2_all) )

# in both datasets
ssb_weak = list( set(ssb1_all) & set(ssb2_all) )
ssb_all = list( set( ssb_weak + ssb1_strong + ssb2_strong) ) 

# all genes of proteins that do not interact with SSB1/2
ssb_none = list( set(yeast) - set(ssb_all) )
ssb_n = list( set(ssb_none) )

# genes of proteins that strongly interact with SSB1/2 in either of two experiments 
ssb_p = list( set(ssb1_strong + ssb2_strong) )

ssb_strong = pd.DataFrame({"ssb":ssb_p})
ssb_strong.to_csv("../data/yeast/ssb_consensus.txt", header=True, index=False)
ssb_not = pd.DataFrame({"nssb":ssb_n})
ssb_not.to_csv("../data/yeast/ssb_cons_not.txt", header=True, index=False)


# identifiers
uniprot_annotation = pd.read_csv("../data/yeast/yeastSwiss-ProtAC.txt", header=0, sep='\t', index_col=None)


# structural classification
domains = pd.read_csv("../data/structural_classification/domainTableDump.txt", header=0, sep='\t', index_col=None)
sfa = pd.read_csv("../data/structural_classification/dir.cla.scope.2.06-stable.txt", skiprows=4, sep='\t', 
	index_col=None, header=None, names  = ["domain_code", "pdb_id", "chain", "sfa", "ID", "descr"] )
sfa['scopID'] = sfa.domain_code.str[:6]

regions = pd.read_csv("../data/structural_classification/domainRegionTableDump.txt", header=0, sep='\t', index_col=None)
matchid = load_matchid("../data/structural_classification/match.id")

swissprot = pd.read_csv("../data/structural_classification/swissprot_yeast.index", skiprows=2, header=4, index_col=None, sep='\t')


#1. filter for single domain
domain_counts = domains['geneID'].value_counts()
single_domain = list((domain_counts[domain_counts == 1]).index )
domains_single = domains[ domains['geneID'].isin(single_domain)]  


#2. filter for reliable structure prediction method
#hiconf_method = ['pfam', 'pdbblast', 'orfeus']
hiconf_method = ['pdbblast', 'orfeus', 'pcons']
domains_single_hiconf = domains_single[ domains_single['method'].isin(hiconf_method)]


scop_sd = pd.DataFrame(columns=['geneID', 'ORF', 'domainID', 'pdb_code', 'scop', 'chain', 'start', 'end'])

i = 0
for geneid in domains_single_hiconf['geneID']:
    
    orf = matchid.get(geneid, 'none')
    
    current_entry = domains_single_hiconf[domains_single_hiconf['geneID']==geneid]

    domainID = int( current_entry['domainID'].values[0] )
    pdbID = str( current_entry['pdb_id'].values[0] ).lower()

    current_region = regions[regions['domainID']==domainID]
    current_domain = current_region[ ['segment', 'start', 'end'] ].values[0]
    current_chain = current_domain[0]
    current_start = current_domain[1]
    current_end = current_domain[2]

    current_pdb = pdbID.replace('_', '')
    if len(current_pdb) == 4:
        current_structure = 'd' + current_pdb + current_chain.lower()
    elif len(current_pdb) == 5:
        current_structure = 'd' + current_pdb
    
    current_sfa =  sfa[sfa['scopID']== current_structure ]
    
    scop = current_sfa['sfa'].values
   
    if len( list(set(scop) ) ) == 1 and orf != 'none':
        scop_sd.loc[i] = [geneid, orf, domainID, current_structure[1:], scop[0], current_chain, current_start, current_end]
        i += 1


scop_sd.to_csv("../data/processed/scop_singledomain.txt", sep='\t', header=True, index=False)



#3. get sfa/fam assignments
#----------------------------------------------------------------------------------------------------
def substrates_by_family(scop_sd, ssb_p, ssb_n):
    """
    map SSB substrates onto protein families
    """
    unique_family = list(set( scop_sd['scop'] ))
    output = pd.DataFrame(columns=['family', 'total', 'SSB+', 'SSB-'])

    counter = 0
    for i in unique_family:
        current_scop = scop_sd[scop_sd['scop'] == i]
        current_orfs = current_scop['ORF'].values

        if len(current_orfs) > 1:

            ssb_int = 0
            ssb_nint = 0
            for orf in current_orfs:
                if orf in ssb_p:
                    ssb_int += 1
                if orf in ssb_n:
                    ssb_nint += 1

            if ssb_int > 0 and ssb_nint > 0 and ssb_int+ssb_nint > 3:
                output.loc[counter] = [i, ssb_int+ssb_nint, ssb_int, ssb_nint]
                counter += 1

    return output


print( substrates_by_family(scop_sd, ssb_p, ssb_n) )

#----------------------------------------------------------------------------------------------------
def substrates_by_superfamily(scop_sd, ssb_p, ssb_n, unip):
    """
    map SSB substrates onto protein superfamilies
    """
    unique_family = list(set( scop_sd['scop'] ))
    unique_sf = list(set( ['.'.join( scopid.split('.')[0:3]) for scopid in unique_family ] ))

    sf_counts = pd.DataFrame(columns=['sf', 'total', 'SSBp', 'SSBn'])

    orflist = pd.DataFrame(columns=['orf', 'uniprotID', 'substrate', 'sf', 'family', 'Length'] ) #, 'fam', 'pdb'])

    counter_sf_counts = 0
    counter_orflist = 0

    for i in unique_sf:
        current_sf = i + '.'
        idx = scop_sd.apply(lambda x: current_sf in x['scop'], axis=1)
        current_orfs = scop_sd['ORF'][idx].values
        
    
        if len(current_orfs) > 1:
            ssb_int = 0
            ssb_nint = 0

            ssb_substr = []
            ssb_nsubstr = []

            for orf in list(set(current_orfs)):
                uniprot_values = unip[unip['OLN'] == orf].values
                current_scopsd = scop_sd[scop_sd['ORF']==orf]
                current_family = current_scopsd['scop'].values  
                current_family = str(current_family[0])
                
                if orf in list(yeast_seqs.keys()):
                    L = len( yeast_seqs[orf].seq )
                else:
                    L = 0


                if len(uniprot_values) > 0:
                    uniprotID = uniprot_values[0][1]
                else:
                    uniprotID = orf

                if orf in ssb_p:
                    ssb_int += 1
                    ssb_substr.append(orf)
                    
                    orflist.loc[counter_orflist] = [orf, uniprotID, "+", i, current_family, L]
                    counter_orflist += 1

                if orf in ssb_n:
                    ssb_nint += 1
                    ssb_nsubstr.append(orf)

                    orflist.loc[counter_orflist] = [orf, uniprotID, "-", i, current_family, L]
                    counter_orflist += 1

            if ssb_int > 1 and ssb_nint > 1:
                sf_counts.loc[counter_sf_counts] = [i, ssb_int+ssb_nint, ssb_int, ssb_nint] #, ssb_substr, ssb_nsubstr
                counter_sf_counts += 1

    topsf = list(sf_counts['sf'].values)
    orflist = orflist[orflist['sf'].isin(topsf)]
    orflist.to_csv("../data/processed/orflist.txt", sep='\t', header=True, index=False)

    return sf_counts, orflist


ssb_sfam, orflist = substrates_by_superfamily(scop_sd, ssb_p, ssb_n, uniprot_annotation)
ssb_sfam.to_csv("../data/processed/ssb_sfam.txt", sep='\t', header=True, index=False)



#----------------------------------------------------------------------------------------------------
def check_orflist_swissmodel(orflist, anno):
    """
    check if and how well covered proteins are in swiss-prot/expasy
    """
    
    for i in orflist['orf']:
        uniprot = anno[anno['OLN'] == i]['uniprotID'].values
        
        current_model = swissprot[swissprot['UniProtKB_ac'].isin(uniprot) ]



        #print(i, current_model.shape[0])
        print(current_model[ ['UniProtKB_ac', 'from', 'to', 'template']])


#check_orflist_swissmodel(orflist, uniprot)


#----------------------------------------------------------------------------------------------------
def get_sequences(orflist, seqs):
    """
    pull out sequences from yeast fasta file
    """

    records = []
    for i in orflist:
        #print(i)
        #if i in yeast_seqs:
        records.append( seqs[i] )

        SeqIO.write(records, "../data/processed/allseqs.fa", "fasta")

    return records



#ssb_seqs = get_sequences(orflist['orf'], yeast_seqs)



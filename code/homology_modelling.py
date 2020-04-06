import sys, os
import numpy as np
import pandas as pd
import glob
import shutil
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from Bio import SeqIO


# read yeast protein sequences
yeast_seqs = SeqIO.to_dict(SeqIO.parse("../data/yeast/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa", "fasta"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_pir(SEQ, seqout):
    """
    parse fasta sequence in pir format
    """

    identifier = str(SEQ.id)
    sequence   = str(SEQ.seq)

    if sequence[-1] != '*':
        sequence = sequence + '*'

    with open(seqout, 'w') as PIR:
        PIR.write('>P1;'+identifier+'\n')
        PIR.write('sequence:'+identifier+':::::::0.00: 0.00'+'\n')
        PIR.write(sequence+'\n')

    #return identifier

#parse_pir(yeast_seqs['YML078W'], 'tmp.pir')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def model_protein(SEQ, tmplts):
    """
    align query sequence to sequence of template PDB structure
    """
    best_model = haspdb
    best_score = 0
    best_rmsd = 0
    best_seqid = 100

    for template in tmplts:
        PDB = template[0:5]

        identifier = SEQ + '_' + PDB
        seqfile = SEQ + '.pir'

        chain  = PDB[-1]
        coords = PDB + '.pdb'

        model_dir = '../data/homology/' + SEQ + '_' + PDB
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
        shutil.copy2('../data/pdb/templates/'+coords, model_dir) 
        os.chdir(model_dir)
        parse_pir(yeast_seqs[SEQ], seqfile)

        env = environ()
        aln = alignment(env)
        mdl = model(env, file=PDB, model_segment=('FIRST:'+chain,'LAST:'+chain))
    
        aln.append_model(mdl, align_codes=PDB, atom_files=coords)
        aln.append(file=seqfile, align_codes=SEQ)
        aln.align2d()
    
        aln.write(identifier + '.ali', alignment_format='PIR')
        aln.write(identifier + '.pap', alignment_format='PAP')

        env = environ()
        a = automodel(env, identifier + '.ali',
                  knowns=PDB, sequence=SEQ,
                  assess_methods=(assess.DOPE,                      
                              assess.GA341))
        a.starting_model = 1
        a.ending_model = 5
        a.make()
            
        homology_models = glob.glob(SEQ+'*.pdb')

        for homology_model in homology_models:
            env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
            env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

            mdl = complete_pdb(env, homology_model)
            s = selection(mdl)  
                
            score = s.assess_dope()      
    
            rmsd, seqid = compute_rmds(homology_model, template)

            if score < best_score:
                best_score = np.copy(score)
                best_rmsd = np.copy(rmsd)
                best_seqid = np.copy(seqid)
                best_model = model_dir + '/' + homology_model

        os.chdir('../../../code/')

    return best_model, best_score, best_rmsd, best_seqid



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def refine_model(pdb, seqid):
    """
    refine and complete exisiting PDB
    """
    pdbcode = pdb.split('.')[0]
    coords = pdb + '.pdb'

    best_model = 'none'
    best_score = 0

    model_dir = '../data/homology/' + seqid + '_' + pdbcode
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
        shutil.copy2('../data/pdb/templates/'+coords, model_dir) 
    os.chdir(model_dir)

    #write reference sequence to file
    parse_pir(yeast_seqs[seqid], seqid+'.seq')

    # Get the sequence of PDB file, and write to an alignment file
    e = environ()
    m = model(e, file=pdbcode)
    aln = alignment(e)
    aln.append_model(m, align_codes=pdbcode)
    aln.append(file=seqid+'.seq', align_codes=seqid)
    aln.salign()
    aln.write(file=seqid+'.ali')

    # check that there are actually gaps in PDB to model
    g = []
    for record in SeqIO.parse(seqid+'.ali', "pir"):
        g.append( record.seq.count('-') )
   
    # if PDB complete, no need to model
    if np.max(np.array(g)) < 1:
        best_model = model_dir + '/' + coords
        best_score = -100000
       
    else:
        # Load the automodel class
        log.verbose()
        env = environ()

        a = automodel(env, seqid+'.ali',
                  knowns=pdbcode, sequence=seqid,
                  assess_methods=(assess.DOPE,                      
                              assess.GA341))
        a.starting_model = 1
        a.ending_model = 3
        a.make()
            
        refined_models = glob.glob(seqid+'*.pdb')

        for refined_model in refined_models:
            env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
            env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters
  
            mdl = complete_pdb(env, refined_model)
            s = selection(mdl)  
            score = s.assess_dope()    

            if score < best_score:
                best_score = np.copy(score)
                best_model = model_dir + '/' + refined_model  
    
    os.chdir('../../../code/')

    return best_model, best_score


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_rmds(modl, tmpl):
    """
    RMSD between template and model
    based on MODELLER example for align3d function
    """
    env = environ()
    mdl = model(env)
    aln = alignment(env)
    for code in [modl,  tmpl]:
        mdl.read(file=code)
        aln.append_model(mdl, align_codes=code, atom_files=code)
    aln.align(gap_penalties_1d=(-600, -400))
    aln.align3d(gap_penalties_3d=(0, 2.0))
    #aln.write(file='tmp.rmsd.ali')

    # superpose the two structures 
    mdl = model(env, file=modl)
    atmsel = selection(mdl).only_atom_types('CA')
    mdl2 = model(env, file=tmpl)
    r = atmsel.superpose(mdl2, aln)

    rmsd = r.drms
    seqid = mdl.seq_id

    return rmsd, seqid





if __name__ == '__main__': 

    pdblist = pd.read_csv("../data/processed/pdblist.txt", header=0, sep=' ', index_col=None)
    templatelist = pd.read_csv("../data/processed/templates.txt", header=0, sep='\t', index_col=None)
    output = pd.DataFrame(columns=['orf', 'model', 'score', 'rmsd', 'seqid'])
    counter = 0
    
    for i in pdblist['orf']:
        best_score = 'none'
        best_rmsd = 'none'
        best_seqid = 'none'

        current_entry = pdblist[pdblist['orf']==i]
        sf = str(current_entry['sf'].values[0])
        haspdb = str(current_entry['pdb'].values[0])

        if haspdb == 'none':
            current_sf =  current_entry['sf']
            current_sf = current_sf.item()
            current_templates = list( templatelist[templatelist['sf'] == current_sf ]['template'] )
            if len(current_templates) > 0:
                best_model, best_score, best_rmsd, best_seqid = model_protein(i, current_templates)
                
            output.loc[counter] = [i, best_model, best_score, best_rmsd, best_seqid]
            if best_score != 'none':
                shutil.copy2(best_model, '../data/pdb/models')

        else:            
            current_pdb = str(current_entry['pdb'].item() )
            best_model, best_score = refine_model(current_pdb, i)
            output.loc[counter] = [i, best_model, best_score, best_rmsd, best_seqid] 
            if best_score != 'none':
                shutil.copy2(best_model, '../data/pdb/models')

        counter += 1

    output.to_csv("../data/processed/best_models.txt", sep='\t', header=True, index=False)




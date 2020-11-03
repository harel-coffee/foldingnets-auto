# Programmed trade-offs in protein folding networks
Structure (2020) DOI: 10.1016/j.str.2020.09.009

### Sebastian Pechmann
Web: pechmannlab.net\
Contact: sebastian@pechmannlab.net



### Analyses:
#### 1. Structural classification of chaperone substrates
structural_classification.py

#### 2. Homology modeling
homology_modelling.py

#### 3. Minimization of PDB structures
preparePDB.py

#### 4. Compute structural phylogeny with 'superpose'
structphylo.py

#### 5. Compute structure-based multiple-sequence alignment with 'Stamp'
structaln.py

#### 6. Contact maps of hydrophobic contacts with 'CSU'
csu.py

#### 7. Analysis of protein contact networks
protein_network.py

#### 8. Analysis of binding motifs and discriminative peptide sequences
binding_meme.pyi\
binding_motifs.py

#### 9. Analysis of codon usage
rarecodons.py

#### 10. Integrated analysis
tradeoff_tree.py


### Dependencies: 
R: ggplot2, reshape2, pROC, viridis, cowplot, , MASS, rpart, ape, ggExtra, matrixStats, ggpubr

Python: numpy, pandas, networkx, subprocess, Biopython, itertools, sklearn, multiprocessing, random, glob, MODELLER, re, shutil, Weblogo 

External: superpose, STAMP, DSSP, CSU, MEME


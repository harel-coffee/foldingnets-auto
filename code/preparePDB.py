from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys, os
import argparse
import warnings
import logging

import numpy as numpy
import pandas as pd

import warnings
import pdbfixer 
import logging
from simtk.openmm.app import PDBFile



def fixPDB(pdb, pdbname):
    """
    prepares the PDB structure for simulation/minimization usingn the openMM PDBfixer
    """
    
    add_hyds = True

    fixer = pdbfixer.PDBFixer(filename=pdb)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater='keep_crystallographic_water')
    
    if add_hyds:
        fixer.addMissingHydrogens(7.0) # only if we want protons!
    
    outfile = open(pdbname, 'w')
    PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
    outfile.close()  







def minimizePDB(pdb, pdbname):
    """
    minimize PDB structure with OpenMM standard minnimization protocol
    """

    pdb = app.PDBFile(pdb)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    # use app.Modeller to add hydrogens and solvent
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometers)

    # app.PDBFile.writeFile(modeller.topology, modeller.positions, open('1o9s_modeller_tip3p.pdb', 'w'))

    # prepare system and integraton
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometers, constraints=HBonds, rigidWater=True,
        ewaldErrorTolerance=0.0005)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator.setConstraintTolerance(0.00001)

    # prepare simulation
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'CpuThreads': '1'}
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)

    # minimize
    logging.info('Minimizing...')
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(pdbname, 1))
    simulation.step(1)



if __name__ == '__main__': 

    pdbmat = pd.read_csv("../data/processed/pdball.txt", header=0, sep='\t', index_col=None)
    pdbmat = pdbmat[pdbmat['type']=='model']

    counter = 0
    for orf in pdbmat['ORF'][29:]:
        type = pdbmat[pdbmat['ORF'] == orf]['type'].item()
        pdbOUT = '../data/pdb/minimized/' + orf + '.pdb'

        if type == 'pdb':
            #current_pdb = pdbmat[pdbmat['ORF'] == orf]['pdb'].item()
            #pdbIN = '../data/pdb/templates/' + current_pdb + '.pdb'
            pdbIN = pdbmat[pdbmat['ORF'] == orf]['pdb'].item()

        else:
            pdbIN = pdbmat[pdbmat['ORF'] == orf]['pdb'].item()


        fixPDB(pdbIN, pdbOUT)
        minimizePDB(pdbOUT, pdbOUT)
        #os.remove('tmp.pdb')
        print(counter, orf)
        counter += 1

        


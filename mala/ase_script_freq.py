# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 16:46:18 2020

@author: sacam5
"""
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.calculators.vasp import Vasp2
from ase.io import read, write, Trajectory

#calc = get_frequencies(method="Frederiksen", direction="central")
#Local definition of atom...atom distance measurement
#
def atom_dist(atoms, i1, i2):
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])

    size = np.sqrt(np.vdot(vec,vec))
    vec = vec/size
    
    return size

molecule = read('POSCAR', calculator=EMT())
BFGS(molecule).run(fmax=0.01)
molecule.set_calculator(calc)

calc.set(xc='PBE',
         encut=400,
         gamma=True,
         prec='Accurate',
         lreal=False,
         ispin=2,)


vib = Vibrations(molecule)
vib.run()
vib.summary()
vib.write_mode(-1) 

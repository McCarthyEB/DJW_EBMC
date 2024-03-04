# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:55:21 2020

@author: dave
"""
import numpy as np
import sys
import csv
import os
#
from ase import Atom
from ase import Atoms

from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength

from ase.io import read, write, Trajectory, castep

import ase.calculators.castep as castep_calc
#
# Local definition of atom...atom distance measurement
#
def atom_dist(atoms, i1, i2):
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])

    size = np.sqrt(np.vdot(vec,vec))
    vec = vec/size
    
    return size

#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
car_file = sys.argv[1]                                           # Directory where we will do the calculations
splstr=car_file.split('.')
stem=splstr[0]
castep_cellfile="%s.cell" % stem
#
print(f"------------------------------------------------------------")
print(f"car file to convert          : %s" % car_file )
print(f"Stem                         : %s" % stem     )
print(f"CASTEP cell file will be     : %s" % castep_cellfile )
print(f"------------------------------------------------------------")
#
atoms=read(car_file)

write(castep_cellfile, atoms ,format='castep-cell')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))

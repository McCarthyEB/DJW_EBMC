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
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.visualize import view

from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength
from ase.optimize import QuasiNewton

from ase.io import read, write, Trajectory
from ase.constraints import FixAtoms

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
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")

# check fix_z:
print("Will Freeze all atoms below z =", float(sys.argv[1]))
fix_z = float(sys.argv[1])

atoms=read('POSCAR')
write('check_start.cif', atoms ,format='cif')

fix_list=[]
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))

    if (atoms.positions[iatom][2] < fix_z):
        print("fix this one")
        fix_list.append(iatom)

print("Fix list: ", fix_list)


cons = FixAtoms(fix_list)
atoms.set_constraint(cons)

write("POSCAR_fixed", atoms, vasp5=True)
    
view(atoms)      


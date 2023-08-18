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

def atom_dist(atoms, i1, i2):
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])

    size = np.sqrt(np.vdot(vec,vec))
    vec = vec/size
    
    return size

atoms=read('POSCAR')
write('check_start.cif', atoms ,format='cif')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
    
#
# Define the bond that you will change
# The bond should move from 1 to 2 so that index_2 is for the atom that
# will be moved.
index_1 = 4
index_2 = 10
#
# Step size in Angstroms
step = 0.1
# Use the EMT calculator for the forces and energies:
atoms.set_calculator(EMT())

# constrain this bond
c = FixBondLength( index_1, index_2 )
atoms.set_constraint(c)

dist=[]
energy=[]
csv_headers=['Step', 'distance / Angs', 'Energy / eV']
# step through geometries
for istep in range(0,3):
#   atoms.positions[index_2] = atoms.positions[index_1]+(1.0 + float(istep)*step)*vec
    atoms.set_distance(index_1, index_2, (1.0 + float(istep)*step), fix=0)
#    Relax the structure
#    relax = QuasiNewton(atoms)
#    relax.run(fmax=0.05)
#
# Record the distance and energy
    dist.append(atom_dist(atoms, index_1, index_2))
    energy.append(atoms.get_potential_energy())

    print(f"step: %d dist: %10.6f Energy: %10.6f"   \
          % ( istep, dist[istep], energy[istep] ) )

    csv_row = [istep, dist[istep], energy[istep]]    
    
    if istep == 0:
# Update movie
        write('movie.xyz', atoms ,format='xyz')
# Update csv record of energy
        csv_file = open('energy_vs_dist.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_headers)
        csv_writer.writerow(csv_row)
        
        csv_file.close()

    else:
        write('movie.xyz', atoms ,format='xyz', append=True)
        # Update csv record of energy
        csv_file = open('energy_vs_dist.csv', 'a', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_row)
        csv_file.close()
      
print(f"")
print(f"Run completed.......")
print(f"")



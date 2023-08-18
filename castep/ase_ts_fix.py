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
# Bubble sort for atom list, this will sort atoms according to their
# co-ordinates in x,y,z according to the dir vector
#
def bubble_atoms(atoms, dir):
#
   natoms=len(atoms)
   if natoms < 0:
      print("ERROR: Cannot bubble sort an empty atom list")
      exit(0)
#
   print("In bubble_atoms with ", natoms, " atoms. Sorting in direction: ", dir)
#
# find the dot product of each position vector with the dir
#
   dot_list=[]
   for iatom in range(0,natoms):
      vec=[atoms.positions[iatom][0],atoms.positions[iatom][1],atoms.positions[iatom][2]]
      dot_list.append(np.dot(vec,dir))
##
   for iii in range(0,natoms):
      for jjj in range(iii+1,natoms):
#
         if dot_list[jjj] < dot_list[iii]:
            temp         =dot_list[iii].copy()
            dot_list[iii]=dot_list[jjj].copy()
            dot_list[jjj]=temp.copy()
##
            atoms.positions[[iii,jjj]]=atoms.positions[[jjj,iii]]
            atoms.symbols[[iii,jjj]]  =atoms.symbols[[jjj,iii]]
##
   return
#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
cell_file1 = sys.argv[1]                                           # Directory where we will do the calculations
cell_file2 = sys.argv[2]                                           # Directory where we will do the calculations
#splstr=cell_file.split('.')
#stem=splstr[0]
#cif_file="%s.cif" % stem
#
print(f"------------------------------------------------------------")
print(f"Comparing cell file          : %s" % cell_file1 )
print(f"with cell file               : %s" % cell_file2 )
print(f"------------------------------------------------------------")
#
atoms1=read(cell_file1)
atoms2=read(cell_file2)
nat1 = len(atoms1)
nat2 = len(atoms2)
#
# Note that species test is not robust, so P would be confused with Pd etc
#
species_list=[]
species_num=[]
for iatom in range(0,nat1):
    if iatom == 0:
       species_list.append(atoms1.symbols[iatom])
       species_num.append(1)
    else:
       have_this=False
       for ispec in range(0,len(species_list)):
           if atoms1.symbols[iatom] in species_list[ispec]:
               species_num[ispec]= species_num[ispec]+1
               have_this = True
#
       if not have_this:
           species_list.append(atoms1.symbols[iatom])
           species_num.append(1)
#
species_count=[]
print("\nFound species:")
for ispec in range(0,len(species_list)):
    species_count.append(0)
    print("%s : %d" % (species_list[ispec],species_num[ispec]))
#

lattice=atoms1.get_cell()
#
aaa=np.sqrt(np.dot(lattice[0],lattice[0]))
bbb=np.sqrt(np.dot(lattice[1],lattice[1]))
ccc=np.sqrt(np.dot(lattice[2],lattice[2]))

print("File 1 contains %d atoms" % nat1)
print("File 2 contains %d atoms" % nat2)
print(lattice)
print("aaa: %10.6f bbb: %10.6f ccc: %10.6f " % (aaa,bbb,ccc))

#write(cif_file, atoms ,format='cif')
#
# reorder atoms according to z-co-ordinate
#
for iatom in range(0,nat1):
#
# use atom as reference for z_fix
    if iatom == 251:
        fix_z =  atoms1.positions[iatom][2]+0.1
#
    for jatom in range(iatom+1,nat1):
        if atoms1.positions[iatom][2] > atoms1.positions[jatom][2]:
##
            atoms1.positions[[iatom,jatom]]=atoms1.positions[[jatom,iatom]]
            atoms1.symbols[[iatom,jatom]]  =atoms1.symbols[[jatom,iatom]]
##
            atoms2.positions[[iatom,jatom]]=atoms2.positions[[jatom,iatom]]
            atoms2.symbols[[iatom,jatom]]  =atoms2.symbols[[jatom,iatom]]
#
print("Will Freeze all non-Pd atoms below z =", fix_z)
#
# Report atom1 -> atom2 distances flag big deviations
fix_list=[]
fix_symbol=[]
fix_spec_no=[]
#
print("Atom list for file 1:")
for iatom in range(0,nat1):
    vec=[]
    fix_this=False
#
# Keep track of species counts
    for ispec in range(0,len(species_list)):
       if atoms1.symbols[iatom] in species_list[ispec]:
          this_spec=ispec
          species_count[ispec] = species_count[ispec] + 1
#
# To allow constraints make atom2 atoms that will be constrained copies of atom1 list
#
    if (atoms1.positions[iatom][2] < fix_z and "Pd" not in atoms1.symbols[iatom]):
        atoms2.positions[iatom] = atoms1.positions[iatom].copy()
        fix_this=True 
        fix_list.append(iatom)
        fix_symbol.append(atoms1.symbols[iatom])
        fix_spec_no.append(species_count[this_spec])
#
# vec is the vector from atom1 to atom2
#
    vec.append(atoms2.positions[iatom][0] - atoms1.positions[iatom][0])
    vec.append(atoms2.positions[iatom][1] - atoms1.positions[iatom][1])
    vec.append(atoms2.positions[iatom][2] - atoms1.positions[iatom][2])

    ddd = np.sqrt(np.dot(vec,vec))

    if ddd > 3.0:

       dcell = []
       dcell.append(np.dot(vec,lattice[0])/(aaa*ddd))
       dcell.append(np.dot(vec,lattice[1])/(bbb*ddd))
       dcell.append(np.dot(vec,lattice[2])/(ccc*ddd))
#
       print(f"%d %s %10.6f %10.6f %10.6f, dist: %10.6f Angs ---- big move....%6.2f %6.2f %6.2f" \
          % (iatom, atoms1.symbols[iatom], atoms1.positions[iatom][0], \
             atoms1.positions[iatom][1], atoms1.positions[iatom][2], ddd,\
             dcell[0],dcell[1],dcell[2]))
      
       shift=[0.0,0.0,0.0]
       for iii in range(0,3):
           absdcell=np.fabs(dcell[iii])
           if absdcell > 0.99:
               sign=dcell[iii]/absdcell
               for jjj in range(0,3):
                   shift[jjj]=shift[jjj]+sign*lattice[iii][jjj]
#
       for iii in range(0,3):
          atoms2.positions[iatom][iii]=atoms2.positions[iatom][iii]-shift[iii]
#
       print(f"Shifted atom2: %s %10.6f %10.6f %10.6f6.2f" \
          % ( atoms2.symbols[iatom], atoms2.positions[iatom][0], \
             atoms2.positions[iatom][1], atoms2.positions[iatom][2]))
#
    else:
       if fix_this:
          print(f"%d %s %10.6f %10.6f %10.6f, dist: %10.6f Angs fixed" \
             % (iatom, atoms1.symbols[iatom], atoms1.positions[iatom][0], \
                atoms1.positions[iatom][1], atoms1.positions[iatom][2], ddd))
       else:
          print(f"%d %s %10.6f %10.6f %10.6f, dist: %10.6f Angs" \
             % (iatom, atoms1.symbols[iatom], atoms1.positions[iatom][0], \
                atoms1.positions[iatom][1], atoms1.positions[iatom][2], ddd))


print("Fix list: ", fix_list)
#
# inter system vector
#
inter_vec=[]

for iatom in range(0,nat1):
    vec=[]
#
# vec is the vector from atom1 to atom2
#
    vec.append(atoms2.positions[iatom][0] - atoms1.positions[iatom][0])
    vec.append(atoms2.positions[iatom][1] - atoms1.positions[iatom][1])
    vec.append(atoms2.positions[iatom][2] - atoms1.positions[iatom][2])
#
    inter_vec.append(vec)
#
# Create CASTEP START, INTERMEDIATE and PRODUCT files
#
step=0.5
#
for istep in range(0,3):
#
   atoms=atoms1.copy()
#
   for iatom in range(0,nat1):
       atoms.positions[iatom][0]= atoms1.positions[iatom][0] + float(istep)*step*inter_vec[iatom][0]
       atoms.positions[iatom][1]= atoms1.positions[iatom][1] + float(istep)*step*inter_vec[iatom][1]
       atoms.positions[iatom][2]= atoms1.positions[iatom][2] + float(istep)*step*inter_vec[iatom][2]
#
   if istep == 0:
       write("new_start.cell", atoms , format='castep-cell')

   if istep == 1:
       write("new_intermed.cell", atoms, format='castep-cell')

   if istep == 2:
       write("new_end.cell", atoms, format='castep-cell')
#
# Add constraints to end cell file
#
cell_pntr=open("new_end.cell",'a')

cell_pntr.write("\n\n\%BLOCK CELL_CONSTRAINTS\n")
cell_pntr.write("   0    0   0\n   0    0   0\n")
cell_pntr.write("\%ENDBLOCK CELL_CONSTRAINTS\n")
#
cell_pntr.write("\n\n\%BLOCK IONIC_CONSTRAINTS\n")
#
iflist=0
for ifix in range(0,len(fix_list)):
   iflist+=1
   cell_pntr.write("   %d     %s    %d  1.0000000000000000    0.0000000000000000    0.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
   iflist+=1
   cell_pntr.write("   %d     %s    %d  0.0000000000000000    1.0000000000000000    0.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
   iflist+=1
   cell_pntr.write("   %d     %s    %d  0.0000000000000000    0.0000000000000000    1.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
#
cell_pntr.write("\%ENDBLOCK IONIC_CONSTRAINTS\n")
#
cell_pntr.close()
#
#
# Add constraints and create a trajectory for a linear interpolation check
#
cons = FixAtoms(fix_list)
atoms1.set_constraint(cons)
atoms2.set_constraint(cons)
#
step=0.1
traj_file="test.traj"
traj = Trajectory(traj_file, 'w')
#
for istep in range(0,11):
#
   atoms=atoms1.copy()
#
   for iatom in range(0,nat1):
       atoms.positions[iatom][0]= atoms1.positions[iatom][0] + float(istep)*step*inter_vec[iatom][0]
       atoms.positions[iatom][1]= atoms1.positions[iatom][1] + float(istep)*step*inter_vec[iatom][1]
       atoms.positions[iatom][2]= atoms1.positions[iatom][2] + float(istep)*step*inter_vec[iatom][2]
#
   traj.write(atoms)

   if istep == 0:
       write("new_start.in", atoms , format='aims')

   if istep == 5:
       write("new_intermed.in", atoms , format='aims')

   if istep == 10:
       write("new_end.in", atoms , format='aims')















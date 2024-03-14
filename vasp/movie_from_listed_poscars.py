# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:55:21 2020

Code to read a structure.castep file, i.e. output from CASTEP 
following an MD run so that a movie of the result can be visualised.

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

from ase.io import read, write, Trajectory
#
# Local definition of atom...atom distance measurement
#
def atom_dist(atoms, i1, i2):
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])

    size = np.sqrt(np.vdot(vec,vec))
    vec = vec/size
#
    return size
# Write arc file to get movie in Materials Studio
#
def write_arc(fp, atoms, iframe):
#
   if ( iframe == 0 ):
      fp.write("!BIOSYM archive 3\n")
      fp.write("PBC=ON\n")
#
   fp.write("Mode follower trajectory movie frame %d\n" % (iframe+1))
   fp.write("!DATE Mon Oct 12 12:17:26 1998\n")
#
# Write Cell Data as a b c alpha beta gamma
#
   abc=atoms.cell.cellpar()
   pbc_line="PBC {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} (P1)\n".\
           format(abc[0],abc[1],abc[2],abc[3],abc[4],abc[5])
   fp.write(pbc_line)
#
# Write atom data
# Au1      9.603440860   10.079051980    9.991388390 VASP X      ?       Au  0.000
#
   charges=[ 0 for i in range(0,len(atoms))]
   for iatom in range(0,len(atoms)):
       atom_line="{:4s}  {:14.9f} {:14.9f} {:14.9f} CAST X      ?       {:2s} {:6.3f}\n".\
                format(atoms.symbols[iatom], atoms.positions[iatom][0], \
                       atoms.positions[iatom][1], atoms.positions[iatom][2], \
                       atoms.symbols[iatom], charges[iatom])
       fp.write(atom_line)

   return
#
# Convert fractional to cartessian co-ordinates
#
def frac_to_cart(fvec, cell):
#
   vec=np.zeros(3)
   vec[0]=fvec[0]*cell[0][0]+fvec[1]*cell[1][0]+fvec[2]*cell[2][0]
   vec[1]=fvec[0]*cell[0][1]+fvec[1]*cell[1][1]+fvec[2]*cell[2][1]
   vec[2]=fvec[0]*cell[0][2]+fvec[1]*cell[1][2]+fvec[2]*cell[2][2]
#
   return vec
#
# Main code begins
#
filter=sys.argv[1] 
dir_list=os.listdir(".")
frame_files=[ fn for fn in dir_list if filter in fn ]
#print(dir_list)
num_files=len(frame_files)
print("Found %d frames" % num_files)

print(frame_files)
#
# Order files by the number given after "_"
#
num=[]
for file in frame_files:
   splstr=file.split("_")
#
   num.append(int(splstr[-1]))
#
#
# bubble num list and move frame_files into order
for iii in range(0,num_files):
   for jjj in range(iii+1,num_files):
#
      if num[jjj] < num[iii]:
         temp    =num[iii]
         num[iii]=num[jjj]
         num[jjj]=temp
##
         temp = frame_files[iii]
         frame_files[iii]=frame_files[jjj]
         frame_files[jjj]=temp

print("After sorting....")
print(num)
print(frame_files)
#
traj_file="movie.traj"
#
print(f"------------------------------------------------------------")
print(f"POSCAR files for movie  : ")
for iii in range(0,num_files):
    print(f"                 %s" % frame_files[iii] )
print(f"Trajectory file              : %s" % traj_file )
print(f"------------------------------------------------------------")
#
traj = Trajectory(traj_file, 'w')

for icell in range(0,num_files):
   print("Adding %s to trajectory" % frame_files[icell])
   atoms=read(frame_files[icell])
   traj.write(atoms)
#

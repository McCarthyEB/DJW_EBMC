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
# Read an MD .castep output to extract structure
#
def read_castep(fp, traj):
#
    lines = fp.readlines()
    nlines = len(lines)
    icell=0
    cell_line=0
    contents_line=0
    cell=np.zeros((3,3))
#
    have_unit_cell=False
    have_atoms=False
    first_struct=True
    atom_offset=10
    atoms=Atom()
    natoms=0
    num_atoms=-1
#
    for iline in range(nlines):
#
# Look for the stress lines
#
        line=lines[iline]
#        print(line)
        linesplit=line.split()
        nwords=len(linesplit)
#
# Read in cell vectors
#
        if ( nwords == 2 ):
           if ( 'Unit' in linesplit[0] and 'Cell' in linesplit[1] ):
             have_unit_cell=True
             cell_line=iline
#             print("Found Unit Cell cell_line = {}".format(cell_line))
#
           elif  ( 'Cell' in linesplit[0] and 'Contents' in linesplit[1] ):
             have_atoms=True
             contents_line=iline
#             print("Found Cell Contents contents_line = {}".format(contents_line))
#
        if ( have_unit_cell ):
           if ( iline-cell_line > 5):
              have_unit_cell=False
              icell=0
           elif ( iline-cell_line > 2):
              cell[icell][0] = float(linesplit[0])
              cell[icell][1] = float(linesplit[1])
              cell[icell][2] = float(linesplit[2])
              icell+=1
#
# Read in atoms first structure
#  
        if ( have_atoms ):
#           print("have_atoms iline-contents = {}".format(iline-contents_line))
#
           if (iline-contents_line == 2):
              if (nwords > 1):
                 first_struct=True
                 atom_offset=10
#                 print("first")
              else:
#                 print("not first")
                 atom_offset=7
                 first_struct=False
#
#           print("Atom offset now {}".format(atom_offset))
#
           if (first_struct and iline-contents_line == 2):
              num_atoms=int(linesplit[7])
#              print("Will read {} atoms".format(num_atoms))
#
           elif ( iline-contents_line == atom_offset):
              label=linesplit[1]
              fracs=[float(linesplit[3]),float(linesplit[4]),float(linesplit[5])]
#              print("label: {}, coords: {}".format(label,fracs))
              coords=frac_to_cart(fracs, cell)
#              print("label: {}, coords: {}".format(label,coords))
#
              new_atom=Atom(label,coords)
              atoms=Atoms([new_atom], cell=cell)
              natoms+=1
#
           elif (iline-contents_line > atom_offset and natoms < num_atoms): 
              label=linesplit[1]
              fracs=[float(linesplit[3]),float(linesplit[4]),float(linesplit[5])]
#              print("label: {}, coords: {}".format(label,fracs))
              coords=frac_to_cart(fracs, cell)
#              print("label: {}, coords: {}".format(label,coords))
              new_atom=Atom(label,coords)
              atoms.append(new_atom)
              natoms+=1
#
           elif ( natoms == num_atoms):
              have_atoms=False
              natoms=0
              first_struct=False
#              print("Writing Frame")
              traj.write(atoms)
#
#    print(cell)           
    
    return atoms
#
# Write arc file to get movie in Materials Studio
#
def write_arc(fp, atoms, iframe):
#
   if ( iframe == 0 ):
      fp.write("!BIOSYM archive 3\n")
      fp.write("PBC=ON\n")
#
   fp.write("MD trajectory movie frame %d\n" % iframe)
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
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
castep_file = sys.argv[1]              # File name supplied from command line. 
splstr=castep_file.split('.')
stem=splstr[0]
traj_file="%s.traj" % stem
movie_file="%s.xyz" % stem
arc_file="%s.arc" % stem
#
print(f"------------------------------------------------------------")
print(f"CASTEP output file to read   : %s" % castep_file )
print(f"Stem                         : %s" % stem     )
print(f"Trajectory file              : %s" % traj_file )
print(f"Movie file of trajectory     : %s" % movie_file )
print(f"Materials Studio version     : %s" % arc_file )
print(f"------------------------------------------------------------")
#
file_pntr=open(castep_file, 'rt')
arc_pntr=open(arc_file, 'wt')
traj = Trajectory(traj_file, 'w')

atoms=read_castep(file_pntr, traj)

#write(castep_cellfile, atoms ,format='castep-cell')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
traj_range="%s@:" % traj_file
configs = read(traj_range) 
print("configs length = {}".format(len(configs)))
#
iframe=0
for iframe in range(0,len(configs)):
   atoms=configs[iframe]
   if (iframe==0):
      write(movie_file, atoms, format='xyz')
   else:
      write(movie_file, atoms, format='xyz', append=True)
   write_arc(arc_pntr, atoms, iframe)
#
file_pntr.close()
arc_pntr.close()

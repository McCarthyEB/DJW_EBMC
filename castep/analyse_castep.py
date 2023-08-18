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
# Read a .castep output to extract structure
# and Forces, calculate constrained RMS force         
# May 2023.
#
def read_castep(fp):
#
    is_trans_state=False
    is_geom_opt   =False
    have_QST=False
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
    have_forces=False
    constrained=False
    atom_offset=10
    atoms=Atom()
    natoms=0
    nforces=0
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
# Check for type of calculation
#
        if ( nwords >= 6 and 'type' in linesplit[0] ):
            if ( 'transition' in linesplit[4] and 'state' in linesplit[5] ):
                print("This is a transition state calculation.")
                is_trans_state = True
       
            elif ( 'geometry' in linesplit[4] and 'optimization' in linesplit[5] ):
                print("This is a geometry optimization.")
                is_geom_opt   = True 
#
# Read in cell vectors and flag start of atom list
#
        if ( nwords == 2 ):
           if ( 'Unit' in linesplit[0] and 'Cell' in linesplit[1] ):
             have_unit_cell=True
             cell_line=iline
             print("Found Unit Cell cell_line = {}".format(cell_line))
#
           elif  ( 'Cell' in linesplit[0] and 'Contents' in linesplit[1] ):
             have_atoms=True
             contents_line=iline
             print("Found Cell Contents contents_line = {}".format(contents_line))
#
# Look for forces
#
        elif ( nwords == 4 ):
          if "Forces" in linesplit[2] :
             force_offset=6
             nforces=0
             nforces_cons=0
             fsq_tot=0.0
             if "Uncon" in linesplit[1] :
                have_forces=True
                constrained=False
             elif "Const" in linesplit[1] :
                have_forces=True
                constrained=True 
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
# For Transition state searches pick up structures after QST maximum found
#
        if (is_trans_state and nwords == 3 and 'QST' in linesplit[0] and 'Found' in linesplit[2]):
            have_atoms=True
            have_QST=True
            contents_line=iline
#
# Read in forces
#
        if have_forces:
          if force_offset > 0:
             force_offset-=1
          else:
#             print("Force line:")
#             print("%s" % line)
             nforces+=1
             fx = float(linesplit[3])
             fy = float(linesplit[4])
             fz = float(linesplit[5])
#
             fsq = fx*fx + fy*fy + fz*fz
             fsq_tot += fsq
               
             print("fsq: %10.6f" % fsq)
#
          if nforces == num_atoms:
             have_forces=False
             fsq_tot=fsq_tot/float(nforces)
             frms = np.sqrt(fsq_tot)
             if constrained:
                print("RMS Force, ( constrained ) : %10.6f" % frms)
             else:
                print("RMS Force, (   all atoms ) : %10.6f" % frms)
                
#
# Read in atoms first structure
#  
        if ( have_atoms ):
#           print("have_atoms iline-contents = {}".format(iline-contents_line))
#
           if (iline-contents_line == 2):
              if (nwords > 1 and 'Total' in linesplit[0] ):
                 first_struct=True
#                 print("first")
                 if is_trans_state:
                   atom_offset=12
                 else:
                   atom_offset=10
              else:
#                 print("not first")
                 if is_trans_state:
                   atom_offset=16
                 else:
                   atom_offset=7
                 first_struct=False
#
#           print("Atom offset now {}".format(atom_offset))
#
           if (first_struct and iline-contents_line == 2):
              num_atoms=int(linesplit[7])
              print(line)
              print("Will read {} atoms".format(num_atoms))
#
           elif ( iline-contents_line == atom_offset):
              print("First atom line:")
              print(line)
#
              if have_QST:
                label=linesplit[0][1:]
                fracs=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]
              else:
                label=linesplit[1]
                fracs=[float(linesplit[3]),float(linesplit[4]),float(linesplit[5])]
#              print("label: {}, coords: {}".format(label,fracs))
#
              coords=frac_to_cart(fracs, cell)
#              print("label: {}, coords: {}".format(label,coords))
#
              new_atom=Atom(label,coords)
              atoms=Atoms([new_atom], cell=cell)
              natoms+=1
#
           elif (iline-contents_line > atom_offset and natoms < num_atoms): 
#
              if have_QST:
                label=linesplit[0][1:]
                fracs=[float(linesplit[2]),float(linesplit[3]),float(linesplit[4])]
              else:
                label=linesplit[1]
                fracs=[float(linesplit[3]),float(linesplit[4]),float(linesplit[5])]
#
#              print("label: {}, coords: {}".format(label,fracs))
              coords=frac_to_cart(fracs, cell)
#              print("label: {}, coords: {}".format(label,coords))
              new_atom=Atom(label,coords)
              atoms.append(new_atom)
              natoms+=1
#
           elif ( natoms == num_atoms):
              have_atoms=False
              have_QST=False
              natoms=0
              first_struct=False
#              print("Writing Frame")
#
#    print(cell)           
    
    return atoms
#
# Main code begins
#
nargs=len(sys.argv)
print(f"Arguements passed {len(sys.argv)}")
#
if nargs > 1: 
   castep_file = sys.argv[1]              # File name supplied from command line. 
   indices=[-1,-1]
   if nargs > 3:
     indices=[int(sys.argv[2]),int(sys.argv[3])]
else:
   print("ERROR: Need at least one arguement giving filename")
   print("ERROR: with three arguments also measure atom separation")
   exit(0)
#
splstr=castep_file.split('.')
stem=splstr[0]
#
print(f"------------------------------------------------------------")
print(f"CASTEP output file to read   : %s" % castep_file )
print(f"Stem                         : %s" % stem     )
#
file_pntr=open(castep_file, 'rt')

atoms=read_castep(file_pntr)
#
if indices[0] > -1 and  indices[1] > -1:
   print(f"Will report distance between atoms : %s %d and %s %d" % \
              ( atoms.symbols[indices[0]], indices[0],\
                atoms.symbols[indices[1]], indices[1]))

print(f"------------------------------------------------------------")
#

#write(castep_cellfile, atoms ,format='castep-cell')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
if indices[0] > -1 and  indices[1] > -1:
   ddd = atoms.get_distance(indices[0],indices[1],mic=True)
   print("Atom %s %d to Atom %s %d = %10.6f Angs." % ( atoms.symbols[indices[0]], indices[0],\
                atoms.symbols[indices[1]], indices[1], ddd))
#
file_pntr.close()

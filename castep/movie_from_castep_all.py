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
# Adapted to also take in Transition State Search info
# May 2023.
#
def read_castep(fp, traj):
#
    is_trans_state=False
    is_geom_opt   =False
    have_QST=False
    ts_list=[]
    bar_react_list=[]
    bar_prod_list=[]
    loc_list=[]
    loc_change=[]
    qst_rms_force=[]
    iframe=0
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
# Read in cell vectors
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
            ts_list.append(iframe)
#
# Read in barriers and RMS force
#
        if have_QST:
           if (nwords == 5 and 'Barrier' in  linesplit[0]):
              if 'reactant' in linesplit[2]:
                 bar_react_list.append(float(linesplit[3]))
              if 'product' in linesplit[2]:
                 bar_prod_list.append(float(linesplit[3]))
           if (nwords == 5 and 'Location' in linesplit[0]):
              loc_list.append(float(linesplit[4]))
#
# Pick up RMS force for TS calc, last item before resetting the flag
#
           if (nwords == 6 and 'Local' in linesplit[0] and 'QST' in linesplit[3]):
              loc_change.append(float(linesplit[5]))
           if (nwords == 7 and 'RMS' in linesplit[0] and 'QST' in linesplit[3]):
              qst_rms_force.append(float(linesplit[5]))
              have_QST=False
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
              natoms=0
              first_struct=False
#              print("Writing Frame")
              iframe+=1
              traj.write(atoms)
#
#    print(cell)           
    print("Leaving read_castep....")
    print("Recorded %d frames" % iframe)
#
    if (is_trans_state):
       num_ts_states= len(ts_list)
       print("Of which %d are transition states" % len(ts_list))
       print("Found %d forward barriers" % len(bar_react_list)) 
       print("Found %d reverse barriers" % len(bar_prod_list)) 
       print("Found %d barrier locations" % len(loc_list)) 
       print("Found %d local change entries" % len(loc_change)) 
       print("Found %d QST RMS force entries" % len(qst_rms_force)) 
#
       if ( len(bar_react_list) != num_ts_states):
          print("ERROR: number of TS state energies does not match, check castep file")
#
       if ( len(bar_prod_list) != num_ts_states):
          print("ERROR: number of TS state energies does not match, check castep file")
#
       if ( len(loc_list) != num_ts_states):
          print("ERROR: number of TS state locations reported does not match, check castep file")
#
       if ( len(loc_change) != num_ts_states):
          print("ERROR: number of local change reported does not match, check castep file")
#
       if ( len(qst_rms_force) != num_ts_states):
          print("ERROR: number of QST RMS force reported does not match, check castep file")
#
# Write csv summary file
#    
       csv_pntr=open("TS_summary.csv", 'wt')
       title_line="TS num, frame num in movie, location, Barrier from react, Barrier from prod,"
       title_line+= "QST RMS force, QST local change\n"
       csv_pntr.write(title_line)
#  
       for iii in range(0, num_ts_states):
          data_line="%d, %d, %10.6f,  %10.6f,  %10.6f,  %10.6f,  %10.6f\n"\
                    % (iii+1, ts_list[iii]+1, loc_list[iii], bar_react_list[iii],\
                       bar_prod_list[iii], qst_rms_force[iii], loc_change[iii])
          csv_pntr.write(data_line)
#   
       csv_pntr.close()
#   
    return atoms, ts_list
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
cell_stem="ts_guess_"
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

atoms,ts_list=read_castep(file_pntr, traj)

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
ts_num=1
for iframe in range(0,len(configs)):
#
   atoms=configs[iframe]
#
   if iframe in ts_list:
      ts_guess_file="%s%d.cell" % (cell_stem, ts_num)
      write(ts_guess_file, atoms , positions_frac=False, format='castep-cell')
      ts_num+=1
#
   if (iframe==0):
      write(movie_file, atoms, format='xyz')
   else:
      write(movie_file, atoms, format='xyz', append=True)
   write_arc(arc_pntr, atoms, iframe)
#
#write("last_TS.cell", atoms , positions_frac=False, format='castep-cell')
#
file_pntr.close()
arc_pntr.close()

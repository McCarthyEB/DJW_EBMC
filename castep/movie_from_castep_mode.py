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
# Read the lattice vectors and atom coordinates from a cell file 
# and make an atoms object
#
def read_cell_djw(fp):
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
   have_cell_cons=False
   have_kpoint_grid=False
   have_fix_com=False
   first_atom=True
   natoms=0
   num_atoms=-1
   num_intermed_atoms=-1
   num_prod_atoms=-1
   num_cons=-1
#
   spec_mass_lab=[]
   spec_mass=[]
#
   spec_pot_lab=[]
   spec_pot=[]
#
   spec_lcao_lab=[]
   spec_lcao=[]
#
#
   read_lines=False
   read_cell=False
   read_atoms=False
   read_cell_cons=False
   read_ion_cons=False
   read_spec_mass=False
   read_spec_pot=False
   read_spec_lcao=False
#
   have_cell=False
   have_atoms=False
   have_cell_cons=False
   have_ion_cons=False
   have_masses=False
   have_pots=False
   have_lcao=False
   have_kpoint_grid=False
   have_fix_com=False
#
   for iline in range(nlines):
#
# Look for the stress lines
#
      line=lines[iline]
      linesplit=line.split()
      nwords=len(linesplit)
#
      if ( nwords >= 2 and  "endblock" in linesplit[0].lower() ):
#         print("Found ENDBLOCK: %s" % line)
         read_lines=False
         read_cell=False
         read_atoms=False
         read_cell_cons=False
         read_ion_cons=False
         read_spec_mass=False
         read_spec_pot=False
         read_spec_lcao=False
#
      if read_cell and nwords == 3:
#         print("Reading cell line: %s" % line)
         cell[icell][0] = float(linesplit[0])
         cell[icell][1] = float(linesplit[1])
         cell[icell][2] = float(linesplit[2])
         icell+=1
#
      elif read_atoms:
#
         label=linesplit[0]
#
         if is_frac:
            fracs=[float(linesplit[1]),float(linesplit[2]),float(linesplit[3])]
            coords=frac_to_cart(fracs, cell)
#
         else:
            coords=[float(linesplit[1]),float(linesplit[2]),float(linesplit[3])]
#
         new_atom=Atom(label,coords)
#
         if first_atom:
#            print("Reading atoms")
#            print("First atom line:")
#            print(line)
            first_atom=False
#
            if is_intermed:
              atoms_intermed=Atoms([new_atom], cell=cell)
              num_intermed_atoms+=1
            elif is_prod:
              atoms_prod=Atoms([new_atom], cell=cell)
              num_prod_atoms+=1
            else:
              atoms=Atoms([new_atom], cell=cell)
              num_atoms+=1
         else:
            if is_intermed:
              atoms_intermed.append(new_atom)
              num_intermed_atoms+=1
            elif is_prod:
              atoms_prod.append(new_atom)
              num_prod_atoms+=1
            else:
              atoms.append(new_atom)
              num_atoms+=1
#
      elif read_cell_cons:
         have_cell_cons=True
#     
         cell_cons.append(int(linesplit[0]))
         cell_cons.append(int(linesplit[1]))
         cell_cons.append(int(linesplit[2]))
#
      elif read_spec_mass and nwords == 2:    
        have_masses=True
#        print("Reading masses.....")
        spec_mass_lab.append(linesplit[0])
        spec_mass.append(float(linesplit[1]))
#
      elif read_spec_pot:
        have_pots=True
        spec_pot_lab.append(linesplit[0])
        spec_pot.append(linesplit[1])
#
      elif read_spec_lcao:
        have_lcao=True
        spec_lcao_lab.append(linesplit[0])
        spec_lcao.append(int(linesplit[1]))
#
#
# Check for BLOCK flags                
#
      if ( nwords >= 2 and 'block' in linesplit[0].lower() and 'end' not in linesplit[0].lower() ):
         read_lines=True
#         print("Looking at BLOCK : %s" % line)
#
# Read second word    
#
         if ( read_lines and "lattice_cart" in linesplit[1].lower() ):
           read_cell=True
           icell=0
#
         elif ( read_lines and "positions" in linesplit[1].lower() ):
           read_atoms=True
           is_frac=False
           is_intermed=False
           is_prod=False
           if "frac" in linesplit[1].lower():
              is_frac=True
           if "intermediate" in linesplit[1].lower():
              is_intermed=True
              first_atom=True
           elif "PRODUCT" in linesplit[1]:
              is_prod=True
              first_atom=True
#
         elif ( read_lines and "cell_constraints" in linesplit[1].lower() ):
           read_cell_cons=True
           cell_cons=[]
#
         elif ( read_lines and "ionic_constraints" in linesplit[1].lower() ):
           read_ion_cons=True
#
         elif ( read_lines and "species_mass" in linesplit[1].lower() ):
#           print("Will read masses...........")
           read_spec_mass=True
#
         elif ( read_lines and "species_pot" in linesplit[1].lower() ):
           read_spec_pot=True
#
         elif ( read_lines and "species_lcao_states" in linesplit[1].lower() ):
           read_spec_lcao=True
#
# Pick up one liners
#
#kpoint_mp_grid 1 1 1
#or kpoint_mp_grid : 1 1 1
#FIX_COM : false
      elif ( nwords >= 2 and 'kpoint_mp_grid' in linesplit[0].lower()):
         have_kpoint_grid=True
         if ":" in linesplit[1]:
            kpoint_mp=[int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
         else:
            kpoint_mp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
#
      elif ( nwords >= 2 and 'fix_com' in linesplit[0].lower()):
         have_fix_com=True
         fix_com=linesplit[2]
#
   return atoms
#
#
# Main code begins
#
filter="cell"
dir_list=os.listdir(".")
cell_files=[ fn for fn in dir_list if fn.endswith(filter) and "-" in fn ]
#print(dir_list)
num_files=len(cell_files)
#print("Found %d mode frames" % num_files)

#print(cell_files)
#
# Order files by the number given between "-" and "."
#
num=[]
for file in cell_files:
   splstr=file.split(".")
#
# deal with negative move frames
   if "-m" in splstr[0]:
      num.append(-int(splstr[0].split("-m")[1]))
   else:
      num.append(int(splstr[0].split("-")[1]))
#
#
# bubble num list and move cell_files into order
for iii in range(0,num_files):
   for jjj in range(iii+1,num_files):
#
      if num[jjj] < num[iii]:
         temp    =num[iii]
         num[iii]=num[jjj]
         num[jjj]=temp
##
         temp = cell_files[iii]
         cell_files[iii]=cell_files[jjj]
         cell_files[jjj]=temp

#print("After sorting....")
#print(num)
#print(cell_files)
#
stem=splstr[0].split("-")[0]
traj_file="%s.traj" % stem
movie_file="%s.xyz" % stem
arc_file="%s.arc" % stem
#
print(f"------------------------------------------------------------")
print(f"CASTEP cell files for movie  : ")
for iii in range(0,num_files):
    print(f"                 %s" % cell_files[iii] )
print(f"Stem                         : %s" % stem     )
print(f"Trajectory file              : %s" % traj_file )
print(f"Movie file of trajectory     : %s" % movie_file )
print(f"Materials Studio version     : %s" % arc_file )
print(f"------------------------------------------------------------")
#
arc_pntr=open(arc_file, 'wt')
traj = Trajectory(traj_file, 'w')

for icell in range(0,num_files):
   print("Adding %s to trajectory" % cell_files[icell])
   file_pntr=open(cell_files[icell], 'rt')
   atoms=read_cell_djw(file_pntr)
   traj.write(atoms)
#
   write_arc(arc_pntr, atoms, icell)

   file_pntr.close()
#
arc_pntr.close()

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:55:21 2020

Code to read a castep cell file, and remove the bottom section.
constraints set for geometry runs are reformatted for frequency calculation 

@author: dave
"""
import numpy as np
import sys
import csv
import os
#
from ase import Atom
from ase import Atoms

from ase.io import read, write, Trajectory
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
print("Starting conversion to freq format")
#
cell_file = sys.argv[1]              # File name supplied from command line.
splstr=cell_file.split('.')
stem=splstr[0]
#
print(f"------------------------------------------------------------")
print(f"CASTEP cell file to read   : %s" % cell_file )
print(f"Stem                       : %s" % stem     )
print(f"------------------------------------------------------------")
#
fp=open(cell_file, 'rt')
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
num_freq_cons=-1
freq_cons_indx=[]
freq_cons_lab=[]
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
   if ( nwords >= 2 and  "ENDBLOCK" in linesplit[0] ):
      print("Found ENDBLOCK: %s" % line)
      read_lines=False
      read_cell=False
      read_atoms=False
      read_cell_cons=False
      read_ion_cons=False
      read_spec_mass=False
      read_spec_pot=False
      read_spec_lcao=False
#
   if read_cell:
      print("Reading cell line: %s" % line)
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
         print("Reading atoms")
         print("First atom line:")
         print(line)
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
   elif read_ion_cons:
      num_cons=int(linesplit[0])
      lab=linesplit[1]
      indx=int(linesplit[2])
#
      if num_freq_cons < 0:
         print("Reading ion constraints")
         num_freq_cons+=1
         freq_cons_indx.append(indx)
         freq_cons_lab.append(lab)
         print("Adding %s %d" % (lab,indx))
      else:
         read_this=True
         for icon in range(0,num_freq_cons+1):
            if indx == freq_cons_indx[icon] and lab in freq_cons_lab[icon]:
               read_this=False 
#
         if read_this:
            num_freq_cons+=1
            freq_cons_indx.append(indx)
            freq_cons_lab.append(lab)
            print("Adding %s %d" % (lab,indx))
#
   elif read_spec_mass:    
     have_masses=True
     print("Reading masses.....")
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
   if ( nwords >= 2 and 'BLOCK' in linesplit[0] and 'END' not in linesplit[0] ):
      read_lines=True
      print("Looking at BLOCK : %s" % line)
#
# Read second word    
#
      if ( read_lines and "LATTICE_CART" in linesplit[1] ):
        read_cell=True
        icell=0
#
      elif ( read_lines and "POSITIONS" in linesplit[1] ):
        read_atoms=True
        is_frac=False
        is_intermed=False
        is_prod=False
        if "FRAC" in linesplit[1]:
           is_frac=True
        if "INTERMEDIATE" in linesplit[1]:
           is_intermed=True
           first_atom=True
        elif "PRODUCT" in linesplit[1]:
           is_prod=True
           first_atom=True
#
      elif ( read_lines and "CELL_CONSTRAINTS" in linesplit[1] ):
        read_cell_cons=True
        cell_cons=[]
#
      elif ( read_lines and "IONIC_CONSTRAINTS" in linesplit[1] ):
        read_ion_cons=True
#
      elif ( read_lines and "SPECIES_MASS" in linesplit[1] ):
        print("Will read masses...........")
        read_spec_mass=True
#
      elif ( read_lines and "SPECIES_POT" in linesplit[1] ):
        read_spec_pot=True
#
      elif ( read_lines and "SPECIES_LCAO_STATES" in linesplit[1] ):
        read_spec_lcao=True
#
# Pick up one liners
#
#kpoint_mp_grid 1 1 1
#FIX_COM : false
   elif ( nwords >= 2 and 'kpoint_mp_grid' in linesplit[0]):
      have_kpoint_grid=True
      kpoint_mp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
#
   elif ( nwords >= 2 and 'FIX_COM' in linesplit[0]):
      have_fix_com=True
      fix_com=linesplit[2]
#
# make num_freq_cons a true count rather than highest index
if num_cons >= 0:
    num_freq_cons+=1
#
if num_atoms >= 0:
    num_atoms+=1
#
if num_intermed_atoms >= 0:
    num_intermed_atoms+=1
#
if num_prod_atoms >= 0:
    num_prod_atoms+=1

print("All read:")
print("%d atoms in standard, %d atoms in intermediate and %d atoms in product" % (num_atoms, num_intermed_atoms, num_prod_atoms))
print("%d atoms with %d geom cons and %d freq cons" % (len(atoms),num_cons,num_freq_cons))

#
# Write frequency cell file
#
test_cell_file="test.cell"
fp_out=open(test_cell_file, 'w')
#
fp_out.write("######################################################\n")
fp_out.write("### Cell file for frequency calc converted from ######\n")
fp_out.write("### constrained geometry optimisation           ######\n")
fp_out.write("######################################################\n\n")
#
fp_out.write("\n%BLOCK LATTICE_CART\n")
for icell in range(0,3):
   fp_out.write("%10.6f  %10.6f  %10.6f\n" % ( cell[icell][0], cell[icell][1], cell[icell][2]) )
fp_out.write("%ENDBLOCK LATTICE_CART\n\n")
#
fp_out.write("\n%BLOCK POSITIONS_ABS\n")
for iatom in range(0,num_atoms):
   fp_out.write("%s  %10.6f  %10.6f  %10.6f\n" % ( atoms.symbols[iatom], atoms.positions[iatom][0],\
                                                   atoms.positions[iatom][1], atoms.positions[iatom][2] ))
fp_out.write("%ENDBLOCK POSITIONS_ABS\n\n")
#
if have_cell_cons:
   fp_out.write("\n%BLOCK CELL_CONSTRAINTS\n")
   fp_out.write("%d    %d    %d\n" % (cell_cons[0],cell_cons[1],cell_cons[2]))
   fp_out.write("%d    %d    %d\n" % (cell_cons[3],cell_cons[4],cell_cons[5]))
   fp_out.write("%ENDBLOCK CELL_CONSTRAINTS\n\n")
#
if num_freq_cons >= 0:
   fp_out.write("\n%BLOCK IONIC_CONSTRAINTS\n")
   for icons in range(0,num_freq_cons):
      fp_out.write("%d   %s   %d   1  1  1\n" % ( icons+1, freq_cons_lab[icons], freq_cons_indx[icons] ))
   fp_out.write("%ENDBLOCK IONIC_CONSTRAINTS\n\n")
#
if have_fix_com:
   fp_out.write("\nFIX_COM : %s\n\n" % fix_com)
#
if have_kpoint_grid:
   fp_out.write("\nkpoint_mp_grid %d %d %d\n\n" % (kpoint_mp[0],kpoint_mp[1],kpoint_mp[2]))
#
if have_masses:
   fp_out.write("\n%BLOCK SPECIES_MASS\n")
   for imass in range(0, len(spec_mass)):
      fp_out.write("     %s    %13.10f\n" % (spec_mass_lab[imass],spec_mass[imass]))
   fp_out.write("%ENDBLOCK SPECIES_MASS\n\n")
      
if have_pots:
   fp_out.write("\n%BLOCK SPECIES_POT\n")
   for ipot in range(0, len(spec_pot)):
      fp_out.write("     %s    %s\n" % (spec_pot_lab[ipot],spec_pot[ipot]))
   fp_out.write("%ENDBLOCK SPECIES_POT\n\n")

if have_lcao:
   fp_out.write("\n%BLOCK SPECIES_LCAO_STATES\n")
   for ilcao in range(0, len(spec_lcao)):
      fp_out.write("     %s    %d\n" % (spec_lcao_lab[ilcao],spec_lcao[ilcao]))
   fp_out.write("%ENDBLOCK SPECIES_LCAO_STATES\n\n")

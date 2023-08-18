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
   return atoms
#
# Write a simple cell file 
#
def write_cell_djw(fp,atoms):
#
   have_cell_cons=False
   have_ion_cons=False
   have_masses=False
   have_pots=False
   have_lcao=False
   have_kpoint_grid=False
   have_fix_com=False
   num_freq_cons=-1
#
   fp.write("######################################################\n")
   fp.write("### Cell file written but write_cell_djw        ######\n")
   fp.write("######################################################\n\n")
#
   cell=atoms.get_cell()
   num_atoms=len(atoms)
   fp.write("\n%BLOCK LATTICE_CART\n")
   for icell in range(0,3):
      fp.write("%10.6f  %10.6f  %10.6f\n" % ( cell[icell][0], cell[icell][1], cell[icell][2]) )
   fp.write("%ENDBLOCK LATTICE_CART\n\n")
#
   fp.write("\n%BLOCK POSITIONS_ABS\n")
   for iatom in range(0,num_atoms):
      fp.write("%s  %10.6f  %10.6f  %10.6f\n" % ( atoms.symbols[iatom], atoms.positions[iatom][0],\
                                                   atoms.positions[iatom][1], atoms.positions[iatom][2] ))
   fp.write("%ENDBLOCK POSITIONS_ABS\n\n")
#
   if have_cell_cons:
      fp.write("\n%BLOCK CELL_CONSTRAINTS\n")
      fp.write("%d    %d    %d\n" % (cell_cons[0],cell_cons[1],cell_cons[2]))
      fp.write("%d    %d    %d\n" % (cell_cons[3],cell_cons[4],cell_cons[5]))
      fp.write("%ENDBLOCK CELL_CONSTRAINTS\n\n")
#
   if num_freq_cons >= 0:
      fp.write("\n%BLOCK IONIC_CONSTRAINTS\n")
      for icons in range(0,num_freq_cons):
         fp.write("%d   %s   %d   1  1  1\n" % ( icons+1, freq_cons_lab[icons], freq_cons_indx[icons] ))
      fp.write("%ENDBLOCK IONIC_CONSTRAINTS\n\n")
#
   if have_fix_com:
      fp.write("\nFIX_COM : %s\n\n" % fix_com)
#
   if have_kpoint_grid:
      fp.write("\nkpoint_mp_grid %d %d %d\n\n" % (kpoint_mp[0],kpoint_mp[1],kpoint_mp[2]))
#
   if have_masses:
      fp.write("\n%BLOCK SPECIES_MASS\n")
      for imass in range(0, len(spec_mass)):
         fp.write("     %s    %13.10f\n" % (spec_mass_lab[imass],spec_mass[imass]))
      fp.write("%ENDBLOCK SPECIES_MASS\n\n")
      
   if have_pots:
      fp.write("\n%BLOCK SPECIES_POT\n")
      for ipot in range(0, len(spec_pot)):
         fp.write("     %s    %s\n" % (spec_pot_lab[ipot],spec_pot[ipot]))
      fp.write("%ENDBLOCK SPECIES_POT\n\n")

   if have_lcao:
      fp.write("\n%BLOCK SPECIES_LCAO_STATES\n")
      for ilcao in range(0, len(spec_lcao)):
         fp.write("     %s    %d\n" % (spec_lcao_lab[ilcao],spec_lcao[ilcao]))
      fp.write("%ENDBLOCK SPECIES_LCAO_STATES\n\n")

   return
#
#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
cell_file = sys.argv[1]                                           # File to be converted
fix_z = float(sys.argv[2])                                        # defined z-co-ordinate level
#splstr=cell_file.split('.')
#stem=splstr[0]
#cif_file="%s.cif" % stem
#
print(f"------------------------------------------------------------")
print(f"Fixing atoms in cell file    : %s" % cell_file )
print(f"will fix atoms with Z below  : %10.6f" % fix_z )
print(f"------------------------------------------------------------")
#
fp=open(cell_file,"r")
atoms=read_cell_djw(fp)
natoms = len(atoms)
#
# Note that species test is not robust, so P would be confused with Pd etc
#
species_list=[]
species_num=[]
for iatom in range(0,natoms):
    if iatom == 0:
       species_list.append(atoms.symbols[iatom])
       species_num.append(1)
    else:
       have_this=False
       for ispec in range(0,len(species_list)):
           if atoms.symbols[iatom] in species_list[ispec]:
               species_num[ispec]= species_num[ispec]+1
               have_this = True
#
       if not have_this:
           species_list.append(atoms.symbols[iatom])
           species_num.append(1)
#
species_count=[]
print("\nFound species:")
for ispec in range(0,len(species_list)):
    species_count.append(0)
    print("%s : %d" % (species_list[ispec],species_num[ispec]))
#

lattice=atoms.get_cell()
#
aaa=np.sqrt(np.dot(lattice[0],lattice[0]))
bbb=np.sqrt(np.dot(lattice[1],lattice[1]))
ccc=np.sqrt(np.dot(lattice[2],lattice[2]))

print("File contains %d atoms" % natoms)
print("lattice:")
print("aaa: %10.6f bbb: %10.6f ccc: %10.6f " % (aaa,bbb,ccc))

#write(cif_file, atoms ,format='cif')
#
# reorder atoms according to z-co-ordinate
#
for iatom in range(0,natoms):
#
    for jatom in range(iatom+1,natoms):
        if atoms.positions[iatom][2] > atoms.positions[jatom][2]:
##
            atoms.positions[[iatom,jatom]]=atoms.positions[[jatom,iatom]]
            atoms.symbols[[iatom,jatom]]  =atoms.symbols[[jatom,iatom]]
#
print("Will Freeze all non-Pd atoms below z =", fix_z)
#
# Report atom1 -> atom2 distances flag big deviations
fix_list=[]
fix_symbol=[]
fix_spec_no=[]
#
print("Atom list for file 1:")
for iatom in range(0,natoms):
    vec=[]
    fix_this=False
#
# Keep track of species counts
    for ispec in range(0,len(species_list)):
       if atoms.symbols[iatom] in species_list[ispec]:
          this_spec=ispec
          species_count[ispec] = species_count[ispec] + 1
#
# To allow constraints make atom2 atoms that will be constrained copies of atom1 list
#
    if (atoms.positions[iatom][2] < fix_z and "Pd" not in atoms.symbols[iatom]):
        fix_this=True 
        fix_list.append(iatom)
        fix_symbol.append(atoms.symbols[iatom])
        fix_spec_no.append(species_count[this_spec])
#
    else:
       if fix_this:
          print(f"%d %s %10.6f %10.6f %10.6f, fixed" \
             % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                atoms.positions[iatom][1], atoms.positions[iatom][2]))
       else:
          print(f"%d %s %10.6f %10.6f %10.6f" \
             % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
print("Fix list: ", fix_list)
#
# Create new cell file
#
fp_out = open("fixed.cell","w")
#
write_cell_djw(fp_out, atoms)
#
# Add constraints to end cell file
#
fp_out.write("\n\n%BLOCK CELL_CONSTRAINTS\n")
fp_out.write("   0    0   0\n   0    0   0\n")
fp_out.write("%ENDBLOCK CELL_CONSTRAINTS\n")
#
fp_out.write("\n\n%BLOCK IONIC_CONSTRAINTS\n")
#
iflist=0
for ifix in range(0,len(fix_list)):
   iflist+=1
   fp_out.write("   %d     %s    %d  1.0000000000000000    0.0000000000000000    0.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
   iflist+=1
   fp_out.write("   %d     %s    %d  0.0000000000000000    1.0000000000000000    0.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
   iflist+=1
   fp_out.write("   %d     %s    %d  0.0000000000000000    0.0000000000000000    1.0000000000000000\n" \
                   %( iflist, fix_symbol[ifix], fix_spec_no[ifix]) )
#
fp_out.write("%ENDBLOCK IONIC_CONSTRAINTS\n")
#
cons = FixAtoms(fix_list)
atoms.set_constraint(cons)
atoms.set_cell(lattice)
write("fixed.in", atoms , format='aims')
#
fp_out.close()
#

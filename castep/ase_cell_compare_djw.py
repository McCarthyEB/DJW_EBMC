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
   fp.write("### Cell file written with write_cell_djw        #####\n")
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
# Find and count the species in an atoms list
#
def find_species(atoms):
#
    species=[]
    nspecs=[]
    num=0
    natoms=len(atoms)

    for iatom in range(0,natoms):
        have_this=False

        for ispec in range(0,num):
            if atoms.symbols[iatom] in species[ispec]:
                have_this=True
                this_spec=ispec

        if have_this:
            nspecs[this_spec]+=1
        else:
            num+=1
            species.append(atoms.symbols[iatom])
            nspecs.append(1)
#
    return species,nspecs
#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
cell_file1 = sys.argv[1]                                           # Directory where we will do the calculations
cell_file2 = sys.argv[2]                                           # Directory where we will do the calculations
#
fp1=open(cell_file1,'r')
fp2=open(cell_file2,'r')
#splstr=cell_file.split('.')
#stem=splstr[0]
#cif_file="%s.cif" % stem
#
print(f"------------------------------------------------------------")
print(f"Comparing cell file          : %s" % cell_file1 )
print(f"with cell file               : %s" % cell_file2 )
print(f"------------------------------------------------------------")
#
atoms1=read_cell_djw(fp1)
atoms2=read_cell_djw(fp2)
#
fp1.close()
fp2.close()
#
nat1 = len(atoms1)
nat2 = len(atoms2)

lattice=atoms1.get_cell()
#
aaa=np.sqrt(np.dot(lattice[0],lattice[0]))
bbb=np.sqrt(np.dot(lattice[1],lattice[1]))
ccc=np.sqrt(np.dot(lattice[2],lattice[2]))
#
species1,nspecs1=find_species(atoms1)
species2,nspecs2=find_species(atoms2)
num_spec1=len(species1)
num_spec2=len(species2)

print("File 1 contains %d atoms" % nat1)
for ispec in range(0,num_spec1):
   print("%s : %d" % ( species1[ispec], nspecs1[ispec]))

print("\n")
#
print("File 2 contains %d atoms" % nat2)
for ispec in range(0,num_spec2):
   print("%s : %d" % ( species2[ispec], nspecs2[ispec]))
# Check species match
#
if num_spec1 != num_spec2:
   print("ERROR: Different number of species in structures, check they are for the same system!")
   exit(0)
#
matched=np.zeros(num_spec1)
for ispec in range(0,num_spec1):
   for jspec in range(0,num_spec2):
       if species1[ispec] in species2[jspec]:
           if nspecs1[ispec] == nspecs2[jspec]:
               matched[ispec]=1
#
if np.sum(matched) != num_spec1:
    print("ERROR: Some species have different numbers in the two structures")
    print("ERROR: Check that they are the same system")
    exit(0)
#
print("Structures match in species and number of each.....")
#
# sort species into alphabetical order
#
for ispec in range(0,num_spec1):
   for jspec in range(ispec+1,num_spec2):
       if species1[ispec] > species2[jspec]:
           spec_tmp=species1[ispec]
           nums_tmp=nspecs1[ispec]
#
           species1[ispec]=species1[jspec]
           nspecs1[ispec]=nspecs1[jspec]
#
           species1[jspec]=spec_tmp
           nspecs1[jspec]=nums_tmp
#
print("Sorted species....")
for ispec in range(0,num_spec1):
   print("%s : %d" % ( species1[ispec], nspecs1[ispec]))
#
# Use species to sort atoms into species blocks
# Note that only species list 1 has been sorted but lists now know to
# be identical so can use this for both structures.
#
atoms1_sorted=Atoms()      
atoms2_sorted=Atoms()      
found1=0 
found2=0 
for ispec in range(0,num_spec1):
    for iatom in range(0,nat1):
        if atoms1.symbols[iatom] in species1[ispec]:
            found1+=1
            atoms1_sorted.append(atoms1[iatom])
#
    for iatom in range(0,nat2):
        if atoms2.symbols[iatom] in species1[ispec]:
            found2+=1
            atoms2_sorted.append(atoms2[iatom])
#
if found1 != nat1:
    print("ERROR: Not all atoms in first structure accounted for when sorting....")
    exit(0)
#
atoms1=atoms1_sorted.copy()       
atoms2=atoms2_sorted.copy()       
#
print(lattice)
print("aaa: %10.6f bbb: %10.6f ccc: %10.6f " % (aaa,bbb,ccc))
#
#write(cif_file, atoms ,format='cif')
#
# reorder atoms according to z-co-ordinate
#
#for iatom in range(0,nat1):
#    for jatom in range(iatom+1,nat1):
#        if atoms1.positions[iatom][2] > atoms1.positions[jatom][2]:
##
#            atoms1.positions[[iatom,jatom]]=atoms1.positions[[jatom,iatom]]
#            atoms1.symbols[[iatom,jatom]]  =atoms1.symbols[[jatom,iatom]]
##
#            atoms2.positions[[iatom,jatom]]=atoms2.positions[[jatom,iatom]]
#            atoms2.symbols[[iatom,jatom]]  =atoms2.symbols[[jatom,iatom]]
#
#
# Report atom1 -> atom2 distances flag big deviations

print("Atom list for file 1:")
for iatom in range(0,nat1):
    vec=[]
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
       print(f"%d %s %10.6f %10.6f %10.6f, dist: %10.6f Angs" \
          % (iatom, atoms1.symbols[iatom], atoms1.positions[iatom][0], \
             atoms1.positions[iatom][1], atoms1.positions[iatom][2], ddd))
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
# Create a trajectory for a linear interpolation
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
       fp_out=open("new_start.cell","w")
       write_cell_djw(fp_out,atoms)

   if istep == 5:
       fp_out=open("new_intermed.cell","w")
       write_cell_djw(fp_out,atoms)

   if istep == 10:
       fp_out=open("new_end.cell","w")
       write_cell_djw(fp_out,atoms)
#

#















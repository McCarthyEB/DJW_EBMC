# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
from scipy import sparse

from ase import Atom
from ase.io import read
from ase import neighborlist

def find_mols(atoms):
#
# Try to get connectivity of atoms
#
    cutoff=neighborlist.natural_cutoffs(atoms)
    neighs=neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
    neighs.update(atoms)
    neigh_matrix=neighs.get_connectivity_matrix()
#
    print(np.shape(neigh_matrix))
#
# Identify molecules as index lists
#
    n_mols, comp_list= sparse.csgraph.connected_components(neigh_matrix)
#
    print("There are %d molecules in the cell" % n_mols)
#
    molinds=[]
    for idx in range(0,n_mols):
         molinds.append([ i for i in range(len(comp_list)) if comp_list[i]==idx ])
#
    return n_mols, molinds, neigh_matrix
#

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
#
# Read in the list of atoms from the structure file into the ase Atoms object.
# the loop prints out symbols and co-ordinates, note how to get these sub-structures
# out of the atom object.
#
    atoms=read("C5H4O2.in")
    natoms=len(atoms)
    for iatom in range(0,natoms):
        print(f"%d %s %10.6f %10.6f %10.6f "                 \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
# The find_mols routine identifies atoms that are chemically bonded to one another based on atom-atom
# distances. These sub-groups of atoms are called "molecules". The routine returns the number of molecules
# found "n_mols" and a two dimensional array, molinds[imol][iatom] gives the index of the iatom atom within
# molecule imol in the atoms list.
# neigh_matrix[i,j] is 1 if atom i is bonded to atom j and zero otherwise
#
    n_mols, molinds, neigh_matrix = find_mols(atoms)
#
# Define atoms that are the bond to be rotated
#
    neigh_lists=[[] for iatom in range(0,natoms)] 
#
    for imol in range(0, n_mols):
        print(f"Atoms in molecule %d" % imol)
        for iatom in range(0,len(molinds[imol])):
            indx=molinds[imol][iatom]
            print(f"%d %s %10.6f %10.6f %10.6f "                 \
                % (indx, atoms.symbols[indx], atoms.positions[indx][0], \
                     atoms.positions[indx][1], atoms.positions[indx][2]))
            for jatom in range(0,len(molinds[imol])):
               jndx=molinds[imol][jatom]
#
# Make a list of lists to hold each atoms chemical neighbours.
#
               if neigh_matrix[indx,jndx] == 1:
                   print(".....bonded to %d %s" % (jndx,atoms.symbols[jndx]))
                   neigh_lists[indx].append(jndx)
#
# dihed_inds defines the bond around which to rotate
# the dihedral angle to be rotated are defined by the elements                
# 0--1--2--3 the fragment containing 2--3 will be moved
#
    dihed_inds=[0,1,5,9]
#
# gather the indices of the two fragments on either side of the bond
#
    frag1=[dihed_inds[1]]
    frag2=[dihed_inds[2]]
#
# Check that the four atoms defining the dihedral are in the same molecule
#
    found=0
    mol_indx=[]
    for imol in range(0, n_mols):
        for iatom in range(0,len(molinds[imol])):
            indx=molinds[imol][iatom]
            for iii in dihed_inds:
               if indx==iii:
                  found=found+1
                  mol_indx.append(imol)
#
    print(found)
    if found == 4:
       print("All dihedral atoms are in molecule %d..." % mol_indx[0])
    else:
       print("ERROR: The four atoms defining the dihedral are not in the same molecule")
       print("ERROR: Molecular indices found as: ", mol_indx)
       exit(0)
#
    ends_frag1=[dihed_inds[1]]
    ends_frag2=[dihed_inds[2],dihed_inds[3]]
#
# Now look for fragments within the molecule referred to by mol_indx[0]                
#
    have_ends = len(ends_frag1) > 0
    found_list=[]
#
    while have_ends:
       indx=ends_frag1[0]
#
# Add neighbours of atom indx to list
#
       print("Latest end atom %d has %d neighbours.." % (indx, len(neigh_matrix[indx])))
#
       for iatom in range(0, len(neigh_lists[indx])):
#
# neigh_matrix is not a list of neighbours its just 1 if bonded and 0 if not.
# need to set up a list of lists for atom neighbours earlier.
#
          ineigh=neigh_lists[indx][iatom]
#
          have_this = False
          for ifnd in found_list:
             if ineigh == ifnd:
                have_this = True
#
          for ifnd in ends_frag1:
             if ineigh == ifnd:
                have_this = True
#
          if not have_this and ineigh != dihed_inds[2]:
             ends_frag1.append(ineigh)
#
       found_list.append(indx)
       del ends_frag1[0]
#
       have_ends = len(ends_frag1) > 0
       print("Ends list: ", ends_frag1)          
#
    print("\nFound frag1 as: ", found_list)
    
          
       
    
# 
#
    
     
    




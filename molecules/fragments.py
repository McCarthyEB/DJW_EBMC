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
    for iatom in range(0,len(atoms)):
        print(f"%d %s %10.6f %10.6f %10.6f "                 \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
# The find_mols routine identifies atoms that are chemically bonded to one another based on atom-atom
# distances. These sub-groups of atoms are called "molecules". The routine returns the number of molecules
# found "n_mols" and a two dimensional array, molinds[imol][iatom] gives the index of the iatom atom within
# molecule imol in the atoms list.
#
    n_mols, molinds, neigh_matrix = find_mols(atoms)
#
# Define atoms that are the bond to be rotated
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
               if neigh_matrix[indx,jndx] == 1:
                   print(".....bonded to %d %s" % (jndx,atoms.symbols[jndx]))
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
       print("All dihedral atoms in the same molecule...")
    else:
       print("ERROR: The four atoms defining the dihedral are not in the same molecule")
       print("ERROR: Molecular indices found as: ", mol_indx)
       exit(0)
#
    ends_frag1=[dihed_inds[0],dihed_inds[1]]
    ends_frag2=[dihed_inds[2],dihed_inds[3]]
    for iatom in range(0,len(molinds[mol_indx[0]])):
       indx=molinds[imol][iatom]
       
    
# 
#
    
     
    




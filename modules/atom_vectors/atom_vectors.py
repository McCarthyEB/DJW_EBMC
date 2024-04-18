import numpy as np
from scipy import sparse
#
from ase import Atom
from ase import Atoms
from ase import neighborlist
#
# Return the distance between an atom and a reference atom list
# if ilist is the word "all" consider all atoms in the reference list.
# otherwise ilist should be a list of indices so only consider atoms in the list. 
#
def atom_dist(ref_atoms, atom2, ilist):
#
   top_atom=len(ref_atoms)
   if ( ilist == "all" ):
       ilist=[ i for i in range(0,top_atom) ]
#
   ref_atoms.append(atom2)
   dist_list = ref_atoms.get_distances(top_atom, ilist, mic=True)
   del ref_atoms[top_atom]
   return dist_list
#
# function to add to a dist_list, these are lists of all the distances between atoms
# in an atom list, used to distinguish configurations.
#
def add_to_dist_list(dist_list, ref_Hatoms, flag_1, flag_2, map_1, map_2, is_self):
#
# Case with list with itself, run second index over numbers greater than
# first to avoid double counting.
#
     if is_self:
        n1_sites=len(flag_1)
#
        for i_1 in range(0, n1_sites):
           if flag_1[i_1] == 1:
#
             ilist=[]
             for jjj in range(i_1+1, n1_sites):
                if flag_1[jjj] == 1:
                   ilist.append(map_1[jjj])

             if len(ilist) > 0:
                 dist_next=atom_dist(ref_Hatoms, ref_Hatoms[map_1[i_1]], ilist)
#
                 for ind in range(0,len(dist_next)):
                          dist_list.append(dist_next[ind])
#
#                print("CODE CHECK:")
#                print("For ind1 = ", i_1," ilist: ", ilist)
#                dum=[round(i,2) for i in dist_list]
#                print("dist_list now: ", dum)
#
# Case with two separate lists do list 1 with all of list 2
#
     else: 
        n1_sites=len(flag_1)
        n2_sites=len(flag_2)
#
        for i_1 in range(0, n1_sites):
           if flag_1[i_1] == 1:
#
             ilist=[]
             for jjj in range(0, n2_sites):
                if flag_2[jjj] == 1:
                   ilist.append(map_2[jjj])

             if len(ilist) > 0:
                 dist_next=atom_dist(ref_Hatoms, ref_Hatoms[map_1[i_1]], ilist)
#
                 for ind in range(0,len(dist_next)):
                       dist_list.append(dist_next[ind])
#
#                 print("CODE CHECK:")
#                 print("For ind1 = ", i_1,"with ind2 ilist: ", ilist)
#                 dum=[round(i,2) for i in dist_list]
#                 print("dist_list now: ", dum)
#
     return
#
# find_mols : works out the connectivity of an atom list and looks for closed sets to define molecules
#
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
# Define the bonded neighbours of each atom in a list from a neigh_matrix
#
def gen_neigh_list(atoms, neigh_matrix):
#
    natoms=len(atoms)
#
    neigh_lists=[[] for iatom in range(0,natoms)]
#
    for iatom in range(0,natoms):
        for jatom in range(iatom+1,natoms):
#
# Make a list of lists to hold each atoms chemical neighbours.
#
               if neigh_matrix[iatom,jatom] == 1:
                   neigh_lists[iatom].append(jatom)
                   neigh_lists[jatom].append(iatom)
#
    return neigh_lists
#
# mol_frag: returns the list of indices for the fragment of molecule molecule imol that has iatm1 in it but not iatm2
#           neigh_lists needs to be a list of lists for the atom neighbours
#
def mol_frag(atoms, imol, iatm1, iatm2, molinds, neigh_lists):
#
# Check that the atoms are in the molecule indicated.                     
#
    atm1_in_mol=False
    atm2_in_mol=False
#
    for iatom in molinds[imol]:
       if iatm1 == iatom:
          atm1_in_mol=True
       if iatm2 == iatom:
          atm2_in_mol=True
#
    if not atm1_in_mol:
       print("ERROR: In mol_frag function, atom %d is not in molecule %d as call suggests..." % \
                                      ( iatm1, imol ))
#
    if not atm2_in_mol:
       print("ERROR: In mol_frag function, atom %d is not in molecule %d as call suggests..." % \
                                      ( iatm2, imol ))
#
    have_ends = True
    found_list=[]
    ends=[iatm1]
#
# Add neighbours of atom indx to list, keep going till tree runs out of ends
#
    while have_ends:
       indx=ends[0]
#
#       print("Latest end atom %d has %d neighbours.." % (indx, len(neigh_lists[indx])))
#
       for ineigh in neigh_lists[indx]:
#
          have_this = ineigh in found_list or ineigh in ends
#
          if not have_this and ineigh != iatm2:
             ends.append(ineigh)
#
       found_list.append(indx)
       del ends[0]
#
       have_ends = len(ends) > 0
#       print("Ends list: ", ends)
#
#    print("\nFound fragment as: ", found_list)
#
    return found_list




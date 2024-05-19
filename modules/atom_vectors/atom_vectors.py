#
# Routines to carry out vector operations for analysis of atoms lists
#
# atom_dist : returns the distances from an atom to a list of other atoms.
#             if a list of indices are supplied only those members of the atom list 
#             are included. If the word "all" is sent in the place of the list all
#             atoms in the list are included.
# add_dist_to_list : "dist_lists" are lists of interatomic distances. This function 
#                    adds to the list.
# find_mols : Looks for chemically bonded atoms in an atom list and groups closed sets 
#             into molecules. Each atom in the list is given a molecule index, the atom
#             order in the list is not affected.
# gen_neigh_list : generates a list of lists for a set of atoms so that each member of the
#                  list gives the indices of the atoms bonded to the corresponding atom.
# mol_frag: returns the list of indices for the fragment of molecule imol that has 
#           iatm1 in it but not iatm2 neigh_lists needs to be a list of lists for 
#           the atom neighbours
# 
#
import numpy as np
from scipy import sparse
#
from ase import Atom
from ase import Atoms
from ase import neighborlist
from ase.constraints import FixAtoms
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
# mol_frag: returns the list of indices for the fragment of molecule imol that has iatm1 in it but not iatm2
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
#
# map atoms to sites
# routine takes a set of atoms and a set of sites ( also atom types ) and
# returns a list of the indices for the nearest lattice site for each atom
# sites should be a list of lists for allowing site types to be defined. 
#
def map_atoms_to_sites(atoms, sites):
#
   natoms=len(atoms)
   nsite_types=len(sites)
#
   print("natoms = %d" % natoms)
   print("nsite_types = %d" % nsite_types)
#
   nsites=0
   for isite in range(0,nsite_types):
      for jsite in range(0,len(sites[isite])):
        nsites= nsites + 1
#
   if len(atoms) > nsites:
      print("ERROR: In map_atoms_to_sites....")
      print("ERROR: Sent more atoms than sites to match")
      print("ERROR: natoms=%d while nsites=%d." % ( natoms, nsites))
      exit(0)
#
   type_list = []
   match_list = []
   min_dist = []
#
   for iatom in range(0,natoms):
      dist_list = []
      min_dist_list = []
      temp_list = []
#
      for isite in range(0,nsite_types):
#
# work out the best match for this site type and record
#
         dist_list.append(atom_dist(sites[isite], atoms[iatom], "all"))
         temp_list.append(np.argmin(dist_list[isite]))
         min_dist_list.append(dist_list[isite][temp_list[isite]])
#
#      print("dist_list    : ", dist_list)
#      print("temp_list    : ", temp_list)
#      print("min_dist_list: ", min_dist_list)
#      exit(0)
#
# Now work out best match over all
#
      type_list.append(np.argmin(min_dist_list))
      match_list.append(temp_list[type_list[-1]])
      min_dist.append(dist_list[type_list[-1]][match_list[-1]])
#
   return type_list, match_list, min_dist
#
# find the closest centre of mass for "atoms" to the underlying "top_atoms", i.e. the
# cofm that is closest to any of the atoms. There are multiple centres of mass as in a
# periodic system there are many ways to take the set of minimum image vectors.
#
# Move the H_bri list to minimum image positions with respect to each H atom in turn
# Choose cofm closest to any of the top layer atoms for consistency between structures
#
def min_cofm(atoms, top_atoms, verbose):
#
     natoms=len(atoms)
     ntop_atoms=len(top_atoms)
     top_list=[ i for i in range(0,ntop_atoms)]
     if verbose:
        print("min_cofm for : ", natoms," atoms")
        print("top atoms num: ", ntop_atoms)
     cofm_vecs=[]
     dist2_cofm=[]
#
# Use the iorig atom as the one to bring others to on minimum image
#
     for iorig in range(0,natoms):
#
        temp_atoms=atoms.copy()
#
        for iatom in range(0,natoms):
           vec=[0.0,0.0,0.0]
           if iatom != iorig:
               vec=temp_atoms.get_distances(iorig,iatom,mic=True,vector=True)[0]
#
           temp_atoms.positions[iatom]=temp_atoms.positions[iorig]+vec
#
        cofm_vecs.append(temp_atoms.get_center_of_mass())
#
# Find cofm closest to any top_atoms and return that
#
     dist_list=[]
     for ic in range(0,len(cofm_vecs)):
        new_atom=Atom("Cl",cofm_vecs[ic])
        top_atoms.append(new_atom)

        ddd=top_atoms.get_distances(ntop_atoms,top_list, mic=True)
        dist_list.append(np.min(ddd))
#
        del top_atoms[ntop_atoms]
#
     imin=np.argmin(dist_list)
#
     if verbose:
        print(len(cofm_vecs)," cofm vectors found:")
        for iii in range(0,len(dist_list)):
           print(iii," : ",cofm_vecs[iii], " dist2: ", dist_list[iii])
        print("Selected case ",imin)
#
     return dist_list[imin]
#
# Bubble sort atom list along a supplied direction.
# Bubble sort for atom list, this will sort atoms according to their
# co-ordinates in x,y,z according to the dir vector, 
# e.g. dir= (0,0,1) for sorting along the z-direction.
#
def bubble_atoms(atoms, dir):
#
   natoms=len(atoms)
   if natoms < 0:
      print("ERROR: Cannot bubble sort an empty atom list")
      exit(0)
#
#   print("In bubble_atoms with ", natoms, " atoms. Sorting in direction: ", dir)
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
##
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
# z-fix function to fix atoms below a certain z-co-ordinate in slab calculations
# if a pre-existing list of fixed indices is passed this is appended to.
#
def zfix(atoms, zfix_string, zfix_list=[-1]):
#
   fix_map=np.zeros(len(atoms))
#
   test = zfix_string.replace('.','',1).isdigit()
   if test:
      zfix = float(zfix_string)
      need_zfix = True
      print("Will fix atoms with z co-ordinate less that %10.6f" % zfix)
   elif zfix_list[0] < 0:
      print("zfix and no list passed so so will not freeze any atom co-ordinates")
      need_zfix = False
   else:
      print("zfix not set but specific list of atom indices passed so will proceed")      
      need_zfix=False
#
# preserve any constraints already set
# Assume these include lower layers if fixed so only implement zfix if POSCAR all free
   have_cons=False
   cons_orig = (atoms.constraints).copy()
#
# make list according to zfix if set
# if no pre-existing list set the zfix_list to empty
   if zfix_list[0] < 0:
     zfix_list=[]
#
   if need_zfix:
     for iatom in range(0,len(atoms)):
        if (atoms.positions[iatom][2] < zfix):
            zfix_list.append(iatom)
            have_cons=True
#
# Remove duplicates ( i.e. atoms that appear in passed list and have z-co-ordinate in range )
#
   remove=[]
   for ifix in range(0,len(zfix_list)):
     for jfix in range(ifix+1,len(zfix_list)):
        if zfix_list[ifix] == zfix_list[jfix]:
            remove.append(zfix_list[jfix])
#
   for iii in range(0,len(remove)):
      zfix_list.remove(remove[iii])
# 
   print("zfix_list :"),
   print(zfix_list)
#
   have_cons_orig=False
#
   if (cons_orig):
#
# Find the indices set in the POSCAR file as fixed
#
     have_cons_orig=True
     sub1=str(cons_orig[0]).split('[')
     sub2=sub1[1].split(']')
     sub1=sub2[0].split(',')
     fixed_list=[ int(sub1[i]) for i in range(0,len(sub1)) ]
     print("Original fixed_list:" )
     print(fixed_list)
#
# Combine lists if you have both, if only zfix copy that to fixed_list
#
     for ilist in range(0, len(zfix_list)):
       found_fix=False
       for jlist in range(0,len(fixed_list)):
          if ( zfix_list[ilist] == fixed_list[jlist] ):
             found_fix=True
       if ( not found_fix ):
          fixed_list.append(zfix_list[ilist])

     fixed_list.sort()
     print(f"Combined list: ")
     print(fixed_list)
#
   elif ( len(zfix_list) > 0 ):
     fixed_list=zfix_list.copy()
#
   del atoms.constraints
   cons_orig=FixAtoms(indices=fixed_list)
   for icon in fixed_list:
      fix_map[icon]=1
   cons = FixAtoms(fixed_list)
   atoms.set_constraint(cons)

   return fix_map
#
#
# find the closest centre of mass for "atoms" to the underlying "top_atoms", i.e. the
# cofm that is closest to any of the atoms. There are multiple centres of mass as in a
# periodic system there are many ways to take the set of minimum image vectors.
#
# Move the H_bri list to minimum image positions with respect to each H atom in turn
# Choose cofm closest to any of the top layer atoms for consistency between structures
#
def min_cofm2(atoms, top_atoms, supercell, verbose=False):
#
     natoms=len(atoms)
#
# return dummy data if there are no atoms in the list
#
     if natoms == 0:
        cofm_info = [-1, -1, -1, -1, -1 ]
        return cofm_info
#
     atom_list=[ i for i in range(0,natoms)]
     ntop_atoms=len(top_atoms)
     top_list=[ i for i in range(0,ntop_atoms)]
#
     if verbose:
        print("min_cofm for : ", natoms," atoms")
        print("top atoms num: ", ntop_atoms)
        print("super cell   : ", supercell)
     cofm_vecs=[]
     dist2_cofm=[]
#
# Use supercell setting to work out possible origins
#
     cell=top_atoms.get_cell().copy()
#
# Note this is set for slabs so expect supercell to be 1 in c direction
     if supercell[2] != 1:
        print("ERROR: min_cofm origin settings expecting slab with c perpendicular to surface")
        exit(0)
#
     finc=[1.0/supercell[0], 1.0/supercell[1], 1.0/supercell[2]] 
     origins=[]
     for iorig in range(0,supercell[0]):
        for jorig in range(0,supercell[1]):
           for korig in range(0,supercell[2]):
              orig= iorig*finc[0]*cell[0] +  jorig*finc[1]*cell[1] \
                          +  korig*finc[2]*cell[2] + [0,0,top_atoms.positions[0][2]]
              origins.append(orig)
#
# Loop over origins, identify most compact structure around cofm
#
     size_min=10000.0
     for iorig in range(0,len(origins)):
#
#       if verbose:
#          print("Origin: ", origins[iorig])
#
# Add an atom as an origin marker
# Use the origin atom as the one to bring others to on minimum image
#
        orig_atoms=atoms.copy()
        scal_atoms=atoms.copy()
        orig_atoms.append(Atom("C",origins[iorig]))
#
        vecs=orig_atoms.get_distances(natoms,atom_list,mic=True,vector=True)
        del orig_atoms[natoms]
#
# Move atoms to minimum image with respect to this origin
#
        for ivec in range(0,natoms):
#         print(ivec, " : ", vecs[ivec])
          scal_atoms.positions[ivec]=vecs[ivec]
          orig_atoms.positions[ivec]=origins[iorig]+vecs[ivec]
#
# On first pass of origin loop identify atoms that have 0.5 fractional co-ordinates
# as we have to try two positions for each of these to be sure to build the minimum
# cluster.
#
        scaled=scal_atoms.get_scaled_positions()
#       print(scaled)
#
# The spec_list will contain lists of atom id on special position, which co-ord, and
# displacement direction to check.
#
        nspecial=0
        spec_list=[]
        for iscal in range(0, len(scaled)):
          for icoord in range(0,2):
             sss = scaled[iscal][icoord] 
             if abs(sss) < 0.501 and abs(sss) > 0.499:
#
# Only bother if atom is aligned with origin
               if icoord == 0:
                  algn = abs(orig_atoms.positions[iscal][1] - origins[iorig][1]) < 0.1
               elif icoord == 1:
                  algn = abs(orig_atoms.positions[iscal][0] - origins[iorig][0]) < 0.1
#
               if algn:
                   nspecial=nspecial + 1
                   dv = np.dot(vecs[iscal],vecs[iscal])
                   vec_test = [vecs[iscal] + idisp*cell[icoord] for idisp in range(-1,2,2)]
                   trans_ddd = [ np.dot(vec_test[id],vec_test[id]) for id in range(0,2)]
                   min_ddd = np.min(trans_ddd)
#
                   if abs(trans_ddd[0] - dv ) < 0.1:
                      spec_list.append([iscal,icoord, -1])
                   elif abs(trans_ddd[1] - dv ) < 0.1:
                      spec_list.append([iscal,icoord,  1])
                   else:
                      print("ERROR: Atom on special point but no equivalent translation found..")
                      exit(0)
#
#                  print("atom %d coord %d is special" % ( iscal, icoord))
#                  print("trans_ddd: ", trans_ddd)
                
#       print("Found %d special point positions" % nspecial)
#       print(spec_list)
#
        spec_combs = [[[-1,-1,-1]]]
#
        if nspecial == 1:
          spec_combs.append([spec_list[0]])
        elif nspecial == 2:
          spec_combs.append([spec_list[0], [-1,-1,-1]])
          spec_combs.append([spec_list[1], [-1,-1,-1]])
          spec_combs.append([spec_list[0], spec_list[1]])
      
        elif nspecial > 2:
          print("ERROR: Cann't cope with more than 2 atoms on special positions...")
          exit(0)

#
#       print("spec_combs:")
#       for ispec in range(0,len(spec_combs)):
#          print(spec_combs[ispec])
#       print("--------------------------------")
#
# size around origin, vecs are still valid for orig_size as distance to origin of all
# images is the same.
#
        orig_size=np.sum([np.dot(vecs[iii],vecs[iii]) for iii in range(0,len(vecs))])
#
# Deal with special positioned atoms, for no move the first element of spec_combs will be -1
#
        for ispec in range(0,len(spec_combs)):
#
#          print("Pass ispec = %d......" % ispec)
           temp_atoms=orig_atoms.copy()
#
# Make any moves required, first entry is always no change
#
           if ispec > 0:
              this_spec=spec_combs[ispec]
              for iii in range(0,nspecial):
                 indx   = this_spec[iii][0] 
# see if anything to do for this atom
                 if indx >= 0:
                    icoord = this_spec[iii][1] 
                    sgn    = this_spec[iii][2] 
#                   print("shifting atom position %d coord %d in direction %d" % ( indx, icoord, sgn))
                    temp_atoms.positions[indx] = temp_atoms.positions[indx] + sgn * cell[icoord]
#
# find cofm for this cluster
#
           cofm=temp_atoms.get_center_of_mass()
#
# Cluster "size" around cofm 
#
           temp_atoms.append(Atom("C",cofm))
#
           vecs=temp_atoms.get_distances(natoms,atom_list,mic=False,vector=True)
           del temp_atoms[natoms]
#
           cofm_size=np.sum([np.dot(vecs[iii],vecs[iii]) for iii in range(0,len(vecs))])
#          print("Now got cofm_size %8.4f" % cofm_size)
#
# Work out longest inter-atomic distance in this setting of the cluster
#
           if natoms == 1:
             dmax = 0
             ddd_inter = [dmax]
           elif natoms == 2:
             dmax = temp_atoms.get_distances(0,1,mic=False)
             ddd_inter = [dmax]
           else:
             ddd_inter=[]
             for iatom in range(0,natoms-1):
                for jatom in range(iatom+1,natoms):
                  ddd_inter.append(temp_atoms.get_distance(iatom,jatom,mic=False))
#            print("ddd_inter: ", ddd_inter)
             dmax = np.max(ddd_inter)   
#
# Use all measures together for "size"
#
           size = orig_size + cofm_size + dmax
#          if verbose:
#             print("orig_size=%8.4f, cofm_size=%8.4f, dmax=%8.4f, tot=%8.4f" % \
#                             (orig_size , cofm_size , dmax, size))
#
           if size < size_min:
#
               size_min      = size     
               orig_size_min = orig_size
               cofm_size_min = cofm_size
               dmax_min      = dmax
#
               ddd_inter_min = ddd_inter.copy()
               origin_min = origins[iorig]
#
# record cofm for this cluster
#
               cofm_min=cofm

               temp_atoms_min=temp_atoms.copy()
# Find closest top_atoms to cofm              
#
     top_atoms.append(Atom("C",cofm_min))

     ddd=top_atoms.get_distances(ntop_atoms,top_list, mic=True)
     dist_min=np.min(ddd)
#
     del top_atoms[ntop_atoms]
#
     cofm_info = [dist_min, orig_size_min, cofm_size_min, dmax_min, size_min]
#
     if verbose:
         print("Most compact cluster with F for cofm")
         for iatom in range(0,natoms):
            print("atom %10.6f  %10.6f  %10.6f %s" % \
                    (temp_atoms_min.positions[iatom][0],temp_atoms_min.positions[iatom][1],\
                     temp_atoms_min.positions[iatom][2],temp_atoms_min.symbols[iatom]))  

         print("atom %10.6f  %10.6f  %10.6f F" % (cofm_min[0],cofm_min[1],cofm_min[2]))  
         print("atom %10.6f  %10.6f  %10.6f C" % (origin_min[0],origin_min[1],origin_min[2]))  
#
         print("cofm_info: dist = %8.4f orig size = %8.4f cofm size = %8.4f dmax = %8.4f tot_measure = %8.4f" % \
                    ( cofm_info[0], cofm_info[1], cofm_info[2], cofm_info[3], cofm_info[4] ))
         print("ddd_inter in cofm: ", ddd_inter_min)
#
# return the distance to the nearest top atom and the compactness measure
#
     return cofm_info
#
# use_list makes a test that the smallest member of a list is within a desired range.
#
def use_list(dists, low_min, high_min=-1):
#
   bother = True
#
   if len(dists) > 0:
     min = np.min(dists)
#
# Check if the minimum in the list is below the low value threshold
#
     if  min < low_min:
       bother = False
#
# Only test for high_min if a positive value has been set
#
     elif high_min > 0:
       if min > high_min:
          bother = False

   return bother
#



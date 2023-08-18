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
# Function to obtain HOME
from os.path import expanduser
#
from ase import Atom
from ase import Atoms
from ase import geometry

from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength

from ase.io import read, write, Trajectory

from ase.calculators.vasp import Vasp2

home = expanduser("~")
set_module="%s/python/modules/slab" % home
sys.path.append(set_module)
set_module="%s/python/modules/atom_settings" % home
sys.path.append(set_module)
set_module="%s/python/modules/num_digi" % home
sys.path.append(set_module)
set_module="%s/python/modules/vectors" % home
sys.path.append(set_module)
#
print("Current system path:")
print(sys.path)
#
# Import our own modules
#
from set_atoms import atom_formal_charges
#
from binary_tools import *
from vectors import *
#
# Local definition of atom...atom distance measurement
#
def atom_dist(atoms, i1, i2):
#
#  Inter-atomic vector
#
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])
#
#  impose minimum image convention
#
    latt      = (atoms.get_cell()).copy()
    recip_latt= (atoms.get_reciprocal_cell()).copy()
#
    print(latt)
    print(recip_latt)
#
    vec= np.array([ atoms.positions[i2,0]-atoms.positions[i1,0],    \
                    atoms.positions[i2,1]-atoms.positions[i1,1],    \
                    atoms.positions[i2,2]-atoms.positions[i1,2] ])

    size = np.sqrt(np.vdot(vec,vec))
    vec = vec/size
    
    return size
#
# Function to find and record layers within a structure assuming
# z-axis is perpendicular to surface
#
def find_layers(atoms, tol):

   z_this=atoms.positions[0][2]
   layers=[]
   flag=np.ones(len(atoms))
   first=True
#
   while np.sum(flag) != 0:
      this_layer=[]
      if first:
        first=False
        this_layer.append(0)
        flag[0]=0
#
      for iatom in range(1,len(atoms)):
         if (flag[iatom]==1):
            if (   (atoms.positions[iatom][2] > (z_this - tol)) \
                 & (atoms.positions[iatom][2] < (z_this + tol)) ):
               flag[iatom]=0
               this_layer.append(iatom)
        
      layers.append(this_layer)
# find next z
      for iatom in range(1,len(atoms)):
         if (flag[iatom]==1):
            z_this=atoms.positions[iatom][2]
            break
#
   return layers      
#
# Function to add a methyl group along z-direction from C markers.
#
def add_methyl(atoms, top_layer, ads_atom):
#
# tetrahedral angle in radians
# sixty degrees in radians
#
   rad_60=np.radians(60.0)
   rad_19p5=np.radians(19.5)
   sin19p5=np.sin(rad_19p5)
   cos19p5=np.cos(rad_19p5)
   CH_ideal=1.0905
#
   natoms=len(atoms)
#
   for iatom in range(0,natoms):
      if atoms.symbols[iatom] == ads_atom:
         coords=atoms.positions[iatom].copy()
#
# Place second carbon at the C-C bond length from test optimisations
# ads_slab_2_0_17: for HCP and ads_slab_2_17_0 : for FCC
# C-C average lengths for hcp or fcc same to 4 d.p.
#
         coords[2]=coords[2]+1.4896
         atoms.append(Atom("C", coords))
#
# Place hydrogens at tetrahedral angle by making 109.5 (90+19.5) angle
# hydrogens form an equilateral triangle staggered with respect to Pd surface atoms.
# hold vector components relative to the top C atom ( at coords ).
         znew=CH_ideal*sin19p5
#
         dists=atoms.get_distances(iatom, top_layer, mic=True)
#
         vec=[]
         print("Finding neighbours to atom %d..." % iatom)
         for idist in range(0,len(dists)):
           if dists[idist] < 2.0:
              print("Atom %d fits the bill" % top_layer[idist])
              vec.append(atoms.get_distance(iatom, top_layer[idist], mic=True, vector=True))
              print(dists[idist])
              print(vec[-1])
#
         if len(vec) != 3:
            print("ERROR: Adsorbate atom %d does not have three surface metal neighbours" % iatom)
            exit(0)
#
# coords is still the co-ordinates of the second C atom
#
         xy=[]
         xy.append([(vec[0][0]+vec[1][0])/2, (vec[0][1]+vec[1][1])/2])
         xy.append([(vec[1][0]+vec[2][0])/2, (vec[1][1]+vec[2][1])/2])
         xy.append([(vec[2][0]+vec[0][0])/2, (vec[2][1]+vec[0][1])/2])
#
# normalise the xy vector and then scale and add in the top C atom coords
#
         for iii in range(0,3):
           scal=CH_ideal*cos19p5/np.linalg.norm(xy[iii])
           xy[iii]= np.multiply(scal,xy[iii])
#
           new_coords=coords+[xy[iii][0],xy[iii][1],znew]
           atoms.append(Atom("H", new_coords))

   return
    
#
# Main code begins
#
# debug: setting "True" allows checking of code without calling fhiaims
#        this can be run in the foreground to check set up before 
#        commiting to a fhiaims calculation. Useful here for magnetic
#        ordering check.
#        setting "False" will run the fhiaims jobs so only do as part of a
#        job in a queue
#
debug = True 
#
if (debug):
  print(f"DEBUG version, checking code, fhiaims will not be run")
  print(f"DEBUG foreground execution with additional printing")
#
  workdir = os.getcwd()                                           # Directory where we will do the calculations
  geom_file=sys.argv[1]                                           # structure file will be passed to script
  cores_per_task = "not defined in debug mode"                    # cores assigned to tasks where available
  jobid = "not defined in debug mode"                             # Unique identifier for this job
  stdout_file= "not defined in debug mode"                        # stdout file will appear in your launch directory
  results_dir= os.getcwd()                                        # directory for results, this should be on your $HOME space
  traj_file = "not defined in debug mode"                         # File for trajectory output
  traj_spin_file = "not defined in debug mode"                    # File for trajectory output
  machine = "not defined in debug mode"                           # string defining machine
else:
  print(f"Arguements passed {len(sys.argv)}")
  for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
# this ase_vasp_opt.py script expects:
# the first arguement to be the directory where we will run vasp, 
# second the number of cores available or a junk string if this is not required
# The third to be the JOB_ID, a unique identifier to add to output information
# The fourth is the path to the directory from which the job was launched
# The fifth is the sub-directory for the structure for this particular run.
# The sixth is the sub-directory for the sub_directory for this type of job, under the structure
# The seventh to be the name of the machine we are running on.
#           Currently this can be "hawk","thomas" or  "archer".
#
#
  workdir = sys.argv[1]                                           # Directory where we will do the calculations
  cores_per_task = sys.argv[2]                                    # cores assigned to tasks where available
  jobid = sys.argv[3]                                             # Unique identifier for this job
  stdout_file= "%s/vasp_%s_%s.stdout" % (sys.argv[4], sys.argv[5], jobid) # stdout file for fhiaims will appear in your launch directory
  results_dir= "%s/%s/%s" % (sys.argv[4], sys.argv[5], sys.argv[6]) # directory for results, this should be on your $HOME space
  traj_file = "%s/%s_%s_spin1.traj" % ( results_dir, sys.argv[5], jobid )      # File for trajectory output
  traj_spin_file = "%s/%s_%s_spin2.traj" % ( results_dir, sys.argv[5], jobid ) # File for trajectory output
  machine = sys.argv[7].lower()                                     # string defining machine
#
print(f"------------------------------------------------------------")
print(f"Work directory               : %s" % workdir )
print(f"Geometry file                : %s" % geom_file )
print(f"stdout file for fhiaims      : %s" % stdout_file )
print(f"results_dir for this script  : %s" % results_dir )
print(f"trajectory file for this run : %s" % traj_file )
print(f"trajectory file for ISPIN 2  : %s" % traj_spin_file )
print(f"Machine running job          : %s" % machine )
print(f"------------------------------------------------------------")
#
# Read in the structure from a cif file   
#
atoms=read(geom_file)
#write('check_start.cif', atoms ,format='cif')
#write('check_start.xyz', atoms ,format='xyz')
#
# Define element symbol of the adsorbate marker
#
ads_elem="C"
#
# Try find_layers
#
layers=find_layers(atoms,0.1)
for ilayer in range(0,len(layers)):
   print("layer %d" % ilayer )
   for iatom in range(0,len(layers[ilayer])):
      indx=layers[ilayer][iatom]
      print("%d %s %10.6f %10.6f %10.6f " % ( indx, atoms.symbols[indx], \
             atoms.positions[indx][0], atoms.positions[indx][1],         \
             atoms.positions[indx][2] ))
#
# In this example we have read in a 5-layer slab with H as layer 6 
# So H atom indices can be copied over from layers[5] ( zero indexed )
#
top_layer=[layers[5][i] for i in range(0,len(layers[5]))]
top_metal_layer=[layers[4][i] for i in range(0,len(layers[4]))]
#
# ----------------------------------------------------------------------
# Divide top_layer atoms into a hcp and an fcc set
# ----------------------------------------------------------------------
#
# Use layer (2) to define the hcp sites
#
top_hcp=[]
for iatom in range(0,len(top_layer)):
   matched=False
   indx=top_layer[iatom]
   for jatom in range(0,len(layers[2])):
      jndx=layers[2][jatom]
      dx=atoms.positions[indx][0]-atoms.positions[jndx][0]
      dy=atoms.positions[indx][1]-atoms.positions[jndx][1]
      vec=np.array([dx,dy])
      if (np.linalg.norm(vec) < 0.1):
         matched=True
   if (matched):
      top_hcp.append(indx)
#
# Use layer (3) to define the fcc sites
#
top_fcc=[]
for iatom in range(0,len(top_layer)):
   matched=False
   indx=top_layer[iatom]
   for jatom in range(0,len(layers[3])):
      jndx=layers[3][jatom]
      dx=atoms.positions[indx][0]-atoms.positions[jndx][0]
      dy=atoms.positions[indx][1]-atoms.positions[jndx][1]
      vec=np.array([dx,dy])
      if (np.linalg.norm(vec) < 0.1):
         matched=True
   if (matched):
      top_fcc.append(indx)
#
print("hcp sites:")
for iatom in range(0,len(top_hcp)):
    indx=top_hcp[iatom]
    print("%d %s %10.6f %10.6f %10.6f " % ( indx, atoms.symbols[indx], \
             atoms.positions[indx][0], atoms.positions[indx][1],         \
             atoms.positions[indx][2] ))
#
print("fcc sites:")
for iatom in range(0,len(top_fcc)):
    indx=top_fcc[iatom]
    print("%d %s %10.6f %10.6f %10.6f " % ( indx, atoms.symbols[indx], \
             atoms.positions[indx][0], atoms.positions[indx][1],         \
             atoms.positions[indx][2]))
#
# Top layer now has two distinct sites defined by top_hcp and top_fcc index arrays
# Check all assigned
#
num_top_hcp=len(top_hcp)
num_top_fcc=len(top_fcc)
num_top_tot=len(top_layer)
#
got_all= num_top_tot-num_top_hcp-num_top_fcc == 0
#
if (got_all):
  print("All top layer atoms assigned as hcp or fcc.....")
else:
  sys.exit("ERROR: Some top layer atoms not assigned to sites, check input")
#
# Some structures to check: clean surface, ads_hcp, ads_fcc
#
clean_surf=atoms.copy()
del clean_surf[top_layer]
#
ads_hcp=atoms.copy()
del ads_hcp[top_fcc]
ads_only_hcp=ads_hcp[[atom.index for atom in ads_hcp if atom.symbol == ads_elem]]
#
ads_fcc=atoms.copy()
del ads_fcc[top_hcp]
#
hcp_file="%s_hcp.cif" % ads_elem
fcc_file="%s_hcp.cif" % ads_elem
ads_file="%s_ads_only.cif" % ads_elem
#
#write('clean_surf.cif', clean_surf, format='cif')
#write(hcp_file, ads_hcp, format='cif')
#write(ads_file, ads_only_hcp, format='cif')
#write(fcc_file, ads_fcc, format='cif')
#
#sys.exit("Checking 5 layer setp up")
#
# work out possible patterns, in current test have 18 sites in all so 9 removals will
# give large number of configurations
#
for num_to_set in range(2,10):
#
# Need to deal with the hcp and fcc sites as distinct. So
# loop through the possible division of adsorbates between the sites.
# loop over possible hcp populations and then set fcc as the differnce
# of these and the total required.
# 
   nconfigs=0
   nunique_tot=0
   unique_list_tot=[]
   degen_unique=[]
   first_pass_outer=True
#
# Deal with cases where num_to_set is greater than the total sites on the 
# hcp sub-lattice
#
   if (num_to_set > num_top_tot):
      sys.exit("ERROR: Request to set more atoms in place than total sites available")
   elif (num_to_set > num_top_hcp):
      max_set_hcp=num_top_hcp
      min_set_hcp=num_to_set-num_top_fcc
   else:
      max_set_hcp=num_to_set
      min_set_hcp=0
#
   print("Indexes for setting %d atoms on %d sites" % (num_to_set, num_top_tot))
   print("There are %d sites of type 1 and %d of type 2" % (num_top_hcp,num_top_fcc))
   print("Loop limits for set 1 : %d to %d" % (min_set_hcp, max_set_hcp))
#
   for num_to_set_hcp in range(min_set_hcp,max_set_hcp+1):
#
# rest must be fcc site atoms
#
      num_to_set_fcc=num_to_set-num_to_set_hcp
#
      print(" ")
      print(" ")
      print("Will set %d of 18 with %d hcp and %d fcc" % (num_to_set, num_to_set_hcp, num_to_set_fcc))
#
# The arrangements can be taken separately for hcp and fcc sites. But for each option 
# on the hcp lattice will have to test all fcc possibilities. Hence total configurations
# is the product of the two.
#
      flag_list_hcp= k_bits_on(num_to_set_hcp,num_top_hcp)
      flag_list_fcc= k_bits_on(num_to_set_fcc,num_top_fcc)
      nlist_hcp=len(flag_list_hcp)
      nlist_fcc=len(flag_list_fcc)
      nlist_tot=nlist_hcp*nlist_fcc
      nconfigs+=nlist_tot
#
      print("This will require %d hcp and %d fcc to be considered, altogether %d configurations" %
                                                            (nlist_hcp, nlist_fcc, nlist_tot))
#
# Now need to split the dists list into those within a set and those between sets
#
      all_dists_hcp=[]
      all_dists_fcc=[]
      all_dists_htf=[]
#
      print
      print("Setting up ", num_to_set, "adsorbates, need to test ", nlist_tot, " arrangements.")
      print
#
# Work out number to expect in dist lists, this is the number of pair distances for that
# arrangement. Again the dist lists need the distances within each set and the distances between
# sets. Within sets we need to avoid double counting, but between sets we do not.
#
# The dist_lists contain the indicies for each atom pair, so for the intra-set lists
# the first entry is for atom 0 and so is a list of all others, the second is for 
# atom 1 and ignores 0 and 1 etc.
# For the inter-set lists we need all the indicies of the other list in each case.
#
      num_in_hcp_dist_list=0
      num_in_fcc_dist_list=0
      num_in_htf_dist_list=0
#
      dist_list_hcp=[]
      dist_list_fcc=[]
      dist_list_htf=[]
#
      unique_hcp_list=[]
      unique_fcc_list=[]
      unique_htf_list=[]
      num_of_unique=[]
#
      for ilist in range(0,num_to_set_hcp):
         num_in_hcp_dist_list += ilist
         if (ilist < num_to_set_hcp-1):
           dist_list_hcp.append([j for j in range(ilist+1,num_to_set_hcp)])
#
# Note that in later code all the H atoms that are to be added will be
# placed in a single list and this will be done hcp followed by fcc so that
# the fcc site indices follow the hcp.
#
      for ilist in range(0,num_to_set_fcc):
         num_in_fcc_dist_list += ilist
         if (ilist < num_to_set_fcc-1):
           dist_list_fcc.append([j+num_to_set_hcp for j in range(ilist+1,num_to_set_fcc)])
#
# No double counting in the hcp to fcc list
      num_in_htf_dist_list=num_to_set_hcp*num_to_set_fcc
#
# Now we need a list of all the fcc sites for each hcp site, so that all inter-adsorbate
# distances will be obtained. Note that this takes account of fcc following hcp 
#
      for ilist in range(0,num_to_set_hcp):
           dist_list_htf.append([j+num_to_set_hcp for j in range(0,num_to_set_fcc)])
#
#
      print("num_in_hcp_dist_list = ",num_in_hcp_dist_list)
      for ilist in range(0,num_to_set_hcp-1):
         print(dist_list_hcp[ilist])
#
      print("num_in_fcc_dist_list = ",num_in_fcc_dist_list)
      for ilist in range(0,num_to_set_fcc-1):
         print(dist_list_fcc[ilist])

      print("num_in_htf_dist_list = ",num_in_htf_dist_list)
      for ilist in range(0,num_to_set_hcp):
         print(dist_list_htf[ilist])

      hcp_dists=np.empty([num_in_hcp_dist_list])
      fcc_dists=np.empty([num_in_fcc_dist_list])
      htf_dists=np.empty([num_in_htf_dist_list])
#
      num_dists=0 
#
# -------------------------------------------------------------
# Create flag lists for the configurations on each sub-set
# Main look for structure creation, now a nested double loop.
# -------------------------------------------------------------
#
      print("Starting loops for this hcp/fcc population setting, first pass set")
      first_pass_inner=True                                      # note first pass of fcc loop
#
      for ihcp in range(0,nlist_hcp):
         hcp_struct_id =  flags_to_int(flag_list_hcp[ihcp])
         print(flag_list_hcp[ihcp],"hcp flag list : ", hcp_struct_id)
#
# Nested so that we consider every configuration on fcc with each hcp arrangement
#
         for ifcc in range(0,nlist_fcc):
            fcc_struct_id =  flags_to_int(flag_list_fcc[ifcc])
            print(flag_list_fcc[ifcc],"fcc flag list : ", fcc_struct_id)
#
# Create the next configuration by adding atoms to the clean surface according to
# flag_list_hcp/fcc
#
            this_atoms=clean_surf.copy()
#
# Add ads atoms to this_atoms as indicated by the flag_list_hcp/fcc settings
#
            for iatom in range(0,num_top_hcp):  
               if ( flag_list_hcp[ihcp][iatom] == 1 ): 
                   indx=top_hcp[iatom]
                   this_atoms.append(atoms[indx])
#
# Add the fcc atoms to the end of the this_ads_atoms set so that they can go as one
# block for the inter-atom distance from "get_distances"
#
            for iatom in range(0,num_top_fcc):  
               if ( flag_list_fcc[ifcc][iatom] == 1 ): 
                   indx=top_fcc[iatom]
                   this_atoms.append(atoms[indx])
#
# Make this_ads_atoms as atom set for just the ads_elem that have been added
#
            this_ads_atoms=this_atoms[[atom.index for atom in this_atoms if atom.symbol == ads_elem]]
#
# Work out the distance lists for this configuration
#
            dists_hcp=[]
            for ilist in range(0,num_to_set_hcp-1):
                dists_hcp=np.append(dists_hcp,this_ads_atoms.get_distances(ilist,dist_list_hcp[ilist],mic=True))
#
            dists_fcc=[]
            fcc_start=num_to_set_hcp
            for ilist in range(0,num_to_set_fcc-1):
                dists_fcc=np.append(dists_fcc,this_ads_atoms.get_distances(fcc_start+ilist,dist_list_fcc[ilist],mic=True))
#
            dists_htf=[]
            if (num_to_set_fcc > 0):
               for ilist in range(0,num_to_set_hcp):
                   dists_htf=np.append(dists_htf,this_ads_atoms.get_distances(ilist,dist_list_htf[ilist],mic=True))
#
# Reject structures with close contacts betweeen surface adsorbates.
# For alkylidene cut-off set to 2.0 A. This rejects occupancy of neighbouring hcp and fcc
# sites but allows fcc - fcc minimum distance.
            bother=True
            if (len(dists_hcp) > 0 and min(dists_hcp) < 2.0):
               bother=False
            elif  (len(dists_fcc) > 0 and min(dists_fcc) < 2.0):
               bother=False
            elif  (len(dists_htf) > 0 and min(dists_htf) < 2.0):
               bother=False
#
#            print(" ")
#            print("dists_hcp: ", dists_hcp)
#            print("dists_fcc: ", dists_fcc)
#            print("dists_htf: ", dists_htf)
#
            if (bother):
# Compare this set of distances with those already seen
# num_dists will hold the number of unique distance patterns found so far.
#
               if (first_pass_inner):
                  first_pass_inner=False
                  all_dists_hcp = dists_hcp
                  all_dists_fcc = dists_fcc
                  all_dists_htf = dists_htf
                  num_dists=1
                  unique_list=[[hcp_struct_id, fcc_struct_id]]
                  num_of_unique=[1]

                  print("First pass all_dists_hcp : ",all_dists_hcp)
                  print("First pass all_dists_fcc : ",all_dists_fcc)
                  print("First pass all_dists_htf : ",all_dists_htf)

                  ads_file="top_ads_atoms_%d_%d_%d.cif" % (num_to_set, hcp_struct_id, fcc_struct_id)
#                  write(ads_file, this_ads_atoms ,format='cif')
#
# Add methyl groups
#
                  add_methyl(this_atoms, top_metal_layer, ads_elem)
#
# Write the slab with adsorbates present
#                  ads_slab_file="ads_slab_%d_%d_%d.cif" % (num_to_set, hcp_struct_id, fcc_struct_id)
                  ads_slab_file_aims="ads_slab_%d_%d_%d.in" % (num_to_set, hcp_struct_id, fcc_struct_id)
                  write(ads_slab_file_aims, this_atoms,format='aims')
#
# Record data for this level of ads_elem overall
#
                  if (first_pass_outer):    
                    first_pass_outer=False
                    unique_list_tot=[[hcp_struct_id, fcc_struct_id]]
                    unique_filenames=[ads_slab_file_aims]
                    degen_unique=[1]
                  else:
                    unique_list_tot.append([hcp_struct_id, fcc_struct_id])
                    unique_filenames.append(ads_slab_file_aims)
                    degen_unique.append(1)
#
# Need to test to see if this configuration has been seen before by comparing the
# dist arrays with those already seen
#
               else:
                  ind_hcp=0
                  ind_fcc=0
                  ind_htf=0
                  matched=False
#  
#  Loop over the dist_lists that we already have
#  
                  for which_dist in range(0,num_dists):
                     test_dist_hcp=[]
                     test_dist_fcc=[]
                     test_dist_htf=[]
#  
#  copy the next distribution from the full list into the test_dist arrays
#  
#                  print("Checking against hcp_dist list with %d members" % num_in_hcp_dist_list)
                     for jdist in range(0,num_in_hcp_dist_list):
#                     print("ind_hcp = %d" % ind_hcp)
                        test_dist_hcp.append(all_dists_hcp[ind_hcp])
                        ind_hcp += 1
#  
#                  print("Checking against fcc_dist list with %d members" % num_in_fcc_dist_list)
                     for jdist in range(0,num_in_fcc_dist_list):
#                     print("ind_fcc = %d" % ind_fcc)
                        test_dist_fcc.append(all_dists_fcc[ind_fcc])
                        ind_fcc += 1
#  
#                  print("Checking against htf_dist list with %d members" % num_in_htf_dist_list)
                     if (num_to_set_fcc > 0):
                        for jdist in range(0,num_in_htf_dist_list):
#                        print("ind_htf = %d" % ind_htf)
                           test_dist_htf.append(all_dists_htf[ind_htf])
                           ind_htf += 1
#
#                  sys,exit("Just testing")
#
# To have seen before all three lists need to match 
#
                     matched_hcp = num_in_hcp_dist_list == 0 or match_vectors(dists_hcp,test_dist_hcp, 0.05)
                     matched_fcc = num_in_fcc_dist_list == 0 or match_vectors(dists_fcc,test_dist_fcc, 0.05)
                     matched_htf = num_in_htf_dist_list == 0 or match_vectors(dists_htf,test_dist_htf, 0.05)
#
                     if ( matched_hcp and  matched_fcc and  matched_htf ):
                        matched=True
                        print("Matched with existing list...")
                        num_of_unique[which_dist] += 1
                        degen_unique[nunique_tot+which_dist] += 1
                        break
#
# If this has not been matched must have a new unique configuration
#
                  if (not matched):
                     num_dists+=1
                     all_dists_hcp=np.append(all_dists_hcp,[dists_hcp])
                     all_dists_fcc=np.append(all_dists_fcc,[dists_fcc])
                     all_dists_htf=np.append(all_dists_htf,[dists_htf])
#
                     unique_list.append([hcp_struct_id,fcc_struct_id])
                     num_of_unique.append(1)
                     unique_list_tot.append([hcp_struct_id, fcc_struct_id])
                     degen_unique.append(1)
#
                     print("New unique configuration adding to lists")
#                  print("Now all_dists_hcp : ",all_dists_hcp)
#                  print("Now all_dists_fcc : ",all_dists_fcc)
#                  print("Now all_dists_htf : ",all_dists_htf)
#
                     ads_file="top_ads_atoms_%d_%d_%d.cif" % (num_to_set, hcp_struct_id, fcc_struct_id)
#                     write(ads_file, this_ads_atoms ,format='cif')
#
# Add methyl groups
#
                     add_methyl(this_atoms, top_metal_layer, ads_elem)
#
#
# Write the slab with adsorbates present
                     ads_slab_file="ads_slab_%d_%d_%d.cif" % (num_to_set, hcp_struct_id, fcc_struct_id)
                     ads_slab_file_aims="ads_slab_%d_%d_%d.in" % (num_to_set, hcp_struct_id, fcc_struct_id)
                     unique_filenames.append(ads_slab_file_aims)
                     write(ads_slab_file_aims, this_atoms,format='aims')
#
# End of bother condition for distance test
            else:
               print("Rejected due to close contact between adsorbates")
#
# End of ihcp loop, the outer of the loops, for this population of hcp and fcc sites
      nunique_tot+=num_dists
#
#
   print(" ")
   print(f"Thats %d patterns altogether of which %d are unique" % (nconfigs, nunique_tot))
#
# Report the data
#
   csv_headers=['index','filename','hcp structure index', 'fcc structure index', 'Number occurances']
   csv_filename="degen_%d.csv" % num_to_set
   csv_file = open(csv_filename, 'w', newline='')
   csv_writer = csv.writer(csv_file)
   csv_writer.writerow(csv_headers)

   for ind in range(0, nunique_tot):
#
      csv_row = [ind, unique_filenames[ind], unique_list_tot[ind][0], unique_list_tot[ind][1], degen_unique[ind]]
#
      csv_writer.writerow(csv_row)
#
   csv_file.close()
#
print(f"End of program")
#

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:55:21 2020

@author: dave
"""
import numpy as np
#from numpy.random import default_rng
import sys
import csv
import os
#
# Function to obtain HOME
from os.path import expanduser
#
from ase import Atom
from ase import Atoms

from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength

from ase.io import read, write, Trajectory

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
# Bonds dictionary
#
standard_bonds={ 'CC': 1.55, 'CO': 1.55, 'CH': 1.2, 'OH': 1.2 }
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
# Generate neighbour list for each atom in the list
#
def generate_neighbours(atoms):
#
    natoms=len(atoms)
    index_list=[i for i in range(0,natoms)]
    neighs=[]
#
    for iatom in range(0,natoms):
       dist_list=atoms.get_distances(iatom, index_list, mic=True)
       this_neighs=[]
#
       for jatom in range(0,natoms):
          if jatom != iatom:
             sym_pair= atoms.symbols[iatom] + atoms.symbols[jatom]
             if sym_pair not in standard_bonds:
                sym_pair= atoms.symbols[jatom] + atoms.symbols[iatom]
#
             if sym_pair in standard_bonds and dist_list[jatom] <= standard_bonds[sym_pair]:
                this_neighs.append(jatom)
       neighs.append(this_neighs)

    return neighs
#
# test for carbons in aromatic rings
#
def C_in_aromatic(atoms, neighs, iatom):
#
   nneighs= len(neighs)
   in_ring = False
# check the environment matches aromatic
   if nneighs == 3:
        in_ring = True
        for ineigh in range(0,nneighs):
           if atoms.symbols[neighs[ineigh]] != "C" and \
              atoms.symbols[neighs[ineigh]] != "H" :
                  in_ring=False
#
   return in_ring

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
atoms=read(geom_file)
natoms=len(atoms)
in_aromatic=[]
for iatom in range(0,natoms):
   in_aromatic.append(False)
#
#atoms.cell=[103.9,103.9,103.9]
atoms.cell=[103.0,103.0,103.0]
atoms.pbc=True
#
csv_filename="ring_atom_indices.csv"

#
neighs= generate_neighbours(atoms)
#
#for iatom in range(0,natoms):
#    print(f"%d %s %10.6f %10.6f %10.6f"                              \
#          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
#             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#    print(neighs[iatom])
#
# Find carbon rings
#
iatom=0
ring_lists=[]
while iatom < natoms:
#  print("looking for rings from atom ", iatom)
#
  if not in_aromatic[iatom] and atoms.symbols[iatom] == "C" :
     in_ring = C_in_aromatic(atoms, neighs[iatom], iatom)
#
# flag all of ring and capture list of atoms
#
     this_ring = []
     if in_ring:
#        print("This is a ring atom")
        in_aromatic[iatom] = True
        this_ring.append(iatom)
#
        for ineigh in range(0,3):
           neigh_indx = neighs[iatom][ineigh]
           if not in_aromatic[neigh_indx] and C_in_aromatic(atoms, neighs[neigh_indx], neigh_indx):       
              in_aromatic[neigh_indx] = True
              this_ring.append(neigh_indx)
#              print("So is it's neighbour:", atoms.symbols[neigh_indx])
# test these neighbours too
              for next_neigh in range(0,3):
                 next_indx=neighs[neigh_indx][next_neigh]
                 if not in_aromatic[next_indx] and C_in_aromatic(atoms, neighs[next_indx], next_indx):       
                    in_aromatic[next_indx] = True
                    this_ring.append(next_indx)
#                    print("So is it's next  neighbour:", atoms.symbols[next_indx])
# and next nearest neighs
                    for next2_neigh in range(0,3):
                       next2_indx=neighs[next_indx][next2_neigh]
                       if not in_aromatic[next2_indx] and C_in_aromatic(atoms, neighs[next2_indx], next2_indx):
                          in_aromatic[next2_indx] = True
                          this_ring.append(next2_indx)
#                          print("So is that ones neighbour:", atoms.symbols[next2_indx])

     if len(this_ring) > 0 :
         ring_lists.append(this_ring)
#
  iatom+=1
#
print("Found ", len(ring_lists), " aromatic rings.")
#for iring in range(0, len(ring_lists)):
#   print("Aromatic ring atoms: ", ring_lists[iring])
print("----------------===============C atoms found==========--------------")
#
# Capture ring atoms only
just_ring_coords=[]
just_ring_symbols=[]
just_small_ring_coords=[]
just_small_ring_symbols=[]
#
csv_file = open(csv_filename, 'w', newline='')
csv_writer = csv.writer(csv_file)
#
# make axis atoms the first two in each ring atom list
#
for iring in range(0, len(ring_lists)):
#   print(ring_lists[iring])
#
   axis_atoms=[]
   for iatom in ring_lists[iring]:
#
     is_axis=True 
     for ineigh in range(0,len(neighs[iatom])):
       neigh_indx = neighs[iatom][ineigh]
       if atoms.symbols[neigh_indx] != "C":
          is_axis=False
          break
     if is_axis:
       axis_atoms.append(iatom)
#  print("Axis atoms: ", axis_atoms)
#
#  print("ring_list before neighs added: ", ring_lists[iring])
   ordered=axis_atoms.copy()
   for iatom in ring_lists[iring]:
     keep=True
     for icheck in ordered:
        if iatom == icheck:
           keep=False
#       
     if keep:
        ordered.append(iatom)
#
     for neigh_indx in neighs[iatom]:
        keep=True
        for icheck in ordered:
           if neigh_indx == icheck:
             keep=False
        if keep:
             ordered.append(neigh_indx)
#
   ring_lists[iring]=ordered.copy()
#   print("csv_row: ", ring_lists[iring])
#   ring_symbols=[]
#   for iatom in ring_lists[iring]:
#      ring_symbols.append(atoms.symbols[iatom])
#   print("labels: ", ring_symbols)
#
   csv_writer.writerow(ring_lists[iring])  
#   
   for ilist in range(0,len(ring_lists[iring])):
     iatom=ring_lists[iring][ilist]
#
     just_ring_coords.append(atoms.positions[iatom])
     just_ring_symbols.append(atoms.symbols[iatom])
#     if len(ring_lists[iring]) < 10:
#         just_small_ring_coords.append(atoms.positions[iatom])
#         just_small_ring_symbols.append(atoms.symbols[iatom])
#
just_ring_atoms = Atoms(just_ring_symbols, just_ring_coords)
just_ring_atoms.set_cell(atoms.get_cell())

write('check.cif', atoms, format="cif")
write('check_rings.cif', just_ring_atoms, format="cif")

csv_file.close()
#
# Now perform the rotations
# get random theta
#
#rng=default_rng()
rnds=np.random.uniform(0,1,len(ring_lists))
#
for iring in range(0, len(ring_lists)):
#   print("rotating ring %d" % iring)
#   print(ring_lists[iring])
#
#  shift all atoms so that first is at origin
   shift_vec=atoms.positions[ring_lists[iring][0]].copy()
   dists=atoms.get_distances(ring_lists[iring][0],ring_lists[iring], mic=False)
   shift_vecs=[]
#
   for idist in range(0,len(dists)):
      if (dists[idist] > 10.0):
         vec_mic=atoms.get_distance(ring_lists[iring][0],ring_lists[iring][idist], mic=True, vector=True)
         vec_no_mic=atoms.get_distance(ring_lists[iring][0],ring_lists[iring][idist], mic=False, vector=True)
         pbc_shift=np.subtract(vec_no_mic,vec_mic)
      else:
         pbc_shift=np.zeros(3)
      shift_vecs.append(np.add(shift_vec,pbc_shift))
#
   if iring == 5:
      check5_symbols=[]
      check5_coords =[]
      print("ring 5 original dists: ", dists)

      for iatom in ring_lists[iring]:
         check5_symbols.append(atoms.symbols[iatom])
         check5_coords.append(atoms.positions[iatom])
#
      check5_ring_atoms = Atoms(check5_symbols, check5_coords)
      check5_ring_atoms.set_cell(atoms.get_cell())
      write('check5_1.cif', check5_ring_atoms, format="cif")
#
   for iii in range(0,len(ring_lists[iring])):
      iatom=ring_lists[iring][iii]
      if iring == 5:
         print("Atom ", iatom, " shift by: ", shift_vecs[iii])
      for jjj in range(0,3):
         atoms.positions[iatom][jjj]=atoms.positions[iatom][jjj]-shift_vecs[iii][jjj]
#
   if iring == 5:
      dists=atoms.get_distances(ring_lists[iring][0],ring_lists[iring], mic=False)
      print("ring 5 dists after shift: ", dists)
      check5_symbols=[]
      check5_coords =[]
      for iatom in ring_lists[iring]:
         check5_symbols.append(atoms.symbols[iatom])
         check5_coords.append(atoms.positions[iatom])
#
      check5_ring_atoms = Atoms(check5_symbols, check5_coords)
      check5_ring_atoms.set_cell(atoms.get_cell())
      write('check5_2.cif', check5_ring_atoms, format="cif")
#
#
# ring lists should be defined so that the first two atoms are the axis atoms
#
   axis=atoms.get_distance(ring_lists[iring][0], ring_lists[iring][1], mic=False, vector=True)
#   print("Axis: ", axis)
#
# Get normalised axis vector
#
   u = axis/np.linalg.norm(axis)
#   print("Normalised Axis: ", u)
#
   theta=2.0*np.pi*(0.5-rnds[iring])
#   print("Using ",rnds[iring], "theta = ", theta/np.pi, " pi")
#
#   quat=[np.cos(theta/2), u[i]*np.sin(theta/2) for i in range(0,3)]
   s=np.sin(theta/2)
   quat=[np.cos(theta/2), u[0]*s, u[1]*s, u[2]*s] 
#   print("Quat:      ", quat)   
#   print("Quat size: ", np.linalg.norm(quat))
#
# define rotation matrix
#
   rotmat=np.array([[ 1-2*(quat[2]*quat[2]+quat[3]*quat[3]),   2*(quat[1]*quat[2]-quat[3]*quat[0]),    2*(quat[1]*quat[3]+quat[2]*quat[0])],\
                    [   2*(quat[1]*quat[2]+quat[3]*quat[0]), 1-2*(quat[1]*quat[1]+quat[3]*quat[3]),    2*(quat[2]*quat[3]-quat[1]*quat[0])],\
                    [   2*(quat[1]*quat[3]-quat[2]*quat[0]),   2*(quat[2]*quat[3]+quat[1]*quat[0]),  1-2*(quat[1]*quat[1]+quat[2]*quat[2])]])
#
#   print("Rotmat: ", rotmat)
#
# Do the rotation
#
   for iii in range(0,len(ring_lists[iring])):
     iatom=ring_lists[iring][iii]
     newpos= rotmat.dot(atoms.positions[iatom])
     atoms.positions[iatom]=np.add(newpos,shift_vecs[iii])
#
   if iring == 5:
      check5_symbols=[]
      check5_coords =[]
      dists=atoms.get_distances(ring_lists[iring][0],ring_lists[iring], mic=False)
      print("ring 5 dists after rot: ", dists)
      for iatom in ring_lists[iring]:
         check5_symbols.append(atoms.symbols[iatom])
         check5_coords.append(atoms.positions[iatom])
#
      check5_ring_atoms = Atoms(check5_symbols, check5_coords)
      check5_ring_atoms.set_cell(atoms.get_cell())
      write('check5_3.cif', check5_ring_atoms, format="cif")
#
#
write('check_rots.cif', atoms, format="cif")

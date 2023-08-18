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

from ase.io import read, write, Trajectory

from ase.calculators.vasp import Vasp2
#
# count the frames in an XDATCAR file
#
def count_xdatcar_frames(xdatcar_file):
#
   fptr=open(xdatcar_file, 'rt')
   lines=fptr.readlines()
#
   nframes=0
   for line in lines:
      if ( 'config' in line ):
         nframes += 1
#
   fptr.close()
#
   return nframes
#
#
# Write the body of a csv file from a list
#
def csv_body_write(fptr, csv_data):

   for row in csv_data:
      for i in range(0,len(row)):
         fptr.write(str(row[i]))
         fptr.write(',')
      fptr.write('\n')
#
   return
#
# Take running average of a list
#
def running_average(data, window):
#
    ndata = len(data)
    if window > ndata:
       print("ERROR: Asked for a window size greater than data length in running_average.")
#
    tot=0
    run_av=[]
    for i in range(0,ndata):
      if i+1 <= window:
        tot += data[i]
        np=i+1
      else:
        tot += data[i] - data[i-window]
        np=window
#
      run_av.append(tot/np)
#
    return run_av
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
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
nargs=len(sys.argv)
#
if nargs != 4:
   print("ERROR:") 
   print("ERROR: Incorrect number of arguments with script.")
   print("ERROR: Syntax is:")
   print("ERROR: python3 ase_vasp_md_dists.py XDATCAR atom1 atom2")
   print("ERROR:") 
   print("ERROR: where \"atom1\\2\" are the indices of atoms defining")
   print("ERROR: the distance to plot.") 
   print("ERROR:") 
   exit(0)
# 
xdatcar_file = sys.argv[1]    # Name of the MD VASP trajectory file that we will analyse, 
#                             # taken from command line
iatom1 = int(sys.argv[2])     # Define the atoms to monitor
iatom2 = int(sys.argv[3])
#
print(f"------------------------------------------------------------")
print(f"VASP trajectory taken from   : %s" % xdatcar_file )
print(f"------------------------------------------------------------")

traj = Trajectory("test.traj", mode='a')

nframes = count_xdatcar_frames(xdatcar_file)
print("Trajectory contains %d frames." % nframes)

fram_inds="0:%d" % nframes
print(fram_inds)
traj=read(xdatcar_file, index=fram_inds, format='vasp-xdatcar')
#
#
dist_list=[]
for atoms in traj:
  dist_list.append(atoms.get_distance(iatom1, iatom2, mic=True)) 
#
# Work out running average
#
npoints=100
run_av=running_average(dist_list, npoints)
#
print("len dist_list %d len run_av %d" % (len(dist_list),len(run_av)))
#
#
# create a csv file for plotting
#
csv_fptr = open('dist_vs_time.csv', 'w', newline='')
csv_fptr.write('Step , distance / Angs, Running average (%d pts) / Angs \n' % npoints)
#
# write the data
csv_data=[]
for i in range(0,len(dist_list)):
#    csv_fptr.write('%d , %10.6f, %10.6f \n' % (i, dist_list[i], run_av[i]))
   csv_data.append([i , dist_list[i], run_av[i]])
#
csv_body_write(csv_fptr, csv_data)
#
csv_fptr.close()
#
print("Average (std) : %10.6f (%10.6f) Angs" % (np.average(dist_list), np.std(dist_list))) 
print("")
print("Run completed.......")
print("")


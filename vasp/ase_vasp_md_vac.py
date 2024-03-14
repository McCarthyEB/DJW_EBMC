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
# Read a VDATCAR file for velocities
#
def read_vdatcar(vdatcar_file, velocities):
#
   fptr=open(vdatcar_file, 'rt')
   lines=fptr.readlines()
#
   ic=0
   nframes=0
   vframe=[]
   ion_num=[]
   ion_types=[]
   new_frame=False
#
   for iline in range(0,len(lines)):
      words=lines[iline].split()
      nwords=len(words)
      if iline == 5:
         for word in words:
            ion_types.append(word)
      if iline == 6:
         for word in words:
            ion_num.append(int(word))
#
      if iline > 7: 
        if nwords == 3:
          vframe.append([float(words[0]),float(words[1]),float(words[2])])
          new_frame=True
          ic+=1
        elif new_frame:
          nframes+=1
          new_frame=False
          velocities.append(vframe)
          vframe=[]
          
#
   print(ion_types)      
   print(ion_num)      
    

   return nframes
#
def vac_calc(atoms,vels,nvframes,label):
#
   vac=[]
   inds=[]
#
# Get indices for atoms that match label
#
   natoms=len(atoms)
   for i in range(0,natoms):
      if atoms.symbols[i] == label:
         inds.append(i)
#
# Only use 80% to get time origin stats
#
   for ifr in range(1,int(0.8*nvframes)):
      nd=0
      dot=0
      for i0 in range(0,nvframes-ifr):
        it=i0+ifr
#
        frame_0=vels[i0]
        frame_t=vels[it]
#
        for iuse in inds:
          dot+=np.dot(frame_0[iuse],frame_t[iuse])
          nd+=1
#
      vac.append(dot/nd) 
#          
   return vac
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
# Main code begins
#
nargs=len(sys.argv)
#
if nargs != 4:
   print("ERROR:") 
   print("ERROR: Incorrect number of arguments with script.")
   print("ERROR: Syntax is:")
   print("ERROR: python3 ase_vasp_md_vac.py XDATCAR VDATCAR elem")
   print("ERROR:") 
   print("ERROR: where \"elem\" is the element for the VAC function.") 
   print("ERROR:") 
   print("ERROR: If you don\'t have a VDATCAR file build using the ") 
   print("ERROR: vtst script xdat2vdat.pl. ") 
   print("ERROR:") 
   exit(0)
# 
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
xdatcar_file = sys.argv[1]    # Name of the MD XDATCAR
vdatcar_file = sys.argv[2]    # Name of the MD VDATCAR, if not present create from XDATCAR
#                             # using the vtstscript xdat2vdat.pl
elem_select = sys.argv[3]     # Element for which the velocity autocorrelation function is
                              # required
#
print(f"------------------------------------------------------------")
print(f"VASP trajectory taken from   : %s" % xdatcar_file )
print(f"VASP velocities taken from   : %s" % vdatcar_file )
print(f"------------------------------------------------------------")
#
atoms=read(xdatcar_file, format='vasp-xdatcar')
#
# Look at velocities file
vels=[]
nvframes=read_vdatcar(vdatcar_file,vels)
print("len vels %d" % len(vels))
#
# calculate velocity autocorrelation for selected atom types
#
vac=vac_calc(atoms,vels,nvframes,elem_select)
print("Len vac: %d\n" % len(vac))
#
csv_fptr = open('vac_vs_time.csv', 'w', newline='')
csv_fptr.write('time , vac \n')
#
for i in range(0,len(vac)):
    csv_fptr.write("%d , %10.6f \n" % (i, vac[i]))   
#
# Fourier transform to get frequency spectrum
#
spectrum=np.fft.fft(vac)
freqs=np.fft.fftfreq(len(vac))
#
# create a csv file for plotting
#
csv_fptr = open('spectrum.csv', 'w', newline='')
csv_fptr.write('index, freq / fs-1 , wavenum cm-1, power \n')
# 
# write the data, use speed of light to get wavenumbers
# also include conversion to s-1 units                 
# only include positive frequencies in spectrum        
cspeed=2.998E-5
for i in range(0,int(len(spectrum)/2)):
    csv_fptr.write("%d, %10.6f, %10.6f, %10.6f \n" % (i, freqs[i], freqs[i]/cspeed, spectrum[i]))   
#
csv_fptr.close()
#
print("")
print("Run completed.......")
print("")


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
# Read data from the OUTCAR file      
#
def read_md_outcar(outcar_file):
#
   fptr=open(outcar_file, 'rt')
   lines=fptr.readlines()
# 
# Set default values
   timestep = -1
   tein = -1
   tebeg = -1
#
# lists
#
   temperature = []
   toten       = []
   ekin        = []
   nose_lines  = []
#
   nframes=0
   for line in lines:
     if ( 'POTIM' in line ):
       words  = line.split()
       nwords = len(words)
       if nwords > 3:
          timestep = float(words[2])
#
     elif ('TEIN' in line ):
       words  = line.split()
       nwords = len(words)
       if nwords > 3:
          tein = float(words[2])
#
     elif ('TEBEG' in line ):
       words  = line.split()
       nwords = len(words)
       if nwords > 3:
          word4val = words[2].split(';')[0]
          tebeg = float(word4val)
#
     elif ('ion-electron' in line and 'TOTEN' in line):
       words  = line.split()
       nwords = len(words)
       if nwords > 4:
          toten.append(float(words[4]))
#
     elif ('kinetic' in line and 'EKIN' in line):
       words  = line.split()
       nwords = len(words)
       if nwords > 4:
          ekin.append(float(words[4]))
#
     elif ('temperature' in line ):
       words  = line.split()
       nwords = len(words)
       if nwords > 6:
          temperature.append(float(words[5]))
#
     elif ('Nose' in line ):
       nose_lines.append(line)
#
   fptr.close()
#
   return timestep, tein, tebeg, temperature, toten, ekin, nose_lines
#
# Write the body of a csv file from a list
#
def csv_body_write(fptr, csv_data):
  
   for row in csv_data:
      for i in range(0,len(row)):
         fptr.write(str(row[i]))
         fptr.write(',')
      fptr.write('\n')
        

   return 
#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
nargs=len(sys.argv)
#
if nargs != 2:
   print("ERROR:") 
   print("ERROR: Incorrect number of arguments with script.")
   print("ERROR: Syntax is:")
   print("ERROR: python3 ase_vasp_md_vac.py OUTCAR")
   print("ERROR:") 
   exit(0)
# 
outcar_file = sys.argv[1]     # Name of the MD VASP trajectory file that we will analyse, 
#                             # taken from command line
#
print(f"------------------------------------------------------------")
print(f"VASP outcar taken as file    : %s" % outcar_file )
print(f"------------------------------------------------------------")

timestep, tein, tebeg, temperature, toten, ekin, nose_lines = read_md_outcar(outcar_file)
print("Timestep for run  : %10.6f femtoseconds." % timestep)
print("Temperature start target  : %10.6f K." % tein)
print("Temperature end target    : %10.6f K." % tebeg)
#
# Writes lists to csv file
#
csv_fptr = open('MD_data.csv', 'w', newline='')
csv_fptr.write('Step , time / fs, PE / eV, KE / eV, temperature / K \n')
#
csv_data=[]
for ind in range(0,len(temperature)):
  csv_data.append([ind, timestep*ind, toten[ind], ekin[ind], temperature[ind]])
#
csv_body_write(csv_fptr, csv_data) 
#
csv_fptr.close()
#
print("Average temperature (std) : %10.6f (%10.6f) K" % (np.average(temperature), np.std(temperature))) 
#
if len(nose_lines) > 0 :
    print("")
    print("Nose Thermostat information: ")
    for line in nose_lines:
       print(line)
#
print("")
print("Run completed.......")
print("")


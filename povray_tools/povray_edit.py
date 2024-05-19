import numpy as np
import sys
import csv
import os
import glob
#
# Function to obtain HOME
from os.path import expanduser
#
# Make modules available
#
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
from binary_tools import *
from vectors import *
#
# Look for position of a keyword in an array representation of a line
#
def key_place(key, list):
#
   nlist=len(list)
   key_place=-1
   for iword in range(0,nlist):
      if key in linesplit[iword]:
          key_place=iword
#
   return key_place
#
# Use atom colours to define element
# The rgb_place is the location of the rgb key word in the line passed in list
#
def rgb_to_elem(rgb_place, list):
#
   if    '0.00' in list[rgb_place+1] \
     and '0.41' in list[rgb_place+2] \
     and '0.52' in list[rgb_place+3]:
        return "Pd"
#
   elif  '0.56' in list[rgb_place+1] \
     and '0.56' in list[rgb_place+2] \
     and '0.56' in list[rgb_place+3]:
        return "C"
#
   elif  '1.00' in list[rgb_place+1] \
     and '1.00' in list[rgb_place+2] \
     and '1.00' in list[rgb_place+3]:
        return "H"
#
   else:
     print("ERROR: rgb settings passed to rgb_to_elem that cannot be matched")
     print("ERROR: to an element: Update code...............................")
     print("ERROR: starting from rgb_place %d in list:" % rgb_place)
     print(list)
     exit(0)
#
# Extract co-ordinates as float triple given start point in list
#
#
def extract_coords(place, list):
#
# x coord
      temp=linesplit[place].split(",")
      coords=[]
      coords.append(float(temp[0]))
#
# y coord
      temp=linesplit[place+1].split(",")
      coords.append(float(temp[0]))
#
# z coord
      temp=linesplit[place+2].split(">")
      coords.append(float(temp[0]))
#
      return coords
#
# Start of main code
#
# define the file to edit on command line, i.e. run with:
# python3 povray_edit.py my_original.pov
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
pov_file=sys.argv[1]
print("Will work with povray file: %s" % pov_file)
#
fileobj=open(pov_file, 'r')
#
temp=pov_file.split(".")
if len(temp) > 1:
  stem=temp[0]
  print("Stem: %s" % stem)
#
pov_file_new=stem + "_new.pov"
print("Will produce an new version: %s" % pov_file_new)
fileobj_new=open(pov_file_new, 'w')
#
# Copy the ini file
#
ini_file= stem + ".ini"
ini_file_new=stem + "_new.ini"
copy_cmd="cp " + ini_file + " " + ini_file_new
os.popen(copy_cmd)
#
lines = fileobj.readlines()
nlines = len(lines)  
print("File has %d lines" % nlines)
##
## This loop will now run through all the lines of the file so need to
## extract information from the lines as appropriate, hence the if blocks
##
ncyl=0
symbols=[]
positions=[]
for iline in range(len(lines)):
   line=lines[iline]
   linesplit = line.split()
   nwords = len(linesplit)
   del_this=False
#
# Deal with atoms
#
   if nwords > 0 and "atom" in linesplit[0]:
#
      rgb_place= key_place('rgb', linesplit)
      symbols.append(rgb_to_elem(rgb_place, linesplit))
      print("Found a %s atom..." % symbols[-1])
      positions.append(extract_coords(1,linesplit))
      print("at coords: ", positions[-1])      
#
# Deal with bonds
#
   if nwords > 0 and "cylinder" in linesplit[0]:
#
# Check not the cell boundaries:
      is_bond=True
      for iword in range(0,nwords):
         if 'cell' in linesplit[iword]:
            is_bond=False
#
      if is_bond:
#
         ncyl = ncyl+1
#
# Remove Pd..Pd bonds, identified with {color rgb <0.00, 0.41, 0.52> transmit 0.0} 
#
         rgb_place= key_place('rgb', linesplit)
#
# Find Pd..Pd bonds and flag for deletion
#
         if  rgb_to_elem(rgb_place, linesplit) == 'Pd':
# Check ends are both Pd
           start=extract_coords(2,linesplit)
           end=extract_coords(6,linesplit)
#
           print("Bond with Pd starts ", start, " ends ", end)
#
           iend=-1
           istart=-1
           min_start_diff=1000.0
           min_end_diff=1000.0
           start_matched=False
           end_matched=False
#
           for iatom in range(0,len(positions)):
#
# If matched then use that, otherwise use closest atom
              match, diff = match_coords(start,positions[iatom],0.01)
#
              if match: 
                 istart=iatom
                 start_matched=True
              elif not start_matched:
                 if 'Pd' not in symbols[iatom] and diff < min_start_diff:
                    min_start_diff = diff
                    istart=iatom
#
# Note that the bond cylinders start at an atom co-ordinate but end short of the other atom
# Also at periodic boundaries get some bonds out to nothing
# So find shortest distance to non-Pd atom
#
              match, diff = match_coords(end,positions[iatom],0.01)
#
              if match: 
                 iend=iatom
                 end_matched=True
              elif not end_matched:
                 if 'Pd' not in symbols[iatom] and diff < min_end_diff:
                    min_end_diff = diff
                    iend=iatom
#
           print("Starts on atom %d (%10.6f) ends on %d (%10.6f)" % ( istart, min_start_diff, iend, min_end_diff))
#
           if istart > 0 and iend > 0:
              if istart == iend:
                 del_this= True
#
              else:
                 dist = vector_dist(positions[istart], positions[iend])
                 if 'Pd' in symbols[istart] and dist > 1.9: 
                    print("Found long Pd..X, flagging for deletion....")
                    del_this=True
#
   if not del_this:
      fileobj_new.write(line) 
#
print("Found %d cylinders...." % ncyl)
#
fileobj.close()
fileobj_new.close()

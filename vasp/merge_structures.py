#
import sys
import os
import numpy as np
from ase.lattice import *
from ase.visualize import view
from ase.build import make_supercell
#
# planewave cut-off scan to identify the best planewave cut-off sampling to use.
# The script expects the structure to be in a POSCAR file in the 
# structure sub-directory.
#
from ase import Atom
from ase import Atoms

from ase.io import read, write, Trajectory

#### Functionality to create all surfaces ####
#from carmm.build.facets import generate

#### Functionality to create surfaces via pymatgen ####
from pymatgen.io.cif import CifParser
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

#
# Bubble sort for atom list, this will sort atoms according to their
# co-ordinates in x,y,z according to the dir vector (0,0,1)
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

#Tried turning off dot product
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
# Read in the structure from a command line supplied file name
#
slab_file=sys.argv[1]
ads_file=sys.argv[2]

print("Reading slab from      : %s" % slab_file)
print("Reading adsorbate from : %s" % ads_file)
print("Will make combined     : POSCAR_slab_ads")
#
atoms= read(slab_file)
ads_atoms = read(ads_file)
#
slab_cent=atoms.get_center_of_mass()
ads_cent=ads_atoms.get_center_of_mass()
#
print("Slab structure: cofm: ", slab_cent)
#
zmax=-1000
zmin=1000
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
    if  atoms.positions[iatom][2] > zmax:
       zmax =  atoms.positions[iatom][2]
#
    if  atoms.positions[iatom][2] < zmin:
       zmin =  atoms.positions[iatom][2]
#
# slab thickness
#
slab_thick = zmax - zmin
#
zdiff = ads_cent[2]-slab_cent[2]
#
# The shift_x should move the adsorbate so that its cofm is ads_height above
# the upper slab surface
#
height=3.0
shift_z = height + 0.5*slab_thick - zdiff 
#
print("Adsorbate structure:")
for iatom in range(0,len(ads_atoms)):
    ads_atoms.positions[iatom][2] = ads_atoms.positions[iatom][2] + shift_z
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, ads_atoms.symbols[iatom], ads_atoms.positions[iatom][0], \
             ads_atoms.positions[iatom][1], ads_atoms.positions[iatom][2]))

    atoms.append(ads_atoms[iatom])

write("POSCAR_slab_ads", atoms, format='vasp') 

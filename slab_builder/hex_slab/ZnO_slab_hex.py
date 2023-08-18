import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from ase.units import kJ
from ase.eos import EquationOfState
from ase.io import aims
from ase.lattice import *
#
# planewave cut-off scan to identify the best planewave cut-off sampling to use.
# The script expects the structure to be in a POSCAR file in the 
# structure sub-directory.
#

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

#### Functionality to create all surfaces ####
from carmm.build.facets import generate


from ase.units import kJ
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
# tol is a tolerance in Angstroms for atoms to be accepted into a layer
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
# Read in the structure from a fhiaims file
#
geom_file="geometry.in"
cif_file="struct.cif"             

if os.path.isfile(geom_file):
   print("reading FHI-aims geometry.in file")
   atoms=read(geom_file)
#
elif os.path.isfile(cif_file):
   print("reading cif file")
   atoms=read(cif_file)
#
else:
   print("ERROR : No structure file, neither struct.car or castep.cell could be found")
   exit()
#
print("Bulk structure:")
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))

##turn Miller-Bravais indices into Miller indices
u=1
v=0
t=-1
w=0

slab_label="slab_"
#
miller_brav=[u,v,t,w]
for mil_ind in miller_brav:
   if mil_ind >= 0:
      slab_label+= "%d" % mil_ind
   else:
      abs_ind = -1*mil_ind
      slab_label+= "%d" % abs_ind + "b"
#
print("Original Miller-Bravais indices:", miller_brav)
print("Will label surface as          :", slab_label)
L = [[2,-1],[-1,2]]
R = [3*u,3*v]

res = np.linalg.inv(L).dot(R)

miller_index=[]
miller_index.append(int(res[0]))
miller_index.append(int(res[1]))
miller_index.append(w)

print("Updated Miller indices:", miller_index)


for ilayer in range(1,2):

#   facet, slab = generate(atoms, layers=ilayer, facets=[(miller_index[0],miller_index[1],miller_index[2])], vacuum=15, save=False)
   facet, slab = generate(atoms, layers=ilayer, facets=[miller_index], vacuum=15, save=False)
   filename=slab_label+"_l"+str(ilayer)+".in"

   write(filename,slab[0])

   print("Slab structure:")
   slab_atoms=slab[0]
   for iatom in range(0,len(slab_atoms)):
       print(f"%d %s %10.6f %10.6f %10.6f"                              \
             % (iatom, slab_atoms.symbols[iatom], slab_atoms.positions[iatom][0], \
                slab_atoms.positions[iatom][1], slab_atoms.positions[iatom][2]))
#
# Detail layers
#
   layers = find_layers(slab_atoms,0.1)
   num_layers=len(layers)
#
   for ilayer in range(0,num_layers):
      print("layer %d" % ilayer )
      for iatom in range(0,len(layers[ilayer])):
         indx=layers[ilayer][iatom]
         print("%d %s %10.6f %10.6f %10.6f " % ( indx, slab_atoms.symbols[indx], \
                slab_atoms.positions[indx][0], slab_atoms.positions[indx][1],         \
                slab_atoms.positions[indx][2] ))


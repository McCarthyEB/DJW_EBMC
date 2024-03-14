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
geom_file=sys.argv[1]

print("Reading geom_file from : %s" % geom_file)

bulk_struct = Structure.from_file(geom_file)

bulk_struct.to(filename="bulk_check.cif")
#
atoms= AseAtomsAdaptor().get_atoms(bulk_struct)       
#
print("Bulk structure:")
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))

####pymatgen test####
###This will generate your surface from the bulk file geometry.in you have already optimised.
### Before running this code, in your command line write
###
### module purge
### module load ase/3.20.1 python/3.7.0 pymatgen/2018.11.30
###
### This will allow the script to access the information it needs to run locally.
### Run the script by writing: python3 (thisscriptname).py
###
### The range parameters tell the script how many layers to create.
### The miller_index parameters can be changed to create different surfaces: for hexagonal cells such as ZnO, just use thefirst, second and fourth indices. 
### 
### With metals like Pd it looks like SlabGenerator will only give slabs that are integer numbers of the repeat
### in Z, so 3 atomic layers for Pd(111). So once slab is built will need own "find_layers" routines to get the
### number of desired layers.
###
miller=(2,1,1)
filename_base="POSCAR_%d%d%d_" % (miller[0],miller[1],miller[2])

slabgen = SlabGenerator(bulk_struct,
                        miller_index=miller,
                        min_slab_size=10.0,
                        lll_reduce=False,
                        center_slab=True,
                        in_unit_planes=False,
                        min_vacuum_size=20,
                        primitive=True,
                        max_normal_search=5,
                        reorient_lattice=True)
#
slab_list = slabgen.get_slabs()

#
for icut, slab in enumerate(slab_list):

   slab = slab_list[icut]
#
# Convert slab to ASE atoms object
#
   symbols=[str(site.specie.symbol) for site in slab]
   positions=[site.coords for site in slab]
   pbc=True
   cell=slab.lattice.matrix
   slab_atoms=Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)
# Sort the atoms along the z-direction
#
   bubble_atoms(slab_atoms,dir=(0,0,1))
#
   print("\nSlab Cell:")
   for icell in range(0,3):
       print("%10.6f %10.6f %10.6f" % \
           (slab_atoms.cell[icell][0],slab_atoms.cell[icell][1],slab_atoms.cell[icell][2]))

   print("Slab atoms:")
   for iatom in range(0,len(slab_atoms)):
       print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, slab_atoms.symbols[iatom], slab_atoms.positions[iatom][0], \
             slab_atoms.positions[iatom][1], slab_atoms.positions[iatom][2]))
#
   write(filename_base + "slab_cut%d" % icut, slab_atoms, format='vasp') 
#
# Find layers so we can trim to the desired thickness
#
   layers=find_layers(slab_atoms,0.1)
   nlayers=len(layers)
#
   print("Found %d layers..." % nlayers)
#
# Define number of atomic layers required
#
   need_layers=15
   del_list=[]
   for ilay in range(0,nlayers):
      print(layers[ilay])
      if ilay > need_layers:
         for jjj in range(0,len(layers[ilay])):
            del_list.append(layers[ilay][jjj])
#
   del slab_atoms[del_list]
#
##To  generate supercells:
   matrix=[[3,0,0],[0,2,0],[0,0,1]]
   sup_atoms=make_supercell(slab_atoms, matrix)
   write(filename_base + "sup_cut%d" % icut, sup_atoms, format='vasp')
#
   print("\nSuper Cell: N=%d" % len(sup_atoms))
   for icell in range(0,3):
       print("%10.6f %10.6f %10.6f" % \
           (sup_atoms.cell[icell][0],sup_atoms.cell[icell][1],sup_atoms.cell[icell][2]))
#
#
#repeat zfix for supercell
#
#      fix_list_sup=[]
#      for iatom in range(0,len(sup)):
#         if sup.positions[iatom][2] < zfix:
#            fix_list_sup.append(iatom)
#      cons_sup = FixAtoms(fix_list_sup)
#      sup.set_constraint(cons_sup)
#      write('supercell_slabheight_%s.cif' % (ilayer), sup, format='cif')
#      write('supercell_slabheight_%s.in' % (ilayer), sup, format='aims')


##Use these instead if you want to write all possible cuts in the z direction, rather than just the last one:
##Creates a lot of files! 
#       write('supercell_slabheight_%s_%s.cif' % (ilayer, n), sup, format='cif') 
#       write('supercell_slabheight_%s_%s.in' % (ilayer, n), sup, format='aims') 



#print(slabs[0])

       #cut0=slabs[0]
       #cut0.to(filename="cut0_slabheight_%s.cif" % (ilayer))
       #cut1=slabs[1]
       #cut1.to(filename="cut1_slabheight_%s.cif" % (ilayer))

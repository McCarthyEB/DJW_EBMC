#
import sys
import csv
import os
import numpy as np
#
from ase import Atom
from ase import Atoms
from ase.io import aims
from ase.lattice import *
from ase import neighborlist

from ase.io import read, write, Trajectory
from ase.io.vasp import *

#### Functionality to create all surfaces ####
#from carmm.build.facets import generate

# Function to obtain HOME
from os.path import expanduser
#
home = expanduser("~")
modules_dir = "%s/python/DJW_EBMC/modules" % home
set_module="%s/atom_vectors" % modules_dir
sys.path.append(set_module)
#
print("Current system path:")
print(sys.path)
#
# Import our own modules
#
from atom_vectors import *
#
# Read in the structure from Alice's bulk optimised car file
#
geom_file=sys.argv[1]
print("Processing file : %s" % geom_file)
#
atoms= read(geom_file)
natoms=len(atoms)
#
zfix_string = "17.0"
#
fix_these=[iii for iii in range(71,80)]
fix_these.append(1)
fix_these.append(3)
fix_these.append(2)
fix_map=zfix(atoms, zfix_string, fix_these)
need_zfix = np.sum(fix_map) > 0
#
for iatom in range(0,len(atoms)):
   if fix_map[iatom] == 1:
        print(f"%d %s %10.6f %10.6f %10.6f fixed"                        \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))
   else:
        print(f"%d %s %10.6f %10.6f %10.6f free"                         \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
write_vasp("POSCAR_fixed", atoms, vasp5=True)

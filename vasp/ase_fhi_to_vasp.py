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

fhi_file= sys.argv[1]                                           # the fhi file to convert
print("Converting file %s" % fhi_file)
#
atoms=read(fhi_file)
write('POSCAR_latest', atoms ,format='vasp')



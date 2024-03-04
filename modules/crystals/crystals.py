import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
#from ase.calculators.aims import Aims
from ase.calculators.vasp import Vasp2
from ase.optimize import BFGS
from ase.units import kJ
from ase.eos import EquationOfState
from ase.io import aims
#
# Code to do cell volume scaling for Murnighan plot.
# The script expects the structure to be in a POSCAR file in the 
# structure sub-directory.
#
import sys
import csv
import os
#
from ase import Atom
from ase import Atoms
#
# Define a matrix following Voigt convention
#
def voigt_mat(eee):
    vmat=[[eee[0], eee[5]/2.0,  eee[4]/2.0 ], \
          [eee[5]/2.0, eee[1],  eee[3]/2.0 ], \
          [eee[4]/2.0, eee[3]/2.0,  eee[2] ] ]

    return vmat
#
# Custom cell scaling
#
def cell_scale(atoms, eee):
#
# If eee is a number just uniformly scale
# If eee is a six membered list form voight matrix directly
# otherwise cause error
#
   this_type=str(type(eee))
   print("eee is of type: ", this_type)
   if 'list' in this_type and len(eee) == 6:
      rrr = eee
#
   elif "float" in this_type:
      rrr = [ eee, eee, eee, 0.0, 0.0, 0.0]
#
   elif type(eee) is "int":
      rrr = [ float(eee), float(eee), float(eee), 0.0, 0.0, 0.0]
#
   else:
      print("ERROR: Unknown directive for cell_scale")
      print("ERROR: cell_scale(atoms, eee) expects atoms structure and")
      print("ERROR: eee as either a six element strain matrix or")
      print("ERROR: a single value for uniform volume scaling of cell.")
      exit(0)
#
# Need a unit matrix to use when operating on cell vectors
#
#   unit_mat = np.zeros((3,3))
#   unit_mat[0][0]=1.0
#   unit_mat[1][1]=1.0
#   unit_mat[2][2]=1.0

   vmat = voigt_mat(rrr)
#
   print("matrix to be used for cell change")
   print(vmat)
#
   cell=atoms.get_cell()
   new_cell =  np.dot(vmat, cell)
#
   print("Original cell:")
   print(cell)
   print("New cell:")
   print(new_cell)
#
   atoms.set_cell(new_cell, scale_atoms=True)
#
   lengths= atoms.cell.lengths()
   angles= atoms.cell.angles()
   print("new cell lengths: ", lengths)
   print("new cell angles : ", angles)
#
   return atoms

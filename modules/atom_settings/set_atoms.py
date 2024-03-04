#
# functions that set the properties for atom lists
#
import numpy as np
from ase import Atom
from ase import Atoms
#
#
def atom_formal_charges(atoms):
#
   num_set = 0
#
# Set formal atom charges
#
   qZn =  2.0
   qO  = -2.0
   qPd =  0.0
   qH  =  1.0
   q_list = []
#
   totq = 0
   for iatom in range(0,len(atoms)):
      if ( atoms.symbols[iatom] == 'Zn' ):
         q_list.append(qZn)
         totq += qZn
      elif ( atoms.symbols[iatom] == 'O' ):
         q_list.append(qO)
         totq += qO
      elif ( atoms.symbols[iatom] == 'Pd' ):
         q_list.append(qPd)
         totq += qPd
      elif ( atoms.symbols[iatom] == 'H' ):
         q_list.append(qH)
         totq += qH
      else:
         print(f"Warning: Atom %d with symbol %s could not be assigned a charge." % (iatom, atoms.symbols[iatom]))
         print(f"Warning: charge set to zero.")
         q_list.append(0)

#
   atoms.set_initial_charges(q_list)


   return totq    


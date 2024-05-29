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
   qSn =  2.0
   qPb =  2.0
   qZn =  2.0
   qO  = -2.0
   qPd =  0.0
   qH  =  1.0
   q_list = []
#
# halides
   qF  = -1.0
   qCl = -1.0
   qBr = -1.0
   qI  = -1.0
#
# Arrange so two letter symbols come first 
#
   totq = 0
   for iatom in range(0,len(atoms)):
      if ( atoms.symbols[iatom] == 'Zn' ):
         q_list.append(qZn)
         totq += qZn
      elif ( atoms.symbols[iatom] == 'Pd' ):
         q_list.append(qPd)
         totq += qPd
      elif ( atoms.symbols[iatom] == 'Sn' ):
         q_list.append(qSn)
         totq += qSn
      elif ( atoms.symbols[iatom] == 'Pb' ):
         q_list.append(qPb)
         totq += qPb
      elif ( atoms.symbols[iatom] == 'Cl' ):
         q_list.append(qCl)
         totq += qCl
      elif ( atoms.symbols[iatom] == 'Br' ):
         q_list.append(qBr)
         totq += qBr
      elif ( atoms.symbols[iatom] == 'O' ):
         q_list.append(qO)
         totq += qO
      elif ( atoms.symbols[iatom] == 'H' ):
         q_list.append(qH)
         totq += qH
      elif ( atoms.symbols[iatom] == 'F' ):
         q_list.append(qF)
         totq += qF
      elif ( atoms.symbols[iatom] == 'I' ):
         q_list.append(qI)
         totq += qI
      else:
         print(f"Warning: Atom %d with symbol %s could not be assigned a charge." % (iatom, atoms.symbols[iatom]))
         print(f"Warning: charge set to zero.")
         q_list.append(0)
#
   atoms.set_initial_charges(q_list)

   return totq    


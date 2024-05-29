#
# Routines to carry out analysis of atoms lists, such as working out chemical formulae etc
#
# composition : returns the chemical formula of an atom set, using only those listed in ilist,
#               if a list of indices is supplied
#
import numpy as np
#
from ase import Atom
from ase import Atoms
from ase import neighborlist
#
# Return the composition of a set of atoms in a list using just particular members
# if specified.
# if ilist is the word "all" consider all atoms in the reference list.
#
def composition(atoms, ilist='all'):
#
   natoms=len(atoms)
#
   if natoms < 0:
     print("ERROR: empty atoms structure sent to composition function") 
     exit(0)
#
   if ( ilist == "all" ):
       ilist=[ i for i in range(0,natoms) ]
#
# Obtain types of atoms present
#
   types=[]
   num=[]
   for iatom in ilist:
      got=False
      for ityp in range(0,len(types)):
#         print("Comparing ",atoms.symbols[iatom],"to type:",types[ityp])
         if atoms.symbols[iatom]==types[ityp]:
            got=True
            this_ityp = ityp
#            print("Matched")
#
      if got:
         num[this_ityp]+=1
      else:
         types.append(atoms.symbols[iatom])
         num.append(1)
#      print("Now types: ", types, " num: ", num)
#      print(" ")
#
# Bring results into a single list
#
   comp = []
   for ityp in range(0,len(types)):
      entry=[types[ityp],num[ityp]]
#      print("entry:",entry)
      comp.append(entry)
#
   return comp      

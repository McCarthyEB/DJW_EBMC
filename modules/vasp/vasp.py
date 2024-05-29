#
# Routines for setting up vasp calculations
#
import numpy as np
#
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength


# Check if we have a number for cores_per_task and use to estimate ncore
#
def best_ncore(cores_per_task):
#
    core_diff=[]
    cores_list=[]
    test=cores_per_task.isdigit()
    if test:
        tot_cores=float(cores_per_task)
        ncore=int(np.sqrt(tot_cores))
        div = tot_cores/ncore
        print("cores_per_task= ", cores_per_task, " ideal ncore = %d" % ncore)
        print("cores_per_task/ncore = ", div)
#
        if tot_cores%ncore == 0:
           print("Already integer so can use ncore = %d" % ncore)
           return ncore
        else:
           print("This needs to be integer....")
#
        for n in range(2,int(np.sqrt(tot_cores))):
            print("n: ", n, " cores_per_task/n = ", tot_cores/n)
            if tot_cores%n == 0:
              m = int(tot_cores/n)
              cores_list.append(n)
              cores_list.append(m)
              core_diff.append(abs(ncore-n))
              core_diff.append(abs(ncore-m))
#
        print(cores_list)
        print(core_diff)
#
        ibest=-1
        diff_best=-1
        for indx in range(0,len(core_diff)):
           if indx==0:
              ibest=indx
              diff_best=abs(core_diff[indx])
           else:
              if abs(core_diff[indx]) < diff_best:
                 diff_best=abs(core_diff[indx])
                 ibest=indx
#
        if ibest < 0:
          print("Cannot divide cores_per_task into factors, putting ncore = 1")
          ncore=1
        else:
          ncore=cores_list[ibest]
          print("Best to use: %d i.e. ncore= %d" % (ibest, ncore))
#
    return ncore
#
# Sort an atom list into species
# Expects the species to have been identified by a call to "composition" in 
# the atom_analysis module first. This supplies the "types" list.
#
def sort_by_species(atoms, comp):
#
# If there are atom constraints need to make sure these follow the sorting too
#
# Assume these include lower layers if fixed so only implement zfix if POSCAR all free
   have_cons=False
   cons_orig = (atoms.constraints).copy()
   del atoms.constraints
#
   if (cons_orig):
#
# Find the indices set in the POSCAR file as fixed
#
     have_cons=True
     sub1=str(cons_orig[0]).split('[')
     sub2=sub1[1].split(']')
     sub1=sub2[0].split(',')
     fixed_list=[ int(sub1[i]) for i in range(0,len(sub1)) ]
     new_fixed_list=[]
#
   indx=0
   temp_atoms=atoms.copy()
#
# refill the temp_atoms list with atoms in types order
#
   for ityp in range(0,len(comp)):
      type=comp[ityp][0]
#
      for iatom in range(0,len(atoms)):
         if atoms.symbols[iatom] == type:
            temp_atoms.positions[indx]=atoms.positions[iatom]
            temp_atoms.symbols[indx]=atoms.symbols[iatom]
            if have_cons and iatom in fixed_list:
               new_fixed_list.append(indx)
            indx+=1
#
   if have_cons:
      cons_new=FixAtoms(indices=new_fixed_list)
      temp_atoms.set_constraint(cons_new)
#
   atoms=temp_atoms.copy()
#
   return atoms
#
# convert a vasp4 format file to a vasp5 format file 
# assuming the atom label line is in the title of the vasp4 file
#
def fix_poscar(filename):
#
   fptr=open(filename,'r')
   lines = fptr.readlines()
   fptr.close()
#
   nlines = len(lines)
#
   fptr=open(filename,'w')
#
   for iline in range(0,nlines):
      fptr.write(lines[iline])
      if iline == 4:
         fptr.write(lines[0])

   fptr.close()
#
   return

   
#





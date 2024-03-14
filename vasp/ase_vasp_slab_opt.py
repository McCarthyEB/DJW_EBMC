# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:55:21 2020

@author: dave
"""
import numpy as np
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

from ase.calculators.vasp import Vasp2
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
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
# this ase_vasp_opt.py script expects:
# the first arguement to be the directory where we will run vasp, 
# the second the number of cores available or a junk string if this is not required
# The third to be the JOB_ID, a unique identifier to add to output information
# The fourth is the path to the directory from which the job was launched
# The fifth is the sub-directory for the structure for this particular run.
# The sixth is the sub-directory for the sub_directory for this type of job, under the structure
# The seventh to be the name of the machine we are running on.
#           Currently this can be "hawk","thomas" or  "archer".
#
#
workdir = sys.argv[1]                                           # Directory where we will do the calculations
cores_per_task = sys.argv[2]                                    # cores assigned to tasks where available
jobid = sys.argv[3]                                             # Unique identifier for this job
stdout_file= "%s/vasp_%s_%s.stdout" % (sys.argv[4], sys.argv[5], jobid) # stdout file for VASP will appear in your launch directory
results_dir= "%s/%s/%s" % (sys.argv[4], sys.argv[5], sys.argv[6]) # directory for results, this should be on your $HOME space
traj_file = "%s/%s_%s_spin1.traj" % ( results_dir, sys.argv[5], jobid )      # File for trajectory output
traj_spin_file = "%s/%s_%s_spin2.traj" % ( results_dir, sys.argv[5], jobid ) # File for trajectory output
machine = sys.argv[7].lower()                                     # string defining machine
#
print(f"------------------------------------------------------------")
print(f"Work directory               : %s" % workdir )
print(f"stdout file for vasp         : %s" % stdout_file )
print(f"results_dir for this script  : %s" % results_dir )
print(f"trajectory file for this run : %s" % traj_file )
print(f"trajectory file for ISPIN 2  : %s" % traj_spin_file )
print(f"Machine running job          : %s" % machine )
print(f"------------------------------------------------------------")
#
# Check if debug mode has been requested by setting DEBUG in JOBID
if ( jobid.find("DEBUG") >= 0):
    debug=True
else:
    debug=False
#
if debug:
   print("This is a DEBUG test, VASP will not be run")
# Check if we have a float for z-fix and set if so
#
zfix_string=sys.argv[8]
test = zfix_string.replace('.','',1).isdigit()
if test:
   zfix = float(zfix_string)
   need_zfix = True
   print("Will fix atoms with z co-ordinate less that %10.6f" % zfix)
else:
   print("zfix not set so will not freeze any atom co-ordinates")
   need_zfix = False
#
# set command for running VASP, stdout file destination
# This is machine dependent, for a new machine add what would have been in
# the job control script to launch the calculation.
#
kpts_needed=(1, 1, 1)
ktot=0
for i in range(0,3):
     ktot+=kpts_needed[i]
#

if 'hawk' in machine:
  if ktot==3:
     cmd_string = "mpirun -n %s vasp_std > %s" % (cores_per_task, stdout_file)
  else:
     cmd_string = "mpirun -n %s vasp_std > %s" % (cores_per_task, stdout_file)
#
elif 'young' in machine:
  if ktot==3:
     cmd_string = "gerun vasp_gam > %s" % stdout_file
  else:
     cmd_string = "gerun vasp_std > %s" % stdout_file
#
elif 'archer' in machine:
  if ktot==3:
     cmd_string = "aprun -n $NPROC vasp_gam > %s" % stdout_file
  else:
     cmd_string = "aprun -n $NPROC vasp_std > %s" % stdout_file
else:
     print(f"ERROR: Unknown machine, do not know command for running VASP!")
     sys.exit("Python script exiting")
#
# Read in the structure from a POSCAR file
# This will have been copied from your structure directory to the
# working directory, it should be a VASP5 file, i.e. with the atom
# symbols above the line with number of ions.
#
#
# Read in the structure from a car file if this is a new run or a POSCAR file
# if a restart.
#
poscar_file=workdir + "/POSCAR"
car_file=workdir + "/struct.car"

if os.path.isfile(poscar_file):
   print("reading POSCAR file")
   atoms=read(poscar_file)
#
elif os.path.isfile(car_file):
   print("reading car file")
   atoms=read(car_file)
#
else:
   print("ERROR : No structure file, neither struct.car or POSCAR could be found")
   exit()
#
write('check_start.cif', atoms ,format='cif')
#
# preserve any constraints already set
# Assume these include lower layers if fixed so only implement zfix if POSCAR all free
have_cons=False
cons_orig = (atoms.constraints).copy()
#
# make list according to zfix if set
zfix_list=[]
if need_zfix:
  for iatom in range(0,len(atoms)):
     if (atoms.positions[iatom][2] < zfix):
         zfix_list.append(iatom)
         have_cons=True
  print("zfix_list :"), 
  print(zfix_list)
#
have_cons_orig=False
if (cons_orig):
#
# Find the indices set in the POSCAR file as fixed
#
  have_cons=True
  have_cons_orig=True
  sub1=str(cons_orig[0]).split('[')
  sub2=sub1[1].split(']')
  sub1=sub2[0].split(',')
  fixed_list=[ int(sub1[i]) for i in range(0,len(sub1)) ]
  print("Original fixed_list:" )
  print(fixed_list)
#
# Combine lists if you have both, if only zfix copy that to fixed_list
#
if (have_cons_orig and need_zfix):
  for ilist in range(0, len(zfix_list)):
    found_fix=False
    for jlist in range(0,len(fixed_list)):
       if ( zfix_list[ilist] == fixed_list[jlist] ):
          found_fix=True
    if ( not found_fix ):
       fixed_list.append(zfix_list[ilist])

  fixed_list.sort()
  print(f"Combined list: ")  
  print(fixed_list)
#
elif ( need_zfix ):
  fixed_list=zfix_list.copy()
#
fix_map=np.zeros(len(atoms))
if ( have_cons ):
  del atoms.constraints
  cons_orig=FixAtoms(indices=fixed_list)
  for icon in fixed_list:
      fix_map[icon]=1
  cons = FixAtoms(fixed_list) 
  atoms.set_constraint(cons)

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
# Define INCAR parameters
#
# Make self-consistent ground state
calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=kpts_needed)
# Use new PAW for Au
calc.set(setups={'Au': '_new'})
# Electronic structure settings
calc.set(prec='Accurate', ediff=1E-8, nelmdl=-4, encut=400, ispin=1, isym=0, ismear=0, sigma=0.2)
# Efficiency settings           
calc.set(algo='Fast', lreal=True ,addgrid=True)
# vdW settings
calc.set(ivdw= 11, vdw_s6=0.75)
#
# Calculate the initial energy, this will also establish the WAVECAR and CHGCAR etc files
#
calc.set(istart=0, prec='Low', ediff=1E-5, encut=400)
atoms.set_calculator(calc)
if (debug):
  e = 99.9
else:
  e=atoms.get_potential_energy()  # Run the calculation
print(f"Initial energy = %10.6f eV" % e )
#
# preserve any constraints already set, e.g. lower layers of slabs
cons_orig = (atoms.constraints).copy()
# constrain the bond defined by the index values supplied
#cons_new = FixBondLength( index_1, index_2 )
if not cons_orig:
    print(f"No original constraints set")
else:
    print(f"Constraints set %s" % cons_orig )

#    Relax the structure using slack options     
calc.set(istart=1)
opt = FIRE(atoms, trajectory=traj_file)

if (not debug):
   opt.run(fmax=0.01)

# Turn up accuracy then redo the energy
calc.set(istart=0, prec='Accurate', lreal='.FALSE.', ediff=1E-8, encut=500)
atoms.set_calculator(calc)
#
if (not debug):
   e=atoms.get_potential_energy()  # Run the calculation
print(f"Energy on raising accuracy = %10.6f eV" % e )

# Full accuracy for final stage
calc.set(istart=1)
opt = BFGS(atoms, trajectory=traj_file)
if (not debug):
   opt.run(fmax=0.01)
#
print(f"Optimisation reached required accuracy")
print(f"---------------------------------------")
print(f"Starting optimisation unrestricted.....")
# Turn on spin then redo the energy, forcing triplet state
calc.set(istart=0, ispin=2, nupdown=0)
atoms.set_calculator(calc)
if (not debug):
   e=atoms.get_potential_energy()  # Run the calculation
print(f"Energy unrestricted ( ISPIN 2 ) = %10.6f eV" % e )
# Full accuracy for final stage unrestricted
calc.set(istart=1)
opt = BFGS(atoms, trajectory=traj_spin_file)
#
if (not debug):
   opt.run(fmax=0.01)
#
print(f"Optimisation reached required accuracy unrestricted")


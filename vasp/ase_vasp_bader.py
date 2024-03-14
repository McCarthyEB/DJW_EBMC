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
elif 'thomas' in machine:
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
#
# Read in the structure from a POSCAR file
# This will have been copied from your structure directory to the
# working directory, it should be a VASP5 file, i.e. with the atom
# symbols above the line with number of ions.
#
atoms=read('POSCAR')
write('check_start.cif', atoms ,format='cif')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
# Define INCAR parameters : MAKE SURE THESE PARAMETERS ARE SET FOR YOUR SYSTEM REQUIREMENTS
#
# Make self-consistent ground state
calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=kpts_needed, directory=sys.argv[1])
# Use new PAW for Au
calc.set(setups={'Au': '_new'})
# Electronic structure settings
calc.set(prec='Accurate', ediff=1E-8, encut=500, ispin=2, isym=0, ismear=0, sigma=0.2)
# Efficiency settings           
calc.set(algo='Fast', lreal='.TRUE.',addgrid='.TRUE.')
# vdW settings
calc.set(ivdw= 11)
#
# Bader settings. NGXF etc sets the grid accuracy, the finer the more accurate. 
#
#calc.set(lcharg=True, laechg=True, ngx=400, ngy=400, ngz=400, ngxf=400, ngyf=400, ngzf=400 )
calc.set(lcharg=True, laechg=True)
#
# Calculate the energy, this will also establish the WAVECAR and CHGCAR etc files for the Bader analysis
#
atoms.set_calculator(calc)
e=atoms.get_potential_energy()  # Run the calculation
print(f"Calculated energy = %10.6f eV" % e )
print(f"Files should be written for Bader" )
#
# 
#


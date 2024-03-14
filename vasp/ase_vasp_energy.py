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
debug=True
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>6}: {arg}")
#
# this ase_vasp_opt.py script expects:
# the first arguement to be the directory where we will run vasp, 
# the second the number of cores available or a junk string if this is not required
# The third to be the JOB_ID, a unique identifier to add to output information
# The fourth is the path to the directory from which the job was launched
# The fifth is the sub-directory for the results from this particular run.
# The sixth to be the name of the machine we are running on.
#           Currently this can be "hawk","thomas" or  "archer".
#
#
mydir = sys.argv[1]                                             # Directory where we will do the calculations
cores_per_task = sys.argv[2]                                    # cores assigned to tasks where available
jobid = sys.argv[3]                                             # Unique identifier for this job
stdout_file= "%s/vasp_%s.stdout" % (sys.argv[4], jobid)         # stdout file for VASP will appear in your launch directory
results_dir= "%s/%s" % (sys.argv[4], sys.argv[5])               # directory for results, this should be on your $HOME space
machine = sys.argv[7].lower()                                   # string defining machine
#
# set command for running VASP, stdout file destination
# This is machine dependent
#
if 'hawk' in machine:
     cmd_string = "mpirun -n %s vasp_std > %s" % (cores_per_task, stdout_file)
elif 'thomas' in machine:
     cmd_string = "gerun vasp_gam > %s" % stdout_file
elif 'archer' in machine:
     cmd_string = "aprun -n $NPROC vasp_gam > %s" % stdout_file
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
traj_file = "%s/ase_job_%s.traj" % ( results_dir, jobid )
print(f"traj file  : %s" % traj_file)
#
# Define INCAR parameters
#
# Make self-consistent ground state
calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=(1, 1, 1), directory=sys.argv[1])
# Use new PAW for Au
calc.set(setups={'Au': '_new'})
# Electronic structure settings
calc.set(prec='Accurate', ediff=1E-8, nelmdl=-4, encut=500, ispin=1, isym=0, ismear=0, sigma=0.2)
# Efficiency settings           
calc.set(algo='Fast', lreal='.FALSE.',addgrid='.TRUE.')
# vdW settings
calc.set(ivdw= 11, vdw_s6=0.75)
#
# Calculate the energy
#
calc.set(istart=0, prec='Accurate', ediff=1E-8, encut=500, ispin=1, isym=0, ismear=0, sigma=0.1)
atoms.set_calculator(calc)
#
if debug:
   e=99.9
else:
   e=atoms.get_potential_energy()  # Run the calculation
#
print(f"Initial energy restricted ( ISPIN 1 ) = %10.6f eV" % e )
#
encut_list=[200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0]

for encut_val in encut_list:
# Turn up accuracy then redo the energy
   calc.set( ispin=2, encut=encut_val )
   atoms.set_calculator(calc)
   if debug:
      e=99.9
   else:
      e=atoms.get_potential_energy()  # Run the calculation

   print("At encut = ", encut_val, " energy = ", e, "eV")
   


















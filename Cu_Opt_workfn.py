import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp2
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.units import kJ
from ase.eos import EquationOfState
from pathlib import Path
#
# planewave cut-off scan to identify the best planewave cut-off sampling to use.
# The script expects the structure to be in a POSCAR file in the 
# structure sub-directory.
#
a = 3.6 # approximate lattice constant
b = a / 2

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

from ase.units import kJ
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
debug = False
#
if (debug):
  print("EMT run")
  csv_filename= "eos_Cu_EMT.csv" # filename for csv record 
#
else:
  print("VASP run")
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
  stdout_file= "%s/vasp_%s_%s.stdout" % (sys.argv[4], sys.argv[5], jobid) # stdout for VASP will appear in your launch directory
  results_dir= "%s/%s/%s" % (sys.argv[4], sys.argv[5], sys.argv[6]) # directory for results, this should be on your $HOME space
  csv_filename= "%s/eos_%s_%s.csv" % ( results_dir, sys.argv[5], jobid ) # filename for csv record 
  machine = sys.argv[7].lower()                                     # string defining machine
#
  print(f"------------------------------------------------------------")
  print(f"Work directory               : %s" % workdir )
  print(f"stdout file for vasp         : %s" % stdout_file )
  print(f"results_dir for this script  : %s" % results_dir )
  print(f"csv file for data            : %s" % csv_filename )
  print(f"Machine running job          : %s" % machine )
  print(f"------------------------------------------------------------")
#
# set command for running VASP, stdout file destination
# This is machine dependent, for a new machine add what would have been in
# the job control script to launch the calculation.
#
  kpts_needed=(3, 3, 1)
  ktot=0
  for i in range(0,3):
        ktot+=kpts_needed[i]
#
  if 'hawk' in machine:
     if ktot==3:
       cmd_string = "mpirun -n %s vasp_gam > %s" % (cores_per_task, stdout_file)
     else:
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

if ( not debug ):
#
# Define INCAR parameters
#
# Make self-consistent ground state
  calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=kpts_needed)
# Use new PAW for Au
  calc.set(setups={'Au': '_new'})
# Electronic structure settings, note that encut is the setting of the planewave cut-off in eV
  calc.set(prec='Accurate', ediff=1E-8, nelmdl=-4, encut=450, ispin=1, isym=0, ismear=0, sigma=0.2)
# Efficiency settings           
  calc.set(algo='Fast', lreal= True, addgrid= True )
# vdW settings
#calc.set(ivdw= 11, vdw_s6=0.75)
#
# Attached the calculator to the atoms
  atoms.set_calculator(calc)
else:
  atoms.set_calculator(EMT())
#
energy = atoms.get_potential_energy()
print("Starting energy : ", energy, " eV")
#
opt = BFGS(atoms, restart='qn.pckl')
opt.run(fmax=0.01)
#
print(f"Optimisation reached required accuracy")
#
if ( not debug ):
  calc.set(addgrid= False, lvhar = True,  lvtot = False ) 
  atoms.set_calculator(calc)
energy = atoms.get_potential_energy()

# Extract energies:

print("Optimised energy : ", energy, " eV")
#

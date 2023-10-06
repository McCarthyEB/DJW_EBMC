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
import time
#
from ase import Atom
from ase import Atoms

from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength

from ase.io import read, write, Trajectory
from ase.vibrations import Vibrations

#from ase.calculators.vasp import Vasp2
from ase.calculators.aims import Aims
from carmm.run.aims_path import set_aims_command

#
# Local definition of atom...atom distance measurement




PROG = "AIMS"
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
debug=False 
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
if debug:
   print("DEBUG DEBUG DEBUG DEBUG......")
   print("DEBUG mode for testing code, prog will not be run")
#
workdir = sys.argv[1]                                           # Directory where we will do the calculations
cores_per_task = sys.argv[2]                                    # cores assigned to tasks where available
jobid = sys.argv[3]                                             # Unique identifier for this job
stdout_file= "%s/vasp_%s_%s.stdout" % (sys.argv[4], sys.argv[5], jobid) # stdout file for VASP will appear in your launch directory
results_dir= "%s/%s/%s" % (sys.argv[4], sys.argv[5], sys.argv[6]) # directory for results, this should be on your $HOME space
traj_file = "%s/%s_%s.traj" % ( results_dir, sys.argv[5], jobid )      # File for trajectory output

if PROG == "VASP":
  poscar_filestem= "POSCAR_%s" % sys.argv[5]                        # stem to use for poscar files output at scan steps
if PROG == "AIMS":
  aims_filestem= "geometry_%s.in" % sys.argv[5]                        # stem to use for poscar files output at scan steps
csv_filename= "%s/dos_v_energy_%s_%s.csv" % ( results_dir, sys.argv[5], jobid ) # filename for csv record of DOS 
movie_filename= "%s/movie_%s_%s.xyz" % ( results_dir, sys.argv[5], jobid )        # filename for xyz movie of bond step 
machine = sys.argv[7].lower()                                     # string defining machine
#
# Check if we have a float for z-fix and set if so
#
print(f"------------------------------------------------------------")
print(f"Work directory                 : %s" % workdir )
print(f"stdout file for prog           : %s" % stdout_file )
print(f"results_dir for this script    : %s" % results_dir )
print(f"trajectory file for this run   : %s" % traj_file )
print(f"Stem for POSCAR files of steps : %s" % traj_file )
print(f"Machine running job            : %s" % machine )
#
zfix_string=sys.argv[8]
test = zfix_string.replace('.','',1).isdigit()
if test:
   zfix = float(zfix_string)
   need_zfix = True
   print("\nWill fix atoms with z co-ordinate less that %10.6f" % zfix)
else:
   print("\nzfix not set so will not freeze any atom co-ordinates unless pre-defined in geometry file.")
   need_zfix = False
#
# Get indices for fixed atoms
# The bond should move from 1 to 2 so that index_2 is for the atom that
# will be moved.
#hess_file=sys.argv[9]
print(f"------------------------------------------------------------")
#
# set command for running VASP, stdout file destination
# This is machine dependent, for a new machine add what would have been in
# the job control script to launch the calculation.
#
if PROG == "VASP":

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
  atoms=read('POSCAR')

if PROG == "AIMS":

  set_aims_command(hpc="hawk", basis_set="light", defaults=2020, nodes_per_instance=sys.argv[2])
#
  print("Running command set to:")
  print(os.environ['ASE_AIMS_COMMAND'])
  atoms=read('geometry.in')
#  if os.path.exists('hessian.aims'):
 #   hessian_data = ('hessian.aims')

write('check_start.cif', atoms ,format='cif')
print("Cell:",atoms.cell)
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))


#
#
if PROG == "VASP":
  # Define INCAR parameters
  #
  # Make self-consistent ground state
  calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=(1, 1, 1))
  # Use new PAW for Au
  calc.set(setups={'Au': '_new'})
  # Electronic structure settings
  calc.set(istart=0, prec='Accurate', ediff=1E-8, encut=500, ispin=1, isym=0, ismear=0, sigma=0.2)
  # Efficiency settings
  calc.set(algo='Fast', lreal='.FALSE.',addgrid='.TRUE.')
  # vdW settings
  calc.set(ivdw= 11, vdw_s6=0.75)
  
  atoms.set_calculator(calc)

if PROG == "AIMS":
#Define control.in parameters
  
# Define control.in parameters
#
# New method that gives a default calculator
# Use dimensions 2 for slab calculations, this will turn on dipole correction
  from carmm.run.aims_calculator import get_aims_and_sockets_calculator

  
  sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=2,
                                    xc_pre=['pbe', '15'],
                                    xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
                                    spin='none',
                                    k_grid=(3, 3, 1),
                                    relativistic=('atomic_zora', 'scalar'),
                                    compute_forces=".true.",
                                    sc_accuracy_etot=1e-3,
                                    sc_accuracy_eev=1e-3,
                                    sc_accuracy_rho=1e-6,
                                    sc_accuracy_forces=1e-4  )

  atoms.set_calculator(sockets_calc)
#
# get initial energy
#
if debug:
   e=9999.0
else:
   #atoms.set_calculator(calc)
   energy=atoms.get_potential_energy()
   # Report the data
   print("Initial energy : ", energy, " eV")
   #
   dyn = BFGS(atoms, trajectory='ZnO.traj', restart='qn.pckl')  #Put optimisation back in
   opt.run(fmax=0.01)
   energy=atoms.get_potential_energy()
   
   write('geometry.in.next_step', atoms, format='aims')
   
#   vib = Vibrations(atoms)
#   vib.run()
#   
#   #get the hessian matrix
#   
#   hessian=vib.get_hessian()
#   
#   
#   print('Getting density of states...')
#   dos=DOS(calc, width=0.2) 
#   d = dos.get_dos()
#   e = dos.get_energies()
#   #want xvsy -> energy, then dos 
#   #dos_data = []
#   csv_headers = ['Energy (eV)', 'Density of States']
#   csv_row = [e[istep],d[istep]]
#
#   for istep in range(len(e)):
#   #    dos_data.append((e[i], d[i]))
#        if istep == 0:
#          csv_file = open(csv_filename, 'w', newline='')
#          csv_writer = csv.writer(csv_file)
#          csv_writer.writerow(csv_headers)
#          csv_writer.writerow(csv_row)
#          csv_file.close()
#  
#        else: 
#          csv_file = open(csv_filename, 'a', newline='')
#          csv_writer = csv.writer(csv_file)
#          csv_writer.writerow(csv_row)
#          csv_file.close()

   print("Calculation finished!")
#


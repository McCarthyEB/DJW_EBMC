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
csv_filename= "%s/energy_vs_dist_%s_%s.csv" % ( results_dir, sys.argv[5], jobid ) # filename for csv record of bond step 
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
ind1_string=sys.argv[9]
ind2_string=sys.argv[10]
# Define the bond that you will change
#
test1=ind1_string.isdigit()
test2=ind2_string.isdigit()
if not test1:
  print("ERROR: No index for atom 1 in a scan job, need to define in script")
if not test2:
  print("ERROR: No index for atom 2 in a scan job, need to define in script")
#
index_1 = int(ind1_string)
index_2 = int(ind2_string) 
#
print("Scan job moving atom %d with reference to atom %d." % (index_2, index_1))
#
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
  print("Running command set to:")
  print(os.environ['ASE_AIMS_COMMAND'])

  atoms=read('geometry.in')


write('check_start.cif', atoms ,format='cif')
print("Cell:",atoms.cell)
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))


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
# Find the indices set in the geometry file as fixed
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
  print(f"Combined original and Z-fix lists")  
#
elif ( need_zfix ):
  print(f"No original constraints so copying Z-fix list to fixed_list")  
  fixed_list=zfix_list.copy()
else:
  print("Zfix not requested so sticking with original list")
#
print("Fixed list now:")
print(fixed_list)
#
fix_map=np.zeros(len(atoms))
if ( have_cons ):
  del atoms.constraints
  print("Re-setting cons_orig.....")
  cons_orig=FixAtoms(indices=fixed_list)
  for icon in fixed_list:
      fix_map[icon]=1
#
print("cons_orig:")
print(cons_orig)
#
start_dist=atom_dist(atoms, index_1, index_2)
#
# Step size in Angstroms and define number of steps to make
step1 = 0.1 
step2 = 0.1
num_steps = 12

for iatom in range(0,len(atoms)):
   if fix_map[iatom] == 1:
        print(f"%d %s %10.6f %10.6f %10.6f fixed"                        \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))
   else:
        if iatom == index_1:
           print(f"%d %s %10.6f %10.6f %10.6f free scan atom ref"      \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))

        elif  iatom == index_2:
           print(f"%d %s %10.6f %10.6f %10.6f free scan atom moving"                         \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))

        else:
           print(f"%d %s %10.6f %10.6f %10.6f free"                         \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2]))

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



if PROG == "AIMS":
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

#
# get initial energy
#
if debug:
   e=9999.0
else:
   atoms.set_calculator(sockets_calc)
   e=atoms.get_potential_energy()  # Run the calculation
#
print(f"Initial energy = %10.6f eV" % e )
#
# preserve any constraints already set, e.g. lower layers of slabs
#`cons_orig = (atoms.constraints).copy()
# constrain the bond defined by the index values supplied
#cons_new = FixBondLength( index_1, index_2 )
#if not have_cons:
#    print(f"No original constraints set")
#
#
dist=[]
energy=[]
new_dist=start_dist
csv_headers=['Step', 'distance / Angs', 'Energy / eV']
# step through geometries
for istep in range(0,num_steps):
#
    print(f"Step %d, for this pass distance set to %10.6f\n" % (istep,new_dist))
    del atoms.constraints
    atoms.set_distance(index_1, index_2, new_dist, fix=0)

    if istep < 4:
      new_dist = new_dist+step1
    else:
      new_dist = new_dist+step2
#
    if not have_cons:
# constrain the bond defined by the index values supplied
       print(f"Constraints set, just new")
       cons_new = FixBondLength( index_1, index_2 )
       atoms.set_constraint(cons_new)
    else:
       print(f"Constraints set, new and old")
       cons_new = FixBondLength( index_1, index_2 )
#       atoms.set_constraint([cons_orig[0], cons_new])
       atoms.set_constraint([cons_orig, cons_new])

    print(f"Constraints set: %s" % atoms.constraints)
#
#
    if debug:
       energy.append(float(istep))
    else:
#    Relax the structure under the set constraint
       opt = FIRE(atoms, trajectory=traj_file, restart='qn.pckl')
       opt.run(fmax=0.05)
       energy.append(atoms.get_potential_energy())

    print(f"Optimisation reached required accuracy")
#
# Record the distance and energy
    dist.append(atom_dist(atoms, index_1, index_2))

    print(f"step: %d dist: %10.6f Energy: %10.6f"   \
          % ( istep, dist[istep], energy[istep] ) )

    csv_row = [istep, dist[istep], energy[istep]]    

    if PROG == "AIMS":
    
      aims_filename= "%s_%d" % ( aims_filestem, istep ) 
      write(aims_filename, atoms ,format='vasp')

    if PROG == "VASP":

# Write POSCAR for this step tagged with job and step number
      poscar_filename= "%s_%d" % ( poscar_filestem, istep ) 
      write(poscar_filename, atoms ,format='vasp')
    
    if istep == 0:
# Update movie
        write(movie_filename, atoms ,format='xyz')
# Update csv record of energy
        csv_file = open(csv_filename, 'w', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_headers)
        csv_writer.writerow(csv_row)
        
        csv_file.close()

    else:
        write(movie_filename, atoms ,format='xyz', append=True)
# Update csv record of energy
        csv_file = open(csv_filename, 'a', newline='')
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(csv_row)
        csv_file.close()
      
print(f"")
print(f"Run completed.......")
print(f"")


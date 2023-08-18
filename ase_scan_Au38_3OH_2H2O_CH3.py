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
code_check=False

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
traj_file = "%s/%s_%s.traj" % ( results_dir, sys.argv[5], jobid )      # File for trajectory output
poscar_filestem= "POSCAR_%s" % sys.argv[5]                        # stem to use for poscar files output at scan steps
csv_filename= "%s/energy_vs_dist_%s_%s.csv" % ( results_dir, sys.argv[5], jobid ) # filename for csv record of bond step 
movie_filename= "%s/movie_%s_%s.xyz" % ( results_dir, sys.argv[5], jobid )        # filename for xyz movie of bond step 
machine = sys.argv[7].lower()                                     # string defining machine
#
print(f"------------------------------------------------------------")
print(f"Work directory                 : %s" % workdir )
print(f"stdout file for vasp           : %s" % stdout_file )
print(f"results_dir for this script    : %s" % results_dir )
print(f"trajectory file for this run   : %s" % traj_file )
print(f"Stem for POSCAR files of steps : %s" % traj_file )
print(f"Machine running job            : %s" % machine )
print(f"------------------------------------------------------------")
#
if code_check:
   print(f"THIS IS A CODE CHECK RUN,")
   print(f"CHECKING SYNTAX WITHOUT HEAVY COMPUTATION")
   print(f"------------------------------------------------------------")
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
# Define INCAR parameters
#
# Make self-consistent ground state
calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=(1, 1, 1), directory=sys.argv[1])
# Use new PAW for Au
calc.set(setups={'Au': '_new'})
# Electronic structure settings
calc.set(prec='Accurate', ediff=1E-8, nelmdl=-4, encut=500, ispin=1, isym=0, ismear=0, sigma=0.2)
# Efficiency settings           
calc.set(algo='Fast', lreal='.TRUE.',addgrid='.TRUE.')
# vdW settings
calc.set(ivdw= 11, vdw_s6=0.75)
#
# Define the bond that you will change
# The bond should move from 1 to 2 so that index_2 is for the atom that
# will be moved.
index_1 = 48
index_2 = 53
start_dist=atom_dist(atoms, index_1, index_2)
#
# Step size in Angstroms and define number of steps to make
step1 = -0.25
step2 = -0.05
num_steps = 22
#
# For new jobs set start_step = 0
# For continuation set at the final, unconverged step number from last scan
#
start_step = 11 
switch = 14

# Lower accuracy for scan step optimisations
calc.set(prec='Low', ediff=1E-5, encut=400)
#
# get initial energy
#
atoms.set_calculator(calc)

if code_check:
   e = 999.9
else:
   e=atoms.get_potential_energy()  # Run the calculation

print(f"Initial energy LREAL TRUE = %10.6f eV" % e )
#
# preserve any constraints already set, e.g. lower layers of slabs
cons_orig = (atoms.constraints).copy()
# constrain the bond defined by the index values supplied
cons_new = FixBondLength( index_1, index_2 )
if not cons_orig:
    print(f"No original constraints set")
#
#
dist=[]
energy=[]
#
# For continuation pad out start of scan with zeroes
#
csv_headers=['Step', 'distance / Angs', 'Energy / eV']
#
if start_step > 0:
#
   csv_file = open(csv_filename, 'w', newline='')
#
   for istep in range(0, start_step):
     dist.append(0.0)
     energy.append(0.0) 
     csv_row = [istep, dist[istep], energy[istep]]    
    
# Update csv record of energy
     csv_writer = csv.writer(csv_file)
     if istep == 0:
        csv_writer.writerow(csv_headers)
#
     csv_writer.writerow(csv_row)
        
   csv_file.close()
#
# Get ready for the scan loop
#
new_dist=start_dist
# step through geometries
for istep in range(start_step, num_steps):
#
    print(f"Step %d, for this pass distance set to %10.6f\n" % (istep,new_dist))
    del atoms.constraints
    atoms.set_distance(index_1, index_2, new_dist, fix=0)

    if istep < switch:
      new_dist = new_dist+step1
    else:
      new_dist = new_dist+step2
#
    if not cons_orig:
# constrain the bond defined by the index values supplied
       print(f"Constraints set, just new")
       cons_new = FixBondLength( index_1, index_2 )
       atoms.set_constraint(cons_new)
    else:
       print(f"Constraints set, new and old")
       cons_new = FixBondLength( index_1, index_2 )
       atoms.set_constraint([cons_orig[0], cons_new])

    print(f"Constraints set: %s" % atoms.constraints)
#
#    Relax the structure under the set constraint
    if not code_check:
       opt = FIRE(atoms, trajectory=traj_file)
       opt.run(fmax=0.05)

    print(f"Optimisation reached required accuracy")
#
# Record the distance and energy
    dist.append(atom_dist(atoms, index_1, index_2))
#
    if not code_check:
       energy.append(atoms.get_potential_energy())
    else:
       energy.append(10.0*istep)

    print(f"step: %d dist: %10.6f Energy: %10.6f"   \
          % ( istep, dist[istep], energy[istep] ) )

    csv_row = [istep, dist[istep], energy[istep]]    

# Write POSCAR for this step tagged with job and step number
    poscar_filename= "%s_%d" % ( poscar_filestem, istep ) 
    write(poscar_filename, atoms ,format='vasp')
    
    if istep == start_step:
# Update movie
        write(movie_filename, atoms ,format='xyz')
# Update csv record of energy
        if start_step == 0:
           csv_file = open(csv_filename, 'w', newline='')
           csv_writer = csv.writer(csv_file)
           csv_writer.writerow(csv_headers)
        else:
           csv_file = open(csv_filename, 'a', newline='')
           csv_writer = csv.writer(csv_file)
#
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


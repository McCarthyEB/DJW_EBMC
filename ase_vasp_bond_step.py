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
        print(f"Argument {i:>6}: {arg}")
#
# this ase_vasp_opt.py script expects:
# the first arguement to be the directory where we will run vasp, 
# the second the number of cores available or a junk string if this is not required
# The third to be the JOB_ID, a unique identifier to add to output information
# The fourth is the path to the directory from which the job was launched
# The fifth is the sub-directory for the results from this particular run.
# The sixth to be the name of the machine we are running on.
#           Currently this can be "hawk" or "thomas".
#
#
# set command for running VASP, stdout file destination
# directory for results back on main
#
mydir = sys.argv[1]                                             # Directory where we will do the calculations
cores_per_task = sys.argv[2]                                    # cores assigned to tasks where available
jobid = sys.argv[3]                                             # Unique identifier for this job
stdout_file= "%s/vasp_%s.stdout" % (sys.argv[4], jobid)         # stdout file for VASP will appear in your launch directory
results_dir= "%s/%s" % (sys.argv[4], sys.argv[5])               # directory for results, this should be on your $HOME space
csv_filename= "%s/energy_vs_dist_%s_%s.csv" % ( results_dir, sys.argv[5], jobid ) # filename for csv record of bond step 
movie_filename= "%s/movie_%s_%s.xyz" % ( results_dir, sys.argv[5], jobid )        # filename for xyz movie of bond step 
poscar_filestem= "%s/POSCAR_%s_" % ( results_dir, jobid )                          # filestem for POSCAR of bond step 
machine = sys.argv[6].lower()                                                     # string defining machine

if 'hawk' in machine:
     cmd_string = "mpirun -n %s vasp_std > %s" % (cores_per_task, stdout_file)
elif 'thomas' in machine:
     cmd_string = "gerun vasp_gam > %s" % stdout_file
elif 'archer' in machine:
     cmd_string = "aprun -n $NPROC vasp_std > %s" % stdout_file
else:
     print(f"ERROR: Unknown machine, do not know command for running VASP!")
     

atoms=read('POSCAR')
write('check_start.cif', atoms ,format='cif')

for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
# Define the bond that you will change
# The bond should move from 1 to 2 so that index_2 is for the atom that
# will be moved.
index_1 = 44
index_2 = 40
start_dist=atom_dist(atoms, index_1, index_2)
#
# Step size in Angstroms and define number of steps to make
step = -0.4
num_steps = 19
#
# Define INCAR parameters
#
# Make self-consistent ground state
calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=(1, 1, 1), directory=sys.argv[1])
# Use new PAW for Au            
calc.set(setups={'Au': '_new'})
# Electronic structure settings
calc.set(prec='Accurate', ediff=1E-7, encut=500, icharg=2, ispin=1, isym=0, ismear=0, sigma=0.02)
# Efficiency settings           
calc.set(algo='Fast', lreal='.FALSE.',addgrid='.TRUE.')
# vdW settings
calc.set(ivdw= 11, vdw_s6=0.75)
 
atoms.set_calculator(calc)
atoms.get_potential_energy()  # Run the calculation
#

traj_file = "%s/ase_job_%s.traj" % ( results_dir, jobid )
#print(f"traj file  : %s" % traj_file)

# Lower accuracy for scan step optimisations
calc.set(prec='Low', ediff=1E-5, encut=400, ispin=1, isym=0, ismear=0, sigma=0.02)

# preserve any constraints already set, e.g. lower layers of slabs
cons_orig = (atoms.constraints).copy()
# constrain the bond defined by the index values supplied
cons_new = FixBondLength( index_1, index_2 )
if not cons_orig:
    print(f"No original constraints set")

# Full accuracy for final stage
#calc.set(prec='Accurate', ediff=1E-8, encut=500, ispin=1, isym=0, ismear=0, sigma=0.02)
#opt = BFGS(atoms, trajectory=traj_file)
#opt.run(fmax=0.01)

dist=[]
energy=[]
new_dist=start_dist
csv_headers=['Step', 'distance / Angs', 'Energy / eV']
# step through geometries
for istep in range(0,num_steps):
#
    del atoms.constraints
    atoms.set_distance(index_1, index_2, new_dist, fix=0)
    print(f"Making step %10.6f for this pass distance set to %10.6f\n" % (step,new_dist))

    new_dist = new_dist+step
#
# refine step as end point gets closer
#
#   if istep == 6: 
#       step = step/2
#   if istep == 14:
#       step = step/2
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
    opt = FIRE(atoms, trajectory=traj_file)
    opt.run(fmax=0.1)

    print(f"Optimisation reached required accuracy")
#    relax = QuasiNewton(atoms)
#    relax.run(fmax=0.05)
#
# Record the distance and energy
    dist.append(atom_dist(atoms, index_1, index_2))
    energy.append(atoms.get_potential_energy())

    print(f"step: %d dist: %10.6f Energy: %10.6f"   \
          % ( istep, dist[istep], energy[istep] ) )

    csv_row = [istep, dist[istep], energy[istep]]    

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



import numpy as np
from scipy import sparse
import sys
import os
# Function to obtain HOME
from os.path import expanduser

#from ase import Atom
from ase import Atoms
from ase.io import read
from ase import neighborlist
#
from ase.io import read, write, Trajectory
#
from ase.optimize import BFGS
from ase.constraints import FixInternals
from ase.calculators.orca import ORCA
import re

import ase.io.orca as io

#from ase.calculators.genericfileio import (BaseProfile, CalculatorTemplate,
 #                                          GenericFileIOCalculator)


from ase.calculators.orca import OrcaProfile
#
sockets=False

if sockets == True:
  from carmm.run.aims_path import set_aims_command
  from carmm.run.aims_calculator import get_aims_and_sockets_calculator
#
else:
  from ase.calculators.aims import Aims
# Make modules available
#
home = expanduser("~")+"/python/DJW_EBMC"
set_module="%s/modules/data_io" % home
sys.path.append(set_module)
set_module="%s/modules/slab" % home
sys.path.append(set_module)
set_module="%s/modules/atom_settings" % home
sys.path.append(set_module)
set_module="%s/modules/num_digi" % home
sys.path.append(set_module)
set_module="%s/modules/vectors" % home
sys.path.append(set_module)
set_module="%s/modules/atom_vectors" % home
sys.path.append(set_module)
#set_module="%s/modules/orca_profile" % home
#sys.path.append(set_module)
#
print("Current system path:")
print(sys.path)
#
# Import our own modules
#
from atom_vectors import *
#from orca_profile import *
#
print("ASE script starting....")
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
PROG="ORCA"
# 
if sockets == True:
  if PROG == "AIMS":
  
    set_aims_command(hpc="hawk", basis_set="light", defaults=2020, nodes_per_instance=sys.argv[2])
  #
    print("Running command set to:")
    print(os.environ['ASE_AIMS_COMMAND'])
#
# Read in the list of atoms from the structure file into the ase Atoms object.
# the loop prints out symbols and co-ordinates, note how to get these sub-structures
# out of the atom object.
#
# Passed variables:
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
hess_file=sys.argv[9]
struct=sys.argv[10]
print(f"------------------------------------------------------------")
##
print(f"------------------------------------------------------------")
#
atoms=read(struct)
natoms=len(atoms)
#
atoms.set_cell([10.0,10.0,10.0])
atoms.set_pbc(pbc=True)
#
# Define control.in parameters
#
# New method that gives a default calculator
# Use dimensions 2 for slab calculations, this will turn on dipole correction
if PROG == "AIMS":
  if sockets == True:  
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=3,
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
  else:
    calc = Aims(xc_pre=['pbe', '10'], #we do here 10 steps with PBE to stabilize the system
           override_warning_libxc="true", #MBEEF is not impplemented directly with FHI, for that we need to call libxc library
#           force_correction="true",
           xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL', #MBEEF
           spin='none',
           k_grid=(4,4,1),
           relativistic=('atomic_zora', 'scalar'),
           compute_forces="true",
           sc_accuracy_etot=1e-3,
           sc_accuracy_eev=1e-3,
           sc_accuracy_rho=1e-6,
           sc_accuracy_forces=1e-6,
           )
    print("Hartree potential force correction turned off")
    atoms.set_calculator(calc)

if PROG == "ORCA":
   MyOrcaProfile = OrcaProfile(command="/apps/chemistry/orca/5.0.0/el7/orca")
   calc = ORCA(profile=MyOrcaProfile, orcasimpleinput='wB97X-V def2-TZVP Opt', orcablocks='%pal nprocs 40 end') 
   #atoms.set_calculator(calc)
   atoms.calc = calc
#
   test=atoms.get_potential_energy()
   print("Test of calculator: initial energy is", test)
   opt = BFGS(atoms, trajectory='rot_check.traj')  #Put optimisation back in
   opt.run(fmax=0.01)
# The find_mols routine identifies atoms that are chemically bonded to one another based on atom-atom
# distances. These sub-groups of atoms are called "molecules". The routine returns the number of molecules
# found "n_mols" and a two dimensional array, molinds[imol][iatom] gives the index of the iatom atom within
# molecule imol in the atoms list.
# neigh_matrix[i,j] is 1 if atom i is bonded to atom j and zero otherwise
#
n_mols, molinds, neigh_matrix = find_mols(atoms)
#
#  Generate a list of lists for each atom to hold its neighbour indexes
neigh_lists= gen_neigh_list(atoms, neigh_matrix)
#
# dihed_inds defines the bond around which to rotate
# the dihedral angle to be rotated are defined by the elements                
# 0--1--2--3 the fragment containing 2--3 will be moved
#
dihed_inds=[0,1,5,9]
#
# Check that the four atoms defining the dihedral are in the same molecule
#
found=0
mol_indx=[]
for imol in range(0, n_mols):
    for iatom in range(0,len(molinds[imol])):
        indx=molinds[imol][iatom]
        for iii in dihed_inds:
           if indx==iii:
              found=found+1
              mol_indx.append(imol)
#
print(found)
if found == 4:
   print("All dihedral atoms are in molecule %d..." % mol_indx[0])
else:
   print("ERROR: The four atoms defining the dihedral are not in the same molecule")
   print("ERROR: Molecular indices found as: ", mol_indx)
   exit(0)
#
# gather the indices of the two fragments on either side of the bond
#
frag1= mol_frag(atoms, mol_indx[0], dihed_inds[1],  dihed_inds[2], molinds, neigh_lists)
frag2= mol_frag(atoms, mol_indx[0], dihed_inds[2],  dihed_inds[1], molinds, neigh_lists)
#
print("Fragment 1: ", frag1)
print("Fragment 2: ", frag2)
#
# Rotate the fragment 2 around the dihedral by 90.0 degrees
#
flags = [ 0 for iatom in range(0, natoms) ]
for iatom in frag2:
   flags[iatom] = 1
#
# 
print("Dihedral defining atoms %d--%d--%d--%d" \
           % (dihed_inds[0], dihed_inds[1], dihed_inds[2], dihed_inds[3]))
print("Rotation mask: ", flags, "\n")
#
for imol in range(0, n_mols):
    print(f"Atoms in molecule %d" % imol)
    for iatom in molinds[imol]:
        nstr="Neighs: "
        for iii in neigh_lists[iatom]:
           nstr=nstr+"%d," % iii
#       
        print(f"%d %s %10.6f %10.6f %10.6f %s"                 \
            % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
                 atoms.positions[iatom][1], atoms.positions[iatom][2], nstr))
#
# Iterate the dihedral angle
#
dihedral_start = atoms.get_dihedral(dihed_inds[0], dihed_inds[1], dihed_inds[2], dihed_inds[3], mic=True)
#
if dihedral_start > 180.0:
   dihedral_start = dihedral_start - 360.0
#
dihedral_steps = [ i for i in range(0,80,20)]
dihedral_steps = dihedral_steps + [ i for i in range(80,120,5)]
dihedral_steps = dihedral_steps + [ i for i in range(120,200,20)]
#
print("\nInitial dihedral angle: %10.6f degrees.\n" % dihedral_start)
#
istep=0


for dihed in dihedral_steps:
#
   new_dihed = dihedral_start + dihed
   print("Dihedral step: %d, %10.6f degrees, actual angle: %10.6f\n" % \
                       ( istep, dihed, new_dihed ))

   atoms.set_dihedral(dihed_inds[0], dihed_inds[1], dihed_inds[2], \
                        dihed_inds[3], new_dihed, mask = flags)

   print("Dihedral indices are", dihed_inds[0], dihed_inds[1], dihed_inds[2], \
                        dihed_inds[3])
istep=0 


for dihed in dihedral_steps:
#
   new_dihed = dihedral_start + dihed
   print("Dihedral step: %d, %10.6f degrees, actual angle: %10.6f\n" % \
                       ( istep, dihed, new_dihed ))
#
   atoms.set_dihedral(dihed_inds[0], dihed_inds[1], dihed_inds[2], \
                        dihed_inds[3], new_dihed, mask = flags)
   fix_dihed = [atoms.get_dihedral(*dihed_inds), dihed_inds]
#
   cons=FixInternals(dihedrals=[fix_dihed])
#
   atoms.set_constraint(cons) 
#
   print("Starting optimisation.........", flush=True)
   opt = BFGS(atoms, trajectory='rot_check.traj')  #Put optimisation back in
   atoms.get_potential_energy()
   #opt.run(fmax=0.01)
   print("Energy obtained!")
   print("Optimisation completed........", flush=True)
#
  # if istep == 0:
  #     traj = Trajectory('rot_check.traj','w',atoms)
  # else:
  #     traj.write(atoms)

  # istep=istep+1

traj.close()


#

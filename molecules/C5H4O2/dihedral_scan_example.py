
import numpy as np
from scipy import sparse
import sys
import os
# Function to obtain HOME
from os.path import expanduser

from ase import Atom
from ase import Atoms
from ase.io import read
from ase import neighborlist
#
from ase.io import read, write, Trajectory
#
from ase.constraints import FixInternals
#
from carmm.run.aims_path import set_aims_command
from carmm.run.aims_calculator import get_aims_and_sockets_calculator
#
# Make modules available
#
home = expanduser("~")
set_module="%s/python/modules/data_io" % home
sys.path.append(set_module)
set_module="%s/python/modules/slab" % home
sys.path.append(set_module)
set_module="%s/python/modules/atom_settings" % home
sys.path.append(set_module)
set_module="%s/python/modules/num_digi" % home
sys.path.append(set_module)
set_module="%s/python/modules/vectors" % home
sys.path.append(set_module)
set_module="%s/python/modules/atom_vectors" % home
sys.path.append(set_module)
#
print("Current system path:")
print(sys.path)
#
# Import our own modules
#
from atom_vectors import *
#
print("ASE script starting....")
#
print(f"Arguements passed {len(sys.argv)}")
for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
#
PROG="AIMS"
# 
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
#
struct=sys.argv[1]
cores_per_task=sys.argv[2]
#
#
print(f"------------------------------------------------------------")
print(f"Work directory                 : %s" % workdir )
print(f"cores per task                 : %s" % cores_per_task )
print(f"------------------------------------------------------------")
#
atoms=read(struct)
natoms=len(atoms)
#
# Define control.in parameters
#
# New method that gives a default calculator
# Use dimensions 2 for slab calculations, this will turn on dipole correction

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
#
   atoms.set_dihedral(dihed_inds[0], dihed_inds[1], dihed_inds[2], \
                        dihed_inds[3], new_dihed, mask = flags)
#
   cons=FixInternals(dihedrals=[new_dihed,dihed_inds])
#
   atoms.set_constraint(cons) 
#
   print("Starting optimisation.........", flush=True)
   dyn = BFGS(atoms)  #Put optimisation back in
   dyn.run(fmax=0.01)
   print("Optimisation completed........", flush=True)
#
   if istep == 0:
       traj = Trajectory('rot_check.traj','w',atoms)
   else:
       traj.write(atoms)

   istep=istep+1

traj.close()
    




write("check_rot.in", atoms, format="aims")
#

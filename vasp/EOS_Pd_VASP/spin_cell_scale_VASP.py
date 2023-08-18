import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
#from ase.calculators.aims import Aims
from ase.calculators.vasp import Vasp2
from ase.optimize import BFGS
from ase.units import kJ
from ase.eos import EquationOfState
from ase.io import aims
#
# Code to do cell volume scaling for Murnighan plot.
# The script expects the structure to be in a POSCAR file in the 
# structure sub-directory.
#
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
# Define a matrix following Voigt convention
#
def voigt_mat(eee):
    vmat=[[eee[0], eee[5]/2.0,  eee[4]/2.0 ], [eee[5]/2.0, eee[1],  eee[3]/2.0 ], [eee[4]/2.0, eee[3]/2.0,  eee[2] ] ]

    return vmat
#
# Main code begins
#
debug=False
if debug:
  print("DEBUG run, just checking code....")
  machine="hawk"
  cores_per_task=40
  stdout_file="dummy.out"
  mydir="."
else:
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
#
if (debug):
  print("EMT run")
  csv_filename= "eos_EMT.csv" # filename for csv record 
#
else:
  print("VASP run")
  csv_filename= "eos_Pd_VASP.csv" 
#
# Read in the structure from a fhiaims file
# or a POSCAR file whichever is present.
# This will have been copied from your structure directory to the
# working directory, if POSCAR, it should be a VASP5 file, i.e. with the atom
# symbols above the line with number of ions.
#
geom_file="POSCAR"
cif_file="struct.cif"             

if os.path.isfile(geom_file):
   print("reading POSCAR file")
   atoms=read(geom_file)
#
elif os.path.isfile(cif_file):
   print("reading cif file")
   atoms=read(cif_file)
#
else:
   print("ERROR : No structure file, neither struct.car or castep.cell could be found")
   exit()
#
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
  calc = Vasp2(command=cmd_string, txt=stdout_file, xc='pbe', kpts=(7, 7, 7))
# Use new PAW for Au
  calc.set(setups={'Au': '_new'})
# Electronic structure settings
  calc.set(prec='Accurate', ediff=1E-8 )
  calc.set(encut=500, ispin=1, isym=0, ismear=0, sigma=0.2)
# Efficiency settings
  calc.set(algo='Fast', lreal='.FALSE.',addgrid='.TRUE.')
# vdW settings
  calc.set(ivdw= 11, vdw_s6=0.75)
#
# vdW settings
# here are several alternatives for the van der Waals interaction
#
# ivdw=12 : sets the D3 method of Grimme with BJ damping 
# calc.set(ivdw= 12)
#
# ivdw=20 : sets the Tkatchenko-Scheffler method
#  calc.set(ivdw=20)
#
# Attached the calculator to the atoms
  atoms.set_calculator(calc)
else:
  atoms.set_calculator(EMT())
#
# record cell vectors as supplied.
#
cell = atoms.get_cell()
print("cell as read:")
print(cell)
traj = Trajectory('bulk_eos.traj', 'w')
#
# Set up csv file to report the data for plotting
#
csv_file = open(csv_filename, 'w', newline='')
csv_writer = csv.writer(csv_file)
csv_headers1=['SpinState','Entry', 'Scale', 'a / Angs', 'b / Angs', 'c / Angs']
csv_headers2=['alpha / deg.', 'beta / deg.', 'gamma / deg.'] 
csv_headers3=['volume / Angs^3', 'Energy / eV','Strain / eV ']
csv_headers = csv_headers1+ csv_headers2+ csv_headers3  
csv_writer.writerow(csv_headers)
csv_headers=[' ', ' ', ' ', ' ',' ',' ',' ',' ',' ', ' ', 'e_1','e_2','e_3','e_4','e_5','e_6']
csv_writer.writerow(csv_headers)
#
index=0
#
# Need a unit matrix to use when operating on cell vectors
#
unit_mat = np.zeros((3,3))
unit_mat[0][0]=1.0
unit_mat[1][1]=1.0
unit_mat[2][2]=1.0
#
# define range of scalings to use
# These will be used to scale the sets of strain vectors that are defined.
#
min_scale= -0.04
max_scale= 0.04
num_scale= 3
#
#
num_points=-1
stress=[]
strain_sets=[]
num_sets=1
for iii in range(0,num_sets):
   strain_sets.append([0.0,0.0,0.0,0.0,0.0,0.0])
#
# define the patterns to use within a set:
#
# Set 0, Volume expansion:
strain_sets[0][0]=1.0
strain_sets[0][1]=1.0
strain_sets[0][2]=1.0
#
#
# Define the strain elements as a six component vectors
#
strain=np.empty((num_sets*num_scale,6))
eee = np.zeros(6)
#
## Attempting to put in spin states...
spin_states=[1,2] 
for spin_state in spin_states: 
   calc.set(ispin=spin_state)
   atoms.set_calculator(calc)
#
   for iset in range(0,num_sets):
      for x in np.linspace(min_scale, max_scale, num_scale):
#
# workout the actual strain to use and form the matrix
#
         for iii in range(0,6):
            eee[iii] = x*strain_sets[iset][iii]
#
         vmat = unit_mat + voigt_mat(eee)
         print("matrix to be used for cell change")
         print(vmat)
#
         new_cell =  np.dot(vmat, cell)
#
         atoms.set_cell(new_cell, scale_atoms=True)
         print("new cell: ", atoms.get_cell())
#
         lengths= atoms.cell.lengths()
         angles= atoms.cell.angles()
         print("new cell lengths: ", lengths) 
         print("new cell angles : ", angles)
#     
         if debug:
            print("No VASP run as debugging")
         else:
            opt = BFGS(atoms)
            opt.run(fmax=0.01)
            print(f"Optimisation reached required accuracy")
#
         num_points+=1
         strain[num_points]=eee
#
# Make name for optimised structure at this scaling level
#
         if x >= 0.0:
           xlab="%5.3f" % x
         else:
           xlab="m%5.3f" % np.absolute(x)
#
         xlab=xlab.replace(".","p")
         xlab=xlab.lstrip()
#
         new_fname="geometry_spin%d_" % spin_state  + xlab + "POSCAR"
#
         print("new_fname: ", new_fname)
#
# here we will write our new cell to a POSCAR file
 #     write(new_fname, atoms ,format='vasp')
#
         print("Strain: ", eee)
#      print("Stress: ", stress[num_points])


# Write data to the csv file
         csv_row1= [spin_state, index, x, lengths[0],lengths[1],lengths[2]]
         csv_row2= [angles[0],angles[1],angles[2],atoms.get_volume()]
         if debug:
            csv_row3= [99.9999,eee[0],eee[1],eee[2],eee[3],eee[4],eee[5]]
         else:
            csv_row3= [atoms.get_potential_energy(),eee[0],eee[1],eee[2],eee[3],eee[4],eee[5]]
         csv_row = csv_row1 + csv_row2 + csv_row3
         csv_writer.writerow(csv_row)

         index+=1
         num_points=-1
#
csv_file.close()


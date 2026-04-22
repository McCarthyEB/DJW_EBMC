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
# Function to obtain HOME
from os.path import expanduser
#
from ase import Atom
from ase import Atoms

from ase.optimize import BFGS
from ase.optimize import FIRE
from ase.constraints import FixAtoms
from ase.constraints import FixBondLength

from ase.io import read, write, Trajectory

from ase.calculators.vasp import Vasp2

#Import database modules
#
from ase.db import connect
from ase.build import molecule
from ase.calculators.emt import EMT

home = expanduser("~")
#
set_module="%s/python/DJW_EBMC/modules/slab" % home
sys.path.append(set_module)
set_module="%s/python/DJW_EBMC/modules/atom_settings" % home
sys.path.append(set_module)
set_module="%s/python/DJW_EBMC/modules/num_digi" % home
sys.path.append(set_module)
set_module="%s/python/DJW_EBMC/modules/vectors" % home
sys.path.append(set_module)
set_module="%s/python/DJW_EBMC/modules/classify_from_bitmask" % home
sys.path.append(set_module)

#
print("Current system path:")
print(sys.path)
#
# Import our own modules
#
from set_atoms import atom_formal_charges
#
from binary_tools import *
from vectors import *

#
# Local definition of atom...atom distance measurement
#
def atom_dist(atoms, i1, i2):
    #
    #  Inter-atomic vector
    #
    vec = np.array([atoms.positions[i2, 0] - atoms.positions[i1, 0], \
                    atoms.positions[i2, 1] - atoms.positions[i1, 1], \
                    atoms.positions[i2, 2] - atoms.positions[i1, 2]])
    #
    #  impose minimum image convention
    #
    latt = (atoms.get_cell()).copy()
    recip_latt = (atoms.get_reciprocal_cell()).copy()
    #
    print(latt)
    print(recip_latt)
    #
    vec = np.array([atoms.positions[i2, 0] - atoms.positions[i1, 0], \
                    atoms.positions[i2, 1] - atoms.positions[i1, 1], \
                    atoms.positions[i2, 2] - atoms.positions[i1, 2]])

    size = np.sqrt(np.vdot(vec, vec))
    vec = vec / size

    return size


#
# Main code begins
#
# debug: setting "True" allows checking of code without calling fhiaims
#        this can be run in the foreground to check set up before
#        commiting to a fhiaims calculation. Useful here for magnetic
#        ordering check.
#        setting "False" will run the fhiaims jobs so only do as part of a
#        job in a queue
#
debug = True
#
if (debug):
    print(f"DEBUG version, checking code, fhiaims will not be run")
    print(f"DEBUG foreground execution with additional printing")
    #
    workdir = os.getcwd()  # Directory where we will do the calculations
    geom_file = "supercell_0_12layer.in"  # structure file will be passed to script
    cores_per_task = "not defined in debug mode"  # cores assigned to tasks where available
    jobid = "not defined in debug mode"  # Unique identifier for this job
    stdout_file = "not defined in debug mode"  # stdout file will appear in your launch directory
    results_dir = os.getcwd()  # directory for results, this should be on your $HOME space
    traj_file = "not defined in debug mode"  # File for trajectory output
    traj_spin_file = "not defined in debug mode"  # File for trajectory output
    machine = "not defined in debug mode"  # string defining machine
else:
    print(f"Arguements passed {len(sys.argv)}")
    for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>7}: {arg}")
    #
    # this ase_vasp_opt.py script expects:
    # the first arguement to be the directory where we will run vasp,
    # second the number of cores available or a junk string if this is not required
    # The third to be the JOB_ID, a unique identifier to add to output information
    # The fourth is the path to the directory from which the job was launched
    # The fifth is the sub-directory for the structure for this particular run.
    # The sixth is the sub-directory for the sub_directory for this type of job, under the structure
    # The seventh to be the name of the machine we are running on.
    #           Currently this can be "hawk","thomas" or  "archer".
    #
    #
    workdir = sys.argv[1]  # Directory where we will do the calculations
    cores_per_task = sys.argv[2]  # cores assigned to tasks where available
    jobid = sys.argv[3]  # Unique identifier for this job
    stdout_file = "%s/vasp_%s_%s.stdout" % (
    sys.argv[4], sys.argv[5], jobid)  # stdout file for fhiaims will appear in your launch directory
    results_dir = "%s/%s/%s" % (
    sys.argv[4], sys.argv[5], sys.argv[6])  # directory for results, this should be on your $HOME space
    traj_file = "%s/%s_%s_spin1.traj" % (results_dir, sys.argv[5], jobid)  # File for trajectory output
    traj_spin_file = "%s/%s_%s_spin2.traj" % (results_dir, sys.argv[5], jobid)  # File for trajectory output
    machine = sys.argv[7].lower()  # string defining machine
#
print(f"------------------------------------------------------------")
print(f"Work directory               : %s" % workdir)
print(f"Geometry file                : %s" % geom_file)
print(f"stdout file for fhiaims      : %s" % stdout_file)
print(f"results_dir for this script  : %s" % results_dir)
print(f"trajectory file for this run : %s" % traj_file)
print(f"trajectory file for ISPIN 2  : %s" % traj_spin_file)
print(f"Machine running job          : %s" % machine)
print(f"------------------------------------------------------------")
#
# set command for running fhiaims, stdout file destination
# This is machine dependent, for a new machine add what would have been in
# the job control script to launch the calculation.
#
if 'hawk' in machine:
    cmd_string = "mpirun -n %s vasp_std > %s" % (cores_per_task, stdout_file)
elif 'thomas' in machine:
    cmd_string = "gerun vasp_gam > %s" % stdout_file
elif 'archer' in machine:
    cmd_string = "aprun -n $NPROC vasp_gam > %s" % stdout_file
else:
    if (debug):
        print(f"DEBUG: No command for fhiaims required or set.")
    else:
        print(f"ERROR: Unknown machine, do not know command for running fhiaims!")
#
# Read in the structure from a cif file
#
atoms = read("supercell_0_12layer.in")
write('check_start.cif', atoms, format='cif')
write('check_start.xyz', atoms, format='xyz')

#
# Set atom charges
# Formal charges just used to estimate the slab dipole
#
totq = atom_formal_charges(atoms)
#
zdip = 0
for iatom in range(0, len(atoms)):
    zdip += atoms[iatom].charge * atoms.positions[iatom][2]
#
zmax = atoms.positions[0][2]
zmin = zmax
for iatom in range(0, len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2], \
             atoms[iatom].charge))
    #
    # Identify top and bottom layers
    #
    if (atoms.positions[iatom][2] > zmax):
        zmax = atoms.positions[iatom][2]
    #
    if (atoms.positions[iatom][2] < zmin):
        zmin = atoms.positions[iatom][2]

print(" ")
print(f"After charge assignment, total system charge = %10.6f and Z-dipole moment = %10.6f" % (totq, zdip))
print(" ")
print(f"Identified zmax = %10.6f and zmin = %10.6f" % (zmax, zmin))
#
# Identify top and bottom layer
#
zdelta = 0.1
top_layer = []
bot_layer = []
for iatom in range(0, len(atoms)):
    if ((atoms.positions[iatom][2] > (zmax - zdelta)) \
            & (atoms.positions[iatom][2] < (zmax + zdelta))):
        top_layer.append(iatom)

    if ((atoms.positions[iatom][2] < (zmin + zdelta)) \
            & (atoms.positions[iatom][2] > (zmin - zdelta))):
        bot_layer.append(iatom)
#
num_top_layer = len(top_layer)
num_bot_layer = len(bot_layer)
#
print(f"Identified %d top layer atoms:" % len(top_layer))
for iii in range(0, num_top_layer):
    iatom = top_layer[iii]
    print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2], \
             atoms[iatom].charge))
#
print(f"Identified %d bottom layer atoms:" % len(bot_layer))
for iii in range(0, num_bot_layer):
    iatom = bot_layer[iii]
    print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2], \
             atoms[iatom].charge))
#
# work out possible patterns
#
for num_to_set in range(7, 8):
    flag_list_top = k_bits_on(num_to_set, num_top_layer)
    nlist = len(flag_list_top)
    all_dists = []
    #
    print
    print("Setting up ", num_to_set, "vacancies, need to test ", nlist, " arrangements.")
    print
    #
    # Work out number to expect in dist lists
    #
    num_in_dist_list = 0
    dist_list = []
    unique_list = []
    num_of_unique = []
    for ilist in range(0, num_to_set):
        num_in_dist_list += ilist
        if (ilist < num_to_set - 1):
            dist_list.append([j for j in range(ilist + 1, num_to_set)])

    print("num_in_dist_list = ", num_in_dist_list)
    for ilist in range(0, num_to_set - 1):
        print(dist_list[ilist])

    dists = np.empty([num_in_dist_list])
    num_dists = 0
    #
    #
    for i in range(0, nlist):
        struct_id = flags_to_int(flag_list_top[i])
        #      print(flag_list_top[i]," : ", struct_id)
        #
        #
        this_atoms = atoms.copy()
        top_vacs = atoms.copy()
        #
        # Remove num_to_set atoms top and bottom of slab
        #
        flag = []
        for iatom in range(0, num_top_layer):
            if (flag_list_top[i][iatom] == 1):
                flag.append(top_layer[iatom])
                #            print(f"Flagging atom %d" % top_layer[iatom])
                flag.append(bot_layer[iatom])
        #
        # delete four top and four bottom atoms
        #
        del this_atoms[flag]
        #
        # Make top_vacs as atom set for just the atoms about to be deleted
        #
        flag = []
        for iatom in range(0, len(atoms)):
            keep = -1
            for jatom in range(0, num_top_layer):
                if (iatom == top_layer[jatom] and flag_list_top[i][jatom] == 1):
                    keep = iatom
            #
            if (keep < 0):
                flag.append(iatom)
        #        else:
        #           print(f"Keeping atom %d in top_vacs" % iatom)
        #
        # delete all except four top atoms
        #
        del top_vacs[flag]
        #
        #      for iatom in range(0,len(top_vacs)):
        #         print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f"                 \
        #             % (iatom, top_vacs.symbols[iatom], top_vacs.positions[iatom][0], \
        #                top_vacs.positions[iatom][1], top_vacs.positions[iatom][2],   \
        #                top_vacs[iatom].charge))
        #
        dists = []
        for ilist in range(0, num_to_set - 1):
            #         print("Preparing list contribution ",ilist," using dist_list ",dist_list[ilist])
            dists = np.append(dists, top_vacs.get_distances(ilist, dist_list[ilist], mic=True))
        #     print("Distances:")
        #     print(dists)
        #
        # Compare this set of distances with those already seen
        #
        if (i == 0):
            all_dists = dists
            num_dists = 1
            unique_list = [struct_id]
            num_of_unique = [1]
            #
            vac_file = "top_vac_atoms_%d_%d.cif" % (num_to_set, struct_id)
            write(vac_file, top_vacs, format='cif')
        #
        else:
            ind = 0
            matched = False
            for which_dist in range(0, num_dists):
                test_dist = []

                for jdist in range(0, num_in_dist_list):
                    test_dist.append(all_dists[ind])
                    ind += 1

                if (match_vectors(dists, test_dist)):
                    matched = True
                    imtch = which_dist
                    break
            #
            if (matched):
                #            print("Matched with existing list...")
                num_of_unique[imtch] += 1
                #
                # Check a small degeneracy set
                if (num_to_set == 4 and unique_list[imtch] == 15):
                    vac_file = "top_vac_atoms_%d_%d_%d.cif" % (num_to_set, struct_id, num_of_unique[imtch])
                    write(vac_file, top_vacs, format='cif')


            else:
                num_dists += 1
                all_dists = np.append(all_dists, [dists])
                unique_list.append(struct_id)
                num_of_unique.append(1)
                #
                # Workout dipole for this set
                #
                print(flag_list_top[i], " : ", struct_id)
                zdip = 0.0
                for iatom in range(0, len(this_atoms)):
                    zdip += this_atoms[iatom].charge * this_atoms.positions[iatom][2]
                print(f"Z-dipole for this : %10.6f" % zdip)
                #
                vac_file = "top_vac_atoms_%d_%d.cif" % (num_to_set, struct_id)
                write(vac_file, top_vacs, format='cif')
                slab_file = "slab_atoms_%d_%d.cif" % (num_to_set, struct_id)
                write(slab_file, this_atoms, format='cif')

    #
    print(" ")
    print(f"Thats %d patterns of which %d are unique" % (nlist, num_dists))
    #
    # Report the data
    #
    csv_headers = ['index', 'structure index', 'Number occurances']
    csv_filename = "degen_5b5_%d.csv" % num_to_set
    csv_file = open(csv_filename, 'w', newline='')
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(csv_headers)

    for ind in range(0, num_dists):
        #
        csv_row = [ind, unique_list[ind], num_of_unique[ind]]
        #
        csv_writer.writerow(csv_row)
    #
    csv_file.close()

    #
    if (debug):
        print(f"DEBUG: Will not run fhiaims..")
        print(f"DEBUG: just checking slab charge and dipole in Z-direction")
    #
    else:
        #
        # Set up for fhiaims calculations
        #
        print(f"Will run fhiaims")


"""
Data input/output routines written for specific tasks
Last modified 27th March 2024

@author: dave
"""
import csv
import sys
import os

# Will use numpy for arrays
import numpy as np
from scipy import special

from ase.io import read, write, Trajectory
from ase import Atom
from ase import Atoms
#
from binary_tools import *
#
# read_config_data, produced for PdH configurations calculations to read the 
#                   energy and configuration information from csv files.
#                   Uses a header line to set which data item is read from which
#                   column of the data file, pointed to by fileobj.
#                   This version is written to use the raw slab data to calculate 
#                   the adsorption energies from the reference for the clean slab
#                   and isolated molecule if supplied.
#
# Update 9th April 2024 to remove site identity so that configs can refer to any system with
# two distinct site types
#
# Note: for Pd111 case sites_1 is hcp sites_2 is fcc
#       label_1 and label_2 should refer to the header strings used in the csv to identify indicies for 
#       site_1 and site_2, e.g. for Pd111 label_1 = "hcp" label_2 = "fcc"
#
def read_config_data(fileobj, nsites_1, nsites_2, label_1, label_2, have_ref, Eref_clean, Eref_molecule):
#
    lines = fileobj.readlines()
    nlines = len(lines)  
#    print('File contains:',nlines,' lines')
#
#     Initiate arrays for collecting information from data file
#     sf is structure filename
#     hcp_i is the decimal index for the hcp sites
#     fcc_i is the decimal index for the fcc sites
#     sub_i is the decimal index for the sub sites where present
#     g_i is the degeneracy for this configuration
#     da is the average displacement for this set of configurations
#     Ea is the calculated adsorption atoms ( for all atoms in original per H2 for csv files )
#
    sf = []
    site1_i = []
    site2_i = []
    sub_i = []
    g_i = []
    da = []
    Ea = []
    Eslab = []
    ZPE = []
    freq_list=[]
    have_index=False
    nconfig=-1
##
## This loop will now run through all the lines of the file so need to
## extract information from the lines as appropriate, hence the if blocks
##
    for iline in range(len(lines)):
        line=lines[iline]
        linesplit = line.split(',')
        nwords = len(linesplit)
        keep = True
#
# Read in the data line : summary file version, Jan 2023 
#
        if have_index:
#            print("line words:")
#            for iii in range(0,nwords):
#                print(linesplit[iii])
            if ("unmatched" in linesplit[1]):
               keep = False
#
            if keep:
               nconfig+=1
               sf.append(linesplit[1])
#
#               print("DEBUG: nconfig: %d filename: %s" % (nconfig, sf[-1]))
# 
               if site1_col > 0:
                  site1_i.append(int(linesplit[site1_col]))
               if site2_col > 0:
                  site2_i.append(int(linesplit[site2_col]))
#
# if site1_col and site2_col not defined try to get indices from file name
#
               if site1_col < 0 and site2_col < 0:
#
                  temp=sf[-1].split(".")
                  if len(temp) > 0:
                    stem=temp[0]
#
                    scored=stem.split("_")
                    num_ads=int(scored[2])
                    site1_i.append(int(scored[3]))
                    site2_i.append(int(scored[4]))
#
               if sub_col > 0:
                  sub_i.append(int(linesplit[sub_col]))
#
               if degen_col > 0:
                   g_i.append(int(linesplit[degen_col]))
#
               if not have_ref:
                  Ea.append(float(linesplit[ea_col]))
               else:
                  if ea_col > 0:
                      E_check=float(linesplit[ea_col] )
#
               if da_col > 0:
                  if ("unmatched" in linesplit[da_col]):
                     da.append(-1.0)
                  else:
                     da.append(float(linesplit[da_col]))
#
               if eslab_col > 0:
                  Eslab.append(float(linesplit[eslab_col]))
#
               if zpe_col > 0:
                  ZPE.append(float(linesplit[zpe_col]))
#           
# Check consistency of structure for coverage
#
               site1_array=int_to_flags(site1_i[nconfig],nsites_1)
               site2_array=int_to_flags(site2_i[nconfig],nsites_2)
#
# If sub-surface atoms present need to add 1 to nads so that the sub-surface atom is included
#
               if sub_col > 0:
                  nads=np.sum(site1_array)+np.sum(site2_array) + 1
               else:
                  nads=np.sum(site1_array)+np.sum(site2_array)
#
               if (nconfig == 0):
                   nads_chk=nads
               else:
                   if (nads != nads_chk):
                       print("ERROR: inconsistent adsorbate number found in configs list")
                       sys.exit(0)
#
               if have_ref and eslab_col > 1:
                  fnadso2=float(nads)/2.0
                  Eads_per_H2 = (Eslab[-1] - Eref_clean - fnadso2*Eref_molecule)/fnadso2
                  Ea.append(Eads_per_H2)
#
                  if ea_col > 0:
                     diff = Eads_per_H2 - E_check
#                     print("%s: Derived Eads_per_H2 as %10.6f, raw data: %10.6f, diff : %10.6f " % \
#                                                    ( sf[-1], Eads_per_H2, E_check, diff ))
                     if (diff > 0.05):
                        print("ERROR: difference between raw data and calculated larger than expected, CHECK")
                        exit(0)
#
# Read in nads*3 sets of vibrational wavenumbers
# check that the frequencies are there first
# 
               have_freq = False
               freq_start=zpe_col+1
               if len(linesplit) > freq_start + 3*nads - 1 :
                  have_freq = True
#
               flist=[]
               if have_freq:
                 for iii in range(0,3*nads):
                    flist.append(float(linesplit[freq_start+iii]))

               freq_list.append(flist) 
#           
        elif nwords>0 and 'index' in linesplit[0]:
            have_index=True
            site1_col=-1
            site2_col=-1
            sub_col=-1
            degen_col=-1
            ea_col=-1
            eslab_col=-1
            zpe_col=-1
            da_col=-1
            for iword in range(0,nwords):
               if ( label_1 in linesplit[iword]):
                   site1_col = iword
               if ( label_2 in linesplit[iword]):
                   site2_col = iword
               if ( "sub" in linesplit[iword]):
                   sub_col = iword
               if ( "degen" in linesplit[iword] or "occurances" in linesplit[iword]):
                   degen_col = iword
               if ("Eads" in  linesplit[iword]):
                   ea_col = iword
               if ("shift" in  linesplit[iword]):
                   da_col = iword
               if ("e_slab" in linesplit[iword].lower()):
                   eslab_col = iword
               if ("zpe" in  linesplit[iword].lower()):
                   zpe_col = iword
            
            print("Found headers line...")
            print(line)
#
        else:
          print("ERROR: Unrecognised data line with nwords= %d" % nwords)
          print(line)
          sys.exit(0)

    return (sf, site1_i, site2_i, sub_i, g_i, da, Ea, Eslab, ZPE, nconfig+1, nads, freq_list)
#

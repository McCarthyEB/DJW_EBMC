#mode_follow-m5.castep
#Final energy, E             =  -118869.5815315     eV
#RMS Force, (   all atoms ) :   0.125254
#RMS Force, ( constrained ) :   0.100022
#Atom H 25 to Atom C 26 =   1.474437 Angs.
# -*- coding: utf-8 -*-
import numpy as np
import sys
import csv
import os
#
# Main code begins
#
print(f"Arguements passed {len(sys.argv)}")
#
if len(sys.argv) < 2:
   print("ERROR : need file name containing mode follow calculation summary")
   exit(0)
else:
   file = sys.argv[1]          # file name passed should contain data list   
#
print("Reading data from file %s" % file)
#
fp=open(file,'r')

lines = fp.readlines()
nlines = len(lines)
#
step=[]
energy=[]
rms_all=[]
rms_con=[]
dist=[]

for iline in range(0,nlines):
#
   line=lines[iline]
   linesplit=line.split()
   nwords=len(linesplit)
#
   if (nwords==1 or nwords == 3) and "mode" in linesplit[0]:
     print("filename: %s" % linesplit[0])
#
     sss=linesplit[0].split('-')
     ppp=sss[1].split('.')
     print("%s" % ppp[0])
#
     if "m" in ppp[0]:
       s_num=ppp[0].replace("m","-")
     else:
       s_num=ppp[0]
#
     step.append(float(s_num))
#
   if nwords==6 and "energy" in linesplit[1]:
     energy.append(float(linesplit[4]))
#
   if nwords==8 and "Force" in linesplit[1] and "all" in linesplit[3]:
        rms_all.append(float(linesplit[7]))
#
   if nwords==7 and "Force" in linesplit[1] and "con" in linesplit[3]:
        rms_con.append(float(linesplit[6]))
#
   if nwords==10 and "Atom" in linesplit[0] and "Angs" in linesplit[9]:
        dist.append(float(linesplit[8]))
#
emin=np.min(energy)
#
# Write csv summary file
#
csv_pntr=open("mode_follower_summary.csv", 'wt')
#
have_dist = False
if  (len(step) == len(dist)):
   have_dist = True
#
if have_dist: 
   title_line="Step, dist / Angs, Energy / eV, Rel. Energy / eV, Force (all), Force (constrained)\n"
else:
   title_line="Step, Energy / eV, Rel. Energy / eV, Force (all), Force (constrained)\n"
#
csv_pntr.write(title_line)
#
for iii in range(0, len(step)):
   erel=energy[iii]-emin
#
   if have_dist: 
      data_line="%6.0f, %10.6f, %10.6f, %10.6f, %10.6f, %10.6f\n"\
             % (step[iii], dist[iii], energy[iii], erel, rms_all[iii], rms_con[iii])
   else:
      data_line="%6.0f, %10.6f, %10.6f, %10.6f, %10.6f\n"\
             % (step[iii], energy[iii], erel, rms_all[iii], rms_con[iii])
#
   csv_pntr.write(data_line)
#
csv_pntr.close()

fp.close()




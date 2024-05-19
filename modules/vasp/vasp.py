#
# Routines for setting up vasp calculations
#
import numpy as np

# Check if we have a number for cores_per_task and use to estimate ncore
#
def best_ncore(cores_per_task):
#
    core_diff=[]
    cores_list=[]
    test=cores_per_task.isdigit()
    if test:
        tot_cores=float(cores_per_task)
        ncore=int(np.sqrt(tot_cores))
        div = tot_cores/ncore
        print("cores_per_task= ", cores_per_task, " ideal ncore = %d" % ncore)
        print("cores_per_task/ncore = ", div)
#
        if tot_cores%ncore == 0:
           print("Already integer so can use ncore = %d" % ncore)
           return ncore
        else:
           print("This needs to be integer....")
#
        for n in range(2,int(np.sqrt(tot_cores))):
            print("n: ", n, " cores_per_task/n = ", tot_cores/n)
            if tot_cores%n == 0:
              m = int(tot_cores/n)
              cores_list.append(n)
              cores_list.append(m)
              core_diff.append(abs(ncore-n))
              core_diff.append(abs(ncore-m))
#
        print(cores_list)
        print(core_diff)
#
        ibest=-1
        diff_best=-1
        for indx in range(0,len(core_diff)):
           if indx==0:
              ibest=indx
              diff_best=abs(core_diff[indx])
           else:
              if abs(core_diff[indx]) < diff_best:
                 diff_best=abs(core_diff[indx])
                 ibest=indx
#
        if ibest < 0:
          print("Cannot divide cores_per_task into factors, putting ncore = 1")
          ncore=1
        else:
          ncore=cores_list[ibest]
          print("Best to use: %d i.e. ncore= %d" % (ibest, ncore))
#
    return ncore


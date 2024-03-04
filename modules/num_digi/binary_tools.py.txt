import numpy as np
import sys
import csv
import os
#
# Function to obtain HOME
from os.path import expanduser
#
import numpy as np

def flags_to_int(array):
#
# Set ndigi_max to the maximum bit precision of integers on
# your machine, usually 32 or 64. This controls the available
# bits for defining integers.
#
   ndigi_max = 32
   ndigi = len(array)
   if (ndigi > ndigi_max):
     print(f"ERROR: array greater than %d sent to flags_to_int exceeds precision level" % ndigi_max)
     return -1
#
   num=0
   pow_of_two=1
   for i in range(0,len(array)):
#
       num=num+array[i]*pow_of_two
       pow_of_two=pow_of_two*2

   return num
#
def int_to_flags(num, ndigi):
#
# Take an integer over to the equivalent binary array
#
  pow_of_two=1
  array=[]
  for i in range(0,ndigi):
#
    array.append(int(np.bitwise_and(num, pow_of_two)/pow_of_two))
    pow_of_two=pow_of_two*2
#
  return array
#
# Routine for generating all n digit binary numbers with only k digits set
# Taken from stackoverflow webpage
#
from itertools import combinations
from functools import reduce # not necessary in python 2.x

def k_bits_on(k,n):
       one_at = lambda v,i:v[:i]+[1]+v[i+1:]
       return [tuple(reduce(one_at,c,[0]*n)) for c in combinations(range(n),k)]




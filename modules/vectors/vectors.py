#
# Routines for vector manipulation or comparison
# Also incorporates some linear algebra
# Last updated 21st March 2024, Dave
#
# Functions defined:
#
#  match_vectors(v1,v2,tol) : Compare vectorsi v1 and v2 looking for matchesi within tol
#                             components are compared but can be in any order so
#                             [3,2,1] matches [1,3,2]
#  match_coords(v1,v2,tol)  : Match vectors requiring components to be in same order
#  vector_dist(v1,v2)       : Return the scalar distance between v1 and v2.
#  best_plane(coords)       : Return the normal vector corresponding to the best plane for
#                             points in the coords list.
#  crossing(l1x1, l1y1,     : Return the crossing pointi (x,y) of lines l1 and l2 
#           l1x2, l1y2,     : defined from points.
#           l2x1, l2y1,     : This is for interpolation on a 2D grid.
#           l2x2, l2y2 )    :
#  voigt_mat(eee)           : define a strain matrix using Voigt convention.
#                             eee is a six component vector. The function returns the 
#                             corresponding 3x3 matrix
#
import numpy as np
from numpy import linalg as LA
#
def unit_vector(vector):
# 
    return vector / np.linalg.norm(vector)
#
# angle function that copes with aligned and anti-aligned cases correctly

def vec_angle(v1, v2):
#
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#
# Compare vectors looking for matches in any order
#
def match_vectors(v1,v2,tol):
#
  length=len(v1)
#
# Check vectors are same length
#
  if (length != len(v2)):
#    print("Warning: Trying to match vectors of unequal length")
    return False, 0.0
#
# Match vectors if both zero length, not if we are here they are already
# the same length.
#
  if length == 0:
    return True, 0.0
#
# tol is the tolerence with which to compare elements
#
#  tol = 0.01
  tol2= tol*tol
#
  diff = (np.sum(v1) - np.sum(v2))/length
  if ( diff * diff > tol2 ):
#    print(f"Array sums differ: diff = %10.6f" % diff)
    return False, 0.0
#  else:
#    print("sums match")
#
  match1=np.zeros(length)
  match2=np.zeros(length)
#
  tdiff2=0.0
  for i in range(0,length):
#
     for j in range(0,length):
#
# Check that v2[j] has not already been matched
#
        if (match2[j] == 0):
           diff = v1[i] - v2[j]
           diff2 = diff*diff
           if (diff2 < tol2 ):
#                print("Matched: ", i, " with ", j)
              match1[i] = 1
              match2[j] = 1
              tdiff2=tdiff2+diff2
              break

#  print("Sum of matches: ", np.sum(match1), np.sum(match2))
  if ( np.sum(match1) == length and np.sum(match2) == length ):
      return True, tdiff2/length
  else:
      return False, 0.0
#
# match vectors requiring components to be in same order
#
def match_coords(v1,v2,tol):
#
  length=len(v1)
#
# Check vectors are same length
#
  if (length != len(v2)):
    print("Warning: Trying to match vectors of unequal length")
    return False, 0.0
#  else:
#    print("lengths match")
#
# tol is the tolerence with which to compare elements
#
  tol2= tol*tol
#
  tdiff2=0.0
  for i in range(0,length):
#
     diff = v1[i] - v2[i]
     diff2 = diff*diff
     tdiff2=tdiff2+diff2

  if (tdiff2 < tol2 ):
     return True, tdiff2/length
  else:
      return False, tdiff2/length
#
# simple vector distance
#
def vector_dist(v1,v2):
#
  length=len(v1)
#
# Check vectors are same length
#
  if (length != len(v2)):
    print("Warning: Trying to match vectors of unequal length")
    return 0.0
#
  tdiff2=0.0
  for i in range(0,length):
#
     diff = v1[i] - v2[i]
     diff2 = diff*diff
     tdiff2=tdiff2+diff2

  return np.sqrt(tdiff2)
#
# Find best plane for a set of 3D cartessian points.
# Expect points as a list of triples.
# [[1,2,3],[2.3.4].....]
#
def best_plane(coords):

  if len(coords) < 3:
    print("ERROR: cannot find best plane with less than 3 points....")
    exit(0)
  else:
    print("Best plane routine working on:")
    for i in range(0,len(coords)):
       print("%d %10.6f  %10.6f  %10.6f\n" % \
             ( i, coords[i][0], coords[i][1], coords[i][2]))

  tmp_A = []
  tmp_b = []
  for i in range(len(coords)):
    tmp_A.append([coords[i][0], coords[i][1],1])
    tmp_b.append(coords[i][2])
#
  b = np.matrix(tmp_b).T
  A = np.matrix(tmp_A)
  fit = (A.T * A).I * A.T * b
  errors = b - A * fit
  residual = np.linalg.norm(errors)
#
  nt=[float(fit[0]), float(fit[1]), -1.0]
  normal=nt / np.linalg.norm(nt)
#
  print("solution:")
  print("%f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
  print("Normal vector: " , normal)
  print("deviations from plane:")
  print(errors)
  print("residual:")
  print(residual)
#
  return normal
#
# find crossing point of lines l1 and l2 defined from points
#
def crossing(l1x1, l1y1, l1x2, l1y2, l2x1, l2y1, l2x2, l2y2 ):
#
# define lines in y=mx + c form from points
    m1= (l1y2 - l1y1)/(l1x2 - l1x1)
    c1= -m1*l1x1 + l1y1
#
    m2= (l2y2 - l2y1)/(l2x2 - l2x1)
    c2= -m2*l2x1 + l2y1
#
    x=(c2-c1)/(m1-m2)
    y=m1*x+c1
#
    return(x,y)
#
# Define a matrix following Voigt convention
#
def voigt_mat(eee):
    vmat=[[eee[0], eee[5]/2.0,  eee[4]/2.0 ], [eee[5]/2.0, eee[1],  eee[3]/2.0 ], [eee[4]/2.0, eee[3]/2.0,  eee[2] ] ]

    return vmat

#
#
#
 










#include <stdio.h>

/* routine to give the vector joining to cartessian co-ordinates */

void join_vectors(double *p_A, double *p_B, double *p_A_to_B)
{

  int iloop;

  for (iloop=0; iloop < 3; iloop++)
      {
         *(p_A_to_B+iloop)= *(p_B+iloop)-*(p_A+iloop);
      }

}


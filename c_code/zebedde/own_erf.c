/******************************************************/
/**** Routine to give values of the error function ****/
/**** Dave Willock May 7th 1996                    ****/
/******************************************************/

#include <stdio.h>

double erf_mac(double x, double accuracy);
 
double erf_asymptotic(double x, double accuracy);
 
double own_erf(double x, double accuracy)
{
#include "own_maths.h"
double erf_ans;

if (x < 4.0)
  {
erf_ans = erf_mac(x, accuracy);
  }
else 
  {
erf_ans = erf_asymptotic(x, accuracy);
  }

return erf_ans;
}

/******************************************************/
/**** Routine to give values of the error function ****/
/**** Dave Willock May 7th 1996                    ****/
/******************************************************/

#include <stdio.h>
#include <math.h>

double erf_asymptotic(double x, double accuracy)
{
#include "own_maths.h"
int sign, power, term_num;

double erf, term, term_abs, x_prod, term_abs_last, erf_last;
double denominator, numerator;

/******************************************************/
/**** Asymptotic form  ********************************/
/**** From Advanced Engineering Mathematics ***********/
/**** Erwin Krejszig Pub. J Wiley & Sons    ***********/
/**** take the sum till terms start to increase *******/ 
/**** see page 723 of the reference!            *******/ 
/******************************************************/

erf= 1/x;
x_prod=x;
term_num= 0;
term_abs= erf;
sign= 1;
numerator= 1;
denominator=1;
power= -1;

do
  {
     term_num++;

     power   += 2;
     sign     = -1*sign;
     numerator= numerator*power;
     denominator= denominator*2;
     x_prod   = x_prod*x*x;

     term_abs_last= term_abs;
     term_abs = numerator/(x_prod*denominator);

if (term_num < 10 && x > 1.99 && x < 2.01)
{
printf("For asymptotic form: numerator= %10.6f, denominator= %10.6f, x_prod= %10.6f\n"
                                  ,numerator, denominator, x_prod);
printf("term_abs = %10.6f\n\n",term_abs);
}
     term     = sign*term_abs;

     erf_last = erf;
     erf += term;

  } while (term_abs < term_abs_last);

erf = 1-erf_normalise*erf_last*0.5*exp(-x*x);

return erf;
}

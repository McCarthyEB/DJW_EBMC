/******************************************************/
/**** Routine to give values of the error function ****/
/**** Dave Willock May 7th 1996                    ****/
/******************************************************/

#include <stdio.h>

double erf_mac(double x, double accuracy)
{
#include "own_maths.h"
int sign, power, term_num;

double erf, term, term_abs, x_prod;

/******************************************************/
/**** NB: factorial has to be a double as it gets *****/
/****     large                                   *****/
/******************************************************/

double factorial;

/******************************************************/
/**** Maclaurian series  ******************************/
/**** From Advanced Engineering Mathematics ***********/
/**** Erwin Krejszig Pub. J Wiley & Sons    ***********/
/******************************************************/

erf= x;
x_prod=x;
term_num= 0;
sign= 1;
factorial= 1;
power= 1;

do
  {
     term_num++;

     power   += 2;
     sign     = -1*sign;
     factorial= factorial*term_num;
     x_prod   = x_prod*x*x;

     term_abs = x_prod / (factorial * power);
     term     = sign*term_abs;

     erf += term;

  } while (term_abs > accuracy);

erf = erf_normalise*erf;

return erf;
}

#include <stdio.h>
/***************************************************************/
/****** To perform c = <a|b> where : ***************************/
/****** a and b are vectors of length dimension ****************/
/***************************************************************/

double bra_x_ket( double *p_a, double *p_b, int dimension )
{

int i;

double c;

c= 0;
for (i=0; i < dimension; i++)
  {
    c += *p_a * *p_b;
    p_a++;
    p_b++;
  }
return c;
}

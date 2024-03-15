#include <stdio.h>
/***************************************************************/
/****** To perform A= A + factor*( |a><b| + |b><a| ) where : ***/
/****** |a> and |b> are kets A is a symmetric matrix ***********/
/****** sent as its top diagonal *******************************/
/***************************************************************/
/****** Dave Willock September 1995 ****************************/
/***************************************************************/

void ket_bra_to_matrix( double *p_matrix, double *p_ai, double *p_bi,
                        double factor, int dimension )
{

int i,j;

double *p_bj, *p_aj, addition;

for (i=0; i < dimension; i++)
  {
    *p_matrix += factor * 2.0 * ( *p_bi  * *p_ai);

    p_matrix++;

    p_bj= p_bi;
    p_aj= p_ai;
    
    for (j=i+1; j < dimension; j++)
      { 
        p_bj++;
        p_aj++;

        addition= factor * ( *p_bi * *p_aj + *p_bj * *p_ai);
        *p_matrix +=  addition;      

        p_matrix++;

      }
    p_ai++;
    p_bi++;
  }
return;
}

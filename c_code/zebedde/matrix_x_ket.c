#include <stdio.h>
/***************************************************************/
/****** To perform |b> = A |a> where : *************************/
/****** |a> and |b> are kets A is a symmetric matrix ***********/
/****** sent as its top diagonal *******************************/
/***************************************************************/
/****** Dave Willock September 1995 ****************************/
/***************************************************************/

void matrix_x_ket( double *p_matrix, double *p_ai, double *p_bi,
                   int dimension )
{

int i,j;

double *p_bj, *p_bzero, *p_aj;

p_bzero= p_bi;
for (i=0; i < dimension; i++)
  {
     *p_bzero = 0;
     p_bzero++;
  }

for (i=0; i < dimension; i++)
  {
    *p_bi += (*p_matrix) * (*p_ai);

    p_matrix++;

    p_bj= p_bi;
    p_aj= p_ai;
    
    for (j=i+1; j < dimension; j++)
      { 
        p_bj++;
        p_aj++;

        *p_bi += (*p_matrix) * (*p_aj);      
        *p_bj += (*p_matrix) * (*p_ai);

        p_matrix++;

      }
    p_ai++;
    p_bi++;
  }
return;
}

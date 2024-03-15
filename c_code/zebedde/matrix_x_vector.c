#include <stdio.h>
/***************************************************************/
/****** To perform b = A a where : *****************************/
/****** a and b are vectors A is a symmetric matrix ************/
/******R sent as its top diagonal ******************************/
/***************************************************************/

void matrix_x_vector( double *p_matrix, double *p_ai, double *p_bi,
                      int dimension )
{

int i,j;

double *p_bj, *p_aj;

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

#include <stdio.h>
#include <math.h>

/* routine to multiply a three component vector by a 3x3 matrix */

void matrix_vec_mult( double *p_matrix, double *p_vector, double *p_answer)
{

  int icoloum;

  double *p_index_vec;

   for (icoloum=0; icoloum< 3; icoloum++)
     {
        *p_answer= 0.0;
        p_index_vec= p_vector;

        *p_answer += (*p_matrix)*(*p_index_vec);
        p_index_vec++;
        p_matrix++;

        *p_answer += (*p_matrix)*(*p_index_vec);
        p_index_vec++;
        p_matrix++;

        *p_answer += (*p_matrix)*(*p_index_vec);
        p_matrix++;

        p_answer++;
     } 

 return;
}


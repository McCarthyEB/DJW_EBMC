#include <stdio.h>
#include "maxima.h"
#include "structures.h"

/* routine to give the vector joining to atoms */

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B)
{
         *p_A_to_B     = p_B->x - p_A->x;
         *(p_A_to_B+1) = p_B->y - p_A->y;
         *(p_A_to_B+2) = p_B->z - p_A->z;
}


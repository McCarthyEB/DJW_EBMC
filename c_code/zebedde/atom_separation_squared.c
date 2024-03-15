#include <stdio.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

/* routine to give the vector square of the separation between two atoms */
/* with minimum image if the pbc flag is set                             */

void min_image( double *x, double *y, double *z);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc)
{

double dx,dy,dz,r2;

    dx  = p_B->x - p_A->x;
    dy  = p_B->y - p_A->y;
    dz  = p_B->z - p_A->z;
   
/***************************************************************************/
/******* if this is a periodic pore use the nearest image convention *******/
/***************************************************************************/

    if (pbc)
      {
         min_image( &dx, &dy, &dz);
      }

    r2= dx*dx + dy*dy + dz*dz;

    return r2;
}


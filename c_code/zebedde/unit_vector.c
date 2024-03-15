#include <stdio.h>
#include <math.h>

/***************************************************************************/
/* routine to give the unit vector in the direction of the supplied vector */
/* modified to also supply the size of the original vector DJW 27 Nov 96   */
/***************************************************************************/

/* Functions required for this program -------------------------*/

double size_vector(double *p_vector);

/*--------------------------------------------------------------*/

void unit_vector(double *p_vector, double *p_size)
{

  int iloop;

  *p_size = size_vector( p_vector );

   for (iloop=0; iloop < 3; iloop++)
    {
         *p_vector= (*p_vector)/(*p_size);
         p_vector++;
    }
}


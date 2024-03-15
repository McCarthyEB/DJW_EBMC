#include <stdio.h>
#include "maxima.h"
#include "structures.h"

/***************************************************/
/* subroutine for copying one type list to another */
/* the from list over rides the to                 */
/*                                                 */
/* expect num_from to be real number not upper     */
/* index. Dec. 2013                                */
/*                                                 */
/* Dave Willock August 1995                        */
/***************************************************/

void copy_types( atom_number *p_to, int *num_to,
                 atom_number *p_from, int num_from )

{

int iloop;

   for (iloop=0; iloop < num_from; iloop++)
      {
         *p_to= *p_from;
         p_to++;
         p_from++;
      }
   *num_to= num_from;

  return;
}





/************************************************************/
/* shuffle_neighbours                                       */
/* shuffle neighbours to cover up for lost hydrogens        */
/*                                                          */
/* Started DJW 14/7/95                                      */
/************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void shuffle_neighbours(atom *p_template, int place_for_one,
                        int place_for_two, int num_template_atoms)
{
#include "header.h"
int iatom, ineigh, this_neigh;

atom *p_atom;

 p_atom= p_template-1;
 for (iatom=0; iatom <= num_template_atoms; iatom++)
   {
     p_atom++;
     for (ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
       {
          this_neigh= p_atom->neighb[ineigh];
          if (this_neigh > place_for_two)
            {
               (p_atom->neighb[ineigh]) -= 2;
            }
          else if (this_neigh > place_for_one)
            {
               (p_atom->neighb[ineigh]) --;
            }
       }
   }

return; 
}


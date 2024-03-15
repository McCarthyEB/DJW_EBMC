
/*************************************************************************/
/* assign_stretch.c                                                      */
/* Sets up the quartic stretch parameter indexing for each atom with all */
/* its neighbours, i.e. neighbour info must be assigned before yee enter */
/* here.                                                                 */
/*                                                                       */
/* Started Dave Willock 16/1/97                                          */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void assign_stretch(atom *p_molecule, int num_atoms)
{
#include "header.h"
int iatom,jpot,ineigh, neigh_index, iequiv;
char atom_pot[4];
char neigh_pot[4];
	
BOOLEAN found_pot_flag;

atom *p_atom, *p_neigh;

/******************************************************************/
/******* assign potential index to each atom in list **************/	
/******************************************************************/

  for (iatom=0; iatom<=num_atoms; iatom++)
    {
      p_atom = p_molecule+iatom;	
      strcpy( atom_pot, p_atom->pot);

/******************************************************************/
/******* Check equivalence table **********************************/
/******************************************************************/
     
      for (iequiv=0; iequiv <= num_equivalences; iequiv++)
        {
          if (strcmp( equivalence_list[iequiv].type, atom_pot ) == 0) 
            {
              strcpy( atom_pot, equivalence_list[iequiv].stretch);
              break;
            }
        }      

/******************************************************************/
/******* Loop over this atom's neighbours *************************/
/******************************************************************/

      for (neigh_index=0; neigh_index <= p_atom->num_neigh; neigh_index++)
        {
          ineigh = p_atom->neighb[neigh_index];
          p_neigh = p_molecule+ineigh;
   
          strcpy( neigh_pot, p_neigh->pot);

/******************************************************************/
/******* Check equivalence table **********************************/
/******************************************************************/

          for (iequiv=0; iequiv <= num_equivalences; iequiv++)
            {
              if (strcmp( equivalence_list[iequiv].type, neigh_pot ) == 0)
                {
                  strcpy( neigh_pot, equivalence_list[iequiv].stretch);
                  break;
                }
            }

          found_pot_flag = FALSE;
          for (jpot=0; jpot<=num_stretches; jpot++)
            {

/******************************************************************/
/***** test if potential matches database potential ***************/
/******************************************************************/

           if ((   strcmp( intra_pair_potent[jpot].atom1, atom_pot  ) == 0 
                && strcmp( intra_pair_potent[jpot].atom2, neigh_pot ) == 0 )
             || (  strcmp( intra_pair_potent[jpot].atom2, atom_pot  ) == 0
                && strcmp( intra_pair_potent[jpot].atom1, neigh_pot ) == 0 ) )

             {
/**************** assign array reference **************************/

               p_atom->neighb_stretch_list[neigh_index] = jpot;
               found_pot_flag = TRUE;
               break;
             }

          } /* end of jpot loop */

      if (!found_pot_flag)
        {

/******* no  potential present for this element *******************/

          if (DEBUG)
           {
            printf("WARNING: Cannot assign intra stretch potential for Atoms %s (pot=%s)",
							p_atom->label,p_atom->pot);
            printf("       and %s (pot=%s)\nCan continue!\n", p_neigh->label,p_neigh->pot);
            }

/*            exit(EXIT_FAILURE); */
        }

     }/* end loop ineigh */

  }/* end loop iatom */

return;
}


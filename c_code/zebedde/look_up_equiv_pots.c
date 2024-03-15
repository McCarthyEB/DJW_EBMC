/*************************************************************************/
/* look_up_equiv_pots.c                                                  */
/* Look up the equivalences for a single atom and its list of            */
/* neighbours.                                                           */
/*                                                                       */
/* Started Dave Willock 21st August 2006                                 */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void look_up_equiv_pots(atom *p_molecule, atom *p_atom, char *p_atom_pot, pot_types *p_neigh_pot_list,
                        int which_type, int need_neighbs)
{
#include "header.h"
int ineigh;
int neigh_index, iequiv;
atom *p_neigh;

/******************************************************************/
/******* Check equivalence table for atom *************************/
/******************************************************************/
     
  strcpy(p_atom_pot, p_atom->pot);

  for (iequiv=0; iequiv <= num_equivalences; iequiv++)
    {
       if (strcmp( equivalence_list[iequiv].type, p_atom_pot ) == 0) 
        {
           if (which_type == EQUIV_NONBOND)  
                   strcpy( p_atom_pot, equivalence_list[iequiv].nonbond);

           if (which_type == EQUIV_STRETCH)
                   strcpy( p_atom_pot, equivalence_list[iequiv].stretch);
               
           if (which_type == EQUIV_ANGLE)
                   strcpy( p_atom_pot, equivalence_list[iequiv].angle);

           if (which_type == EQUIV_TORSION)  
                   strcpy( p_atom_pot, equivalence_list[iequiv].torsion);

           if (which_type == EQUIV_OOP)
                   strcpy( p_atom_pot, equivalence_list[iequiv].oop);
           break;
        }
    }      

/******************************************************************/
/******* Loop over this atom's neighbours *************************/
/******* and look them up in the equivalence table ****************/
/******************************************************************/

  if (need_neighbs)
    {
      for (neigh_index=0; neigh_index <= p_atom->num_neigh; neigh_index++)
        {
          ineigh = p_atom->neighb[neigh_index];
          p_neigh = p_molecule+ineigh;

          strcpy( p_neigh_pot_list->pot, p_neigh->pot);

/******************************************************************/
/******* Check equivalence table **********************************/
/******************************************************************/

          for (iequiv=0; iequiv <= num_equivalences; iequiv++)
            {
              if (strcmp( equivalence_list[iequiv].type, p_neigh_pot_list->pot ) == 0)
                {
                  if (which_type == EQUIV_NONBOND)  
                      strcpy(  p_neigh_pot_list->pot, equivalence_list[iequiv].nonbond);

                  if (which_type == EQUIV_STRETCH)
                      strcpy(  p_neigh_pot_list->pot, equivalence_list[iequiv].stretch);
               
                  if (which_type == EQUIV_ANGLE)
                      strcpy(  p_neigh_pot_list->pot, equivalence_list[iequiv].angle);

                  if (which_type == EQUIV_TORSION)  
                      strcpy(  p_neigh_pot_list->pot, equivalence_list[iequiv].torsion);

                  if (which_type == EQUIV_OOP)
                      strcpy(  p_neigh_pot_list->pot, equivalence_list[iequiv].oop);
                  break;
                }
            }
           p_neigh_pot_list++;
         }  /* end loop ineigh */
      }

return;
}


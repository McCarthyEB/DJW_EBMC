/*************************************************************************/
/* assign_hbond.c                                                        */
/* Sets up the hbond indexes for atoms when the AMBER forcefield is used */
/*                                                                       */
/* Started Dave Willock 4th March 1999                                   */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void assign_hbond(atom *p_molecule, int num_atoms)
{
#include "header.h"

int iatom,jpot;

BOOLEAN found_pot_flag;
atom *p_atom;

/******************************************************************/
/******* assign potential index to each atom in list **************/
/******************************************************************/

  for (iatom=0; iatom<= num_atoms; iatom++)
    {
      p_atom = p_molecule+iatom;
      p_atom->hb_list= -1;
      found_pot_flag = FALSE;

      for (jpot=0; jpot<=num_h_pot_types; jpot++)
        {

/******************************************************************/
/***** test if potential matches database potential ***************/
/******************************************************************/

           if (strcmp( h_pot_types[jpot].name, p_atom->pot ) == 0 )
             {
/**************** assign array reference **************************/

               p_atom->hb_list = jpot;
               found_pot_flag = TRUE;
               break;
             }
        }
    }
return;
}


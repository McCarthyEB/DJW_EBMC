/*************************************************************************/
/* assign_bond_inc.c                                                     */
/* Sets up the bond_inc indexes for atoms :  AMBER forcefield            */
/*                                                                       */
/* Started Dave Willock 14th April 1999                                  */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void assign_bond_inc(atom *p_molecule, int num_atoms)
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
      p_atom->bi_list= -1;
      found_pot_flag = FALSE;

      for (jpot=0; jpot<=num_bond_inc_types; jpot++)
        {

/******************************************************************/
/***** test if potential matches database potential ***************/
/******************************************************************/

           if (strcmp( bond_inc_types[jpot].name, p_atom->pot ) == 0 )
             {
/**************** assign array reference **************************/

               p_atom->bi_list = jpot;
               found_pot_flag = TRUE;
               break;
             }
        }
    }
return;
}


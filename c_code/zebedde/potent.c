
/*************************************************************************/
/* potent.c                                                               */
/* Sets up the non-bonding parameters and calculates intermediate values */
/* Mode of calculation determined by input file                          */
/*                                                                       */
/*                                                                       */
/* Started DWL 27/11/94                                                  */
/*************************************************************************/


/***** I rewrote this on 13/4/95 to take pointers ******/
/***** steric_setup assigns vdw radii*****/
/***** setup_nonbond calcs the intermediate bits for NB calcs ****/
/***** assign_nonbond assigns the pot type to a structure of atoms *****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void steric_setup(atom *p_molecule, int num_atoms)
{
#include "header.h"
atom *p_atom;

int i,j;
BOOLEAN found_atom_flag;

/****** Loop through the atom structure                **********/
/******     calculate the non-bonding parameters       **********/
/******     and assign van der Waals radii             **********/

/****** Loop through ***********************************/

for (i=0; i<num_atoms;i++)
	{

	found_atom_flag = FALSE;
	p_atom = p_molecule + i;

	for (j=0;j<NUM_ELEMENTS;j++)
		{
		if ((strcmp(p_atom->elem, period_table[j].elem)) == 0)
			{
			p_atom->vdw = period_table[j].vdw;
			found_atom_flag = TRUE;	
			break;
			}
		}
	/*** if we haven't found the right element then Serious Error ***/
	if (found_atom_flag == FALSE)
		{
		fprintf(output_fp,"ILLEGAL Element\n");
		fprintf(output_fp,"Atom number %i, element given as %s\n",
							 i,p_atom->elem);
		fprintf(output_fp,"Aborting......\n");
                fflush(output_fp);
                fflush(stdout);
		exit(1);
		}
	}

return;
}

void assign_nonbond(atom *p_molecule, int num_atoms)
{
#include "header.h"
int iatom,jpot, iequiv;
	
BOOLEAN found_pot_flag;
atom *p_atom;

/******************************************************************/
/******* assign potential index to each atom in list **************/	
/******************************************************************/

  for (iatom=0; iatom<= num_atoms; iatom++)
    {
      p_atom = p_molecule+iatom;	
      found_pot_flag = FALSE;

      for (jpot=0; jpot<=num_potential_types; jpot++)
        {

/******************************************************************/
/***** test if potential matches database potential ***************/
/******************************************************************/

           if (strcmp( potent[jpot].pot, p_atom->pot ) == 0 )
             {
/**************** assign array reference **************************/

               p_atom->nb_list = jpot;
               found_pot_flag = TRUE;
               break;
             }

        }

      if (!found_pot_flag)
        {

/******************************************************************/
/******** go through and try to match according to equivalences ***/
/******************************************************************/

              for (iequiv=0; iequiv <= num_equivalences; iequiv++)
                 {
                     if (strcmp( equivalence_list[iequiv].type, p_atom->pot ) == 0)
                        {
                            for (jpot=0; jpot<=num_potential_types; jpot++)
                               {
                                  if (strcmp( equivalence_list[iequiv].nonbond, potent[jpot].pot ) == 0)
                                    {
                                       p_atom->nb_list = jpot;
                                       found_pot_flag = TRUE;
                                       break;
                                    }
                               }
                        }
                   }
        }

      if (!found_pot_flag)
        {
/******* no  potential present for this element *******************/

          fprintf(output_fp,"ERROR: Cannot assign potential for Atom number %d %s (Elem=%s)\n",
						iatom,p_atom->label,p_atom->elem);
          fprintf(output_fp,"ERROR: Atom has potential set in file as >>%s<<\n", p_atom->pot);
          fflush(output_fp);
          fflush(stdout);
          exit(EXIT_FAILURE);
        }

    }/* end loop i*/

return;
}


/***************************************************************/
/* find_hydrogens.c : scan template and makr up all hydogens   */
/*                    i.e. free valencies                      */
/* started 22/3/94 Dave and Dewi                               */
/*                                                             */
/* Updated August 06 Dave Willock                              */
/*                    accept HA and HB as valid element types  */
/*                    for hydrogen to allow specific sites to  */
/*                    be used along with rules.                */
/*                                                             */
/* If HA and/or HB are present only accept these and reject    */
/* normal H atoms.                                             */
/*                                                             */
/***************************************************************/

#include <stdio.h>     
#include "global_values.h"
#include "maxima.h"
#include "structures.h"

void find_hydrogens(atom *p_molecule, int *p_hydrogen_list, 
                    int *p_num_hydrogens, int num_atoms,
                    int *p_have_AB )
{
int i;

atom *p_atom;

*p_num_hydrogens = -1;
*p_have_AB=FALSE;

p_atom=p_molecule;
for (i=0; i<= num_atoms; i++)
   {
	if (   (p_atom->elem[0] == 'H' && p_atom->elem[1] == 'A')
            || (p_atom->elem[0] == 'H' && p_atom->elem[1] == 'B')) 
           {
              *p_have_AB=TRUE;
           }
       p_atom++;
   }

p_atom=p_molecule;
if (*p_have_AB)
  {
     for (i=0; i<= num_atoms; i++)
	{
	if (   (p_atom->elem[0] == 'H' && p_atom->elem[1] == 'A')
            || (p_atom->elem[0] == 'H' && p_atom->elem[1] == 'B')) 
              {
		(*p_num_hydrogens)++;
		*(p_hydrogen_list+*p_num_hydrogens) = i;
              }
      p_atom++;
        }
  }
else
  {
     for (i=0; i<= num_atoms; i++)
	{
	if ( (p_atom->elem[0] == 'H' && p_atom->elem[1] == '\0') )
              {
		(*p_num_hydrogens)++;
		*(p_hydrogen_list+*p_num_hydrogens) = i;
              }
      p_atom++;
        }
  }
return;
}

		

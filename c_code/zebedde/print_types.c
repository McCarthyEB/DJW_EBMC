/**********************************************************************/
/* print_types.c : prints types of atoms in a molecule                */
/* NOTE: neighb is references from ZERO                               */
/* started Dave and Dewi 16/10/95                                     */
/**********************************************************************/


#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_types( atom_number *p_types, int num_types, 
                  int *p_hyd_list, int num_hyds, FILE *fp)
{
int ihyd, itypes;

/******** print out the number of each type of atom *************/

fprintf(fp, "\nThis molecule has %d types, formula: ", num_types+1);

for (itypes=0; itypes <= num_types; itypes++)
   {
     fprintf(fp, "%s%d ", p_types->atom_type,  p_types->num+1); 
     p_types++;
   }
fprintf(fp,"\n");

/******** print out the details of H positions ******************/

fprintf( fp,"\n%d modifiable atoms (hydrogens!).\n", num_hyds+1);

/* loop over hydrogens  */

fprintf(fp, "Hydrogen Indices: "); 

for (ihyd = 0; ihyd <= num_hyds; ihyd++)
  {
  if (ihyd == num_hyds)
    {
       fprintf( fp," %d.",*(p_hyd_list+ihyd)+1);
    }
  else
    {
       fprintf( fp," %d,",*(p_hyd_list+ihyd)+1);
    }
  } 

fprintf( fp,"\n");
return;
}

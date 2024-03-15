/*************************************************************************************/
/* print_neighbours_heavies.c : prints neighbour list for molecule ommitting H atoms */
/* NOTE: neighb is references from ZERO                                              */
/* Created for debugging Dave Willock, July 09                                       */
/*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_neighbours_heavies( atom *p_molecule, int num_atoms, FILE *fp)
{
int iatom1,iatom2;
int idummy1;

atom *p_atom1;

/* loop over atoms  */
printf("printing neighbour information\n");

for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
  p_atom1= p_molecule+iatom1;

  if ( strcmp(p_atom1->elem, "H") != 0 )
    {
      fprintf(fp, "Neighbours for atom %d %s, mol %d: ", iatom1, p_atom1->label, p_atom1->mol); 
      if (p_atom1->num_neigh < 0) fprintf(fp,"\n");

          for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
    		{
                iatom2= (p_atom1->neighb[idummy1]);
        
		fprintf(fp, " %d:",iatom2);
		fprintf(fp, "%s ",(p_molecule+iatom2)->label); 
		}
	  fprintf(fp,"\n");
     } 
  } 
  return;
}

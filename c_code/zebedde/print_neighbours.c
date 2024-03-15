/**********************************************************************/
/* print_neighbours.c : prints neighbour list for molecule            */
/* NOTE: neighb is references from ZERO                               */
/* started Dave and Dewi 23/3/95                                      */
/**********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp)
{
int iatom1,iatom2;
int idummy1;

atom *p_atom1;

/* loop over atoms  */
printf("printing neighbour information\n");

for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
  p_atom1= p_molecule+iatom1;
  fprintf(fp, "Neighbours for atom %d %s (elem %s), mol %d: ", iatom1, p_atom1->label, 
                                                               p_atom1->elem, p_atom1->mol); 
  if (p_atom1->num_neigh < 0) fprintf(fp,"\n");

      for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
		{
        iatom2= (p_atom1->neighb[idummy1]);
        
		fprintf(fp, " %d:",iatom2);
		fprintf(fp, "%s ",(p_molecule+iatom2)->label); 
		}
	  fprintf(fp,"\n");
  } 
  return;
}

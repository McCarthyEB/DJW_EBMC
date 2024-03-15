/**********************************************************************/
/* print_dist_mat.c : prints a distance matrix for all atoms to their */
/*                    neighbours.                                     */
/* NOTE: neighb is references from ZERO                               */
/* started Dave February 2006                                         */
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void min_image( double *x, double *y, double *z);

void print_dist_mat( atom *p_molecule, int num_atoms, int use_pbc, FILE *fp)
{
int iatom1,iatom2;
int idummy1;

double dx, dy, dz, dist;

atom *p_atom1, *p_atom2;

/* loop over atoms  */
printf("printing distance matrix for neighbours\n");

for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
  p_atom1= p_molecule+iatom1;
  fprintf(fp, "Neighbours for  %s: ",p_atom1->label); 
  if (p_atom1->num_neigh < 0) fprintf(fp,"\n");

      for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
		{
                iatom2= (p_atom1->neighb[idummy1]);

                p_atom2=p_molecule+iatom2; 

                dx = p_atom2->x - p_atom1->x;
                dy = p_atom2->y - p_atom1->y;
                dz = p_atom2->z - p_atom1->z;
        
                if (use_pbc) min_image(&dx, &dy, &dz); 

                dist= sqrt(dx*dx + dy*dy + dz*dz);
                
		fprintf(fp, " %d:",iatom2);
		fprintf(fp, "%s %10.6f ",(p_atom2)->label, dist); 
   
                if (dist > 3.0) fprintf(fp, " long? ");
		}
	  fprintf(fp,"\n");
  } 
  return;
}

/****************************************************************/
/* print_molecule.c : prints a molecule (.car format            */
/* Adapted to do D to H conversion Aug.98 DJW                   */
/****************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H)
{
int i;
atom *p_atom;
char elem[5];

if (D_to_H)
 {
     for (i=0;i<=num_atoms;i++)
      {
        p_atom = p_molecule+i;

        strcpy(elem, p_atom->elem);
        if (strcmp(elem,"D")==0) strcpy(elem, "H");

/*** Convert HA and HB to simple H atoms **/
        if (strcmp(elem,"HA")==0) strcpy(elem, "H");
        if (strcmp(elem,"HB")==0) strcpy(elem, "H");
   
        fprintf(output_fp,"%-5s %14.9f %14.9f %14.9f %-4s %-4s   %-3s     %-2s %6.3f\n" ,
                 p_atom->label, p_atom->x,
                 p_atom->y, p_atom->z,p_atom->group,p_atom->group_no,
                 p_atom->pot,elem,p_atom->part_chge);
      }

 }
else
 {
   for (i=0;i<=num_atoms;i++)
     {
	p_atom = p_molecule+i;

	fprintf(output_fp,"%-5s %14.9f %14.9f %14.9f %-4s %-4s   %-3s     %-2s %6.3f\n" ,
            p_atom->label, p_atom->x,
            p_atom->y, p_atom->z,p_atom->group,p_atom->group_no,
            p_atom->pot, p_atom->elem, p_atom->part_chge);
     }
 }
return;
}

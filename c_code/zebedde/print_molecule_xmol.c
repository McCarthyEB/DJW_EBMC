/****************************************************************/
/* print_molecule_xmol.c : prints a molecule in Xmol xyz format */ 
/*                         6/6/96 DWL                           */ 
/****************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_molecule_xmol(atom *p_molecule, int num_atoms, FILE *output_fp)
{
int i;
atom *p_atom;

for (i=0;i<=num_atoms;i++)
	{
	p_atom = p_molecule+i;
    fprintf(output_fp,"%-5s  %14.9f %14.9f %14.9f %6.3f\n",
			p_atom->elem, p_atom->x, p_atom->y, p_atom->z, p_atom->part_chge);
	}
return;
}

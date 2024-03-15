/*********************************************************************/
/* move_molecule.c ; re-position any molecule fragment               */
/*                   adapted to cope with multiple molecules in      */
/*                   list and only move that requested.              */
/*                   March 07 Dave Willock.                          */
/*********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void move_molecule(atom *p_molecule, int num_atoms, double *move_vec,
                   int which_mol)

{
int i; 
double *vec_comp;
atom *p_atom;

/*** if which_mol is negative shift entire list ***/
if (which_mol < 0)
    {
      p_atom = p_molecule;
      for (i=0;i<= num_atoms; i++)
	{
	vec_comp = move_vec;

	p_atom->x += *vec_comp;
	vec_comp++;
	p_atom->y += *vec_comp;
	vec_comp++;
	p_atom->z += *vec_comp;
        p_atom++;
	}
    }

/*** otherwise only move atoms belonging to the requested molecule ***/

else
    {
      p_atom = p_molecule;
      for (i=0;i<= num_atoms; i++)
	{
          if (p_atom->mol == which_mol)
            {
              vec_comp = move_vec;

	      p_atom->x += *vec_comp;
	      vec_comp++;
	      p_atom->y += *vec_comp;
	      vec_comp++;
	      p_atom->z += *vec_comp;
           }
        p_atom++;
	}
    }

return;
}

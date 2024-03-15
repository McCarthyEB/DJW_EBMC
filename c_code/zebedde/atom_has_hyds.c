/****************************************************************/
/* atom_has_hyds.c : look at a particular atom to see if it has */
/*                   hydrogens and so can be bonded to          */
/* started 10/2/96 Dave Willock                                 */
/****************************************************************/

#include <stdio.h>     
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

int atom_has_hyds(atom *p_molecule, int atom_num)
{
int ineigh;

atom *p_neigh, *p_atom;

p_atom= p_molecule + atom_num;

for (ineigh=0; ineigh <= (p_atom->num_neigh); ineigh++)
   {
     p_neigh=  p_molecule + p_atom->neighb[ineigh];
     if (p_neigh->elem[0] == 'H' && p_neigh->elem[1] == '\0') return TRUE;
   }
return FALSE;
}

		

/**********************************************************************/
/* fix_potentials.c      : for Discover run, reset potentials at      */
/*                       : link atoms                                 */
/* started Dewi 16/6/95                                               */
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int  fix_potentials_cff91( atom *p_molecule, int this_atom);

int fix_potentials( atom *p_molecule, int this_atom, int num_atoms)
{
#include "header.h"

return fix_potentials_cff91( p_molecule, this_atom);
}

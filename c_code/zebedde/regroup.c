/******************************************************************/
/* regroup.c :  puts all atoms into unique (Biosym) residues      */
/*              Each residue is call res%i where i is an int      */
/*              likewise residue group_number.                    */
/* dewi 19/6                                                      */
/******************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void regroup(atom *p_molecule, int num_to_group, int group_number)
{
int i;
atom *p_atom;
char dummy[10];
char dummy_no[10];

if (group_number < 10)
    sprintf(dummy,"FR0%i",group_number); /* residue name */
else
    sprintf(dummy,"FR%i",group_number); /* residue name */

sprintf(dummy_no,"%i",group_number); /* residue name */

for (i=0;i<=num_to_group; i++)
	{
	p_atom = p_molecule+i;
	strcpy( p_atom->group,dummy);
	strcpy( p_atom->group_no,dummy_no);
	}
return;
}

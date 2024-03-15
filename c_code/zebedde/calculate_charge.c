/**************************************************************************/
/* calculate_charge.c Estimates formal overall charge on molecule by      */
/*                    bond counting                                       */
/*started DWl 21/6/95                                                     */
/*                                                                        */
/* For now it cycles the molecule, finds N and bumps the charge by 1 if   */
/* it has 4 neighbours!                                                   */
/* This will only get more clever if it is needed!!!                      */
/**************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
 
int calculate_charge(atom *p_molecule, int num_atoms)
{
#include "header.h"
int i;
int charge =0;
atom *p_atom;

for (i=0;i<=num_atoms;i++)
	{
	p_atom = p_molecule+i;

	if ((strcmp(p_atom->elem,"N")==0) && (p_atom->num_neigh == 3)) charge++;

	}

return(charge);
}


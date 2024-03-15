/************************************************************/
/* shuffle_atom_list                                        */
/* shuffles the atom list to remove unwanted hydrogens after*/
/* bond formation                                           */
/*                                                          */
/* Started DJW 10/2/96                                      */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int shuffle_atom_list(atom *p_molecule, 
                      int num_atoms, int hyd_index1, int hyd_index2)
{
#include "header.h"
int num_new_atoms, iloop;

atom *p_atom;

num_new_atoms=-1;
p_atom= p_molecule-1;

if (hyd_index2 < hyd_index1)
  {
    printf("ERROR: entered shuffle_atom_list with index2 (%d)  greater than index1 (%d)...?\n",
                 hyd_index2, hyd_index1);
    exit(0);
  }
/***************************************************************************/
/**** Loop over all atoms in the molecule and throw out the two positions **/
/**** passed.                                                             **/
/***************************************************************************/

for (iloop=0; iloop <= num_atoms; iloop++)
   {
      p_atom++;
      if ((iloop != hyd_index1) && (iloop != hyd_index2))
        {
           num_new_atoms++;
           *(p_molecule+ num_new_atoms) = *p_atom ;
        }
   }

return num_new_atoms;
}

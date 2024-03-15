/******************************************************************************/
/**** Find_line_atoms finds the first occurance of the atoms in           *****/
/**** the template which match the elements defining the line DJW July 98 *****/
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
#include "header.h"
#include "global_values.h"

void find_line_atoms(int line, atom *p_mol, int num_atoms, int *p_iA, int *p_iB)
{
int iatom, got_A;

/******************************************************************************/
/**** Find the first occurance of the atoms ***********************************/
/******************************************************************************/

got_A= FALSE;
for (iatom=0; iatom <= num_atoms; iatom++)
   {
      if (!got_A && strcmp(p_mol->elem, line_atom[line].A) == 0) 
         {
            *p_iA = iatom;
            got_A = TRUE;
         }
      else if (got_A && strcmp(p_mol->elem, line_atom[line].B) == 0) 
         {
            *p_iB = iatom;
            break;
         }

      p_mol++;
   }

return;
}


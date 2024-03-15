/*************************************************************/
/* demarcate_molecule.c use molecular indices assigned to   **/
/* a list containing one or more molecules to shift the     **/
/* atoms so that the molecules are continuos. Use the       **/
/* demarc[] structure to record start, end and num for each **/
/* molecule subsection in the list.                         **/
/******                     started Oct 2010 Dave Willock   **/
/*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void demarcate_molecule(atom *p_molecule, int num_atoms, list_partition *p_demarc)
{
#include "header.h"

int swapped, iii, first; 

atom temp_atom;
atom *p_atom, *p_next_atom;
/**************************************/
swapped=TRUE;
first=0;
while (swapped)
  {
    swapped=FALSE;

/*** Note this loop runs to the next to last atom since ***/
/*** We compare the current atom with the next          ***/

    for (iii=first; iii < num_atoms; iii++)
      {
        p_atom=p_molecule+iii;
        p_next_atom=p_atom+1;
        
        if (p_atom->mol > p_next_atom->mol)
          {
            if (!swapped) first = iii+1;

            temp_atom=*p_atom;
            *p_atom=*p_next_atom;  
            *p_next_atom=temp_atom;  

            swapped=TRUE;
          }
      }
  }

/*** Should now have atoms arranged in blocks corresponding to molecules ****/

    p_atom=p_molecule;
    p_next_atom=p_atom+1;
    p_demarc->start=0;
    for (iii=0; iii < num_atoms; iii++)
      {
        if (p_atom->mol != p_next_atom->mol)
          {
             p_demarc->end = iii;
             p_demarc->num = p_demarc->end - p_demarc->start+1;

             if (iii+1 != num_atoms)
               {
                  p_demarc++; 
                  p_demarc->start = iii+1;
               }
          }
        p_atom++;
        p_next_atom++;
      }
/*** deal with last molecule ****/
p_demarc->end = num_atoms;
p_demarc->num = p_demarc->end - p_demarc->start+1;

return;
}

		

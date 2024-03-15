/***************************************************************/
/****** neighbour_order.c : To return the order of neighbour  **/
/******                     that atom2 is to atom1            **/
/******                                                       **/
/******     started Jan 96 Dave                               **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int neighbour_order(atom *p_molecule, int num_atoms, int atom1, int atom2)
{
#include "header.h"
int ineigh,  num_ends, atoms_at_ends[MAX_ENDS], seen[MAXTEMPLATE];
int this_neigh, iatom, iends, order, index;
int found, num_new_ends;

atom *p_current_atom, *p_neigh;

/*************************************************************/
/******** set atom1 as the first member of the ends list *****/
/*************************************************************/

num_new_ends= 0;
found= FALSE;
order= 0;
atoms_at_ends[0]= atom1;

for (iatom = 0; iatom <= num_atoms; iatom++)
   {
     seen[iatom]= FALSE;
   }
	 /* to look for ring then search thro atom1 as well  */
     /* seen[atom1]= TRUE; */

/*************************************************************/
/******** seearch outwards for atom 2 ************************/
/*************************************************************/

while (!found)
  {
     num_ends = num_new_ends;
     order++;

     for (iends=0; iends <= num_ends; iends++)
       {
          p_current_atom= p_molecule + atoms_at_ends[iends];

          for (ineigh = 0; ineigh <=  (p_current_atom->num_neigh); ineigh++)
            {
               this_neigh=  p_current_atom->neighb[ineigh];
               p_neigh= p_molecule + this_neigh;

/***************************************************************************/
/***** add this neighbour to the ends list if it is not :               ****/
/*****   a) atom1                                                       ****/
/*****   b) already been looked at                                      ****/
/*****   c) a hydrogen or dueterium: i.e. has no neighbours itself!     ****/
/***************************************************************************/

               if ( !seen[this_neigh] )
                  { 
                     if (this_neigh == atom2)
                       {
                         found = TRUE;
                         break;
                       }

                     seen[this_neigh]= TRUE;
                    
                     if (p_neigh->num_neigh >   0) 
                       {
                          num_new_ends++;
                          atoms_at_ends[num_new_ends]= this_neigh;
                       }
                  }
             }
          if (found) break;
        }

/******* shuffle atoms_at_ends list to cover the ones we have dealt with ****/

     if (!found)
       {
         index =-1;
         for (iends= num_ends+1; iends <= num_new_ends; iends++)
           {
              index++;
              atoms_at_ends[index] =  atoms_at_ends[iends];   
           }

/***************************************************************************/
/****** Altered 1st March 1996 to allow multi-molecular templates DJW   ****/
/****** If you run out of ends ignore other molecule!! return -ve order ****/
/***************************************************************************/

         if (index == -1)
           {
              printf("\nWarning: out of ends in neighbour_order!!!\n\n");
              return -1;
           }
         num_new_ends= index;
       }
  }

return order;
}

		

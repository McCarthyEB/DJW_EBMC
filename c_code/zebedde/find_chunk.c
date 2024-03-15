/***************************************************************/
/* find_chunk.c : scan template and flag all atoms in the     **/
/******           same chunk of molecule as atom2 but not     **/
/******           atom1                                       **/
/******     started May 95 Dave and Dewi                      **/
/****** Build the chunck so that all atoms are min-image w.r.t**/
/****** the neighbours that got them involved.                **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);

void min_image( double *x, double *y, double *z);

int find_chunk(atom *p_molecule, int num_atoms, int *p_flag_list,
                int atom1, int atom2)
{
#include "header.h"
int ineigh, iflag, num_ends, atoms_at_ends[MAX_ENDS];
int *p_flag;
int this_neigh, iatom, iends;
double vec[3];
atom *p_current_atom, *p_atom;
atom *p_atom1, *p_atom2, *p_neigh;

/********* zero the flags ************/

p_flag= p_flag_list; 

for (iflag=0; iflag<= num_atoms; iflag++)
    {
      *p_flag = FALSE;
      p_flag++;
    }

/******** set atom2 as the first member of the chunk *********/

p_flag =  p_flag_list+  atom2;
*p_flag = 1;

num_ends= 0;
atoms_at_ends[0]= atom2;

/******** Move atom2 to min-image w.r.t. atom1 ***********************/

p_atom1= p_molecule+atom1;
p_atom2= p_molecule+atom2;

join_atoms( p_atom1, p_atom2, &vec[0]);

/*** take min_image and shift atom2 to it ****************************/

min_image(&vec[0], &vec[1], &vec[2]);

p_atom2->x = p_atom1->x + vec[0];
p_atom2->y = p_atom1->y + vec[1];
p_atom2->z = p_atom1->z + vec[2];

/******* set the flags for all end atoms and generate new ends ******/
/******** till we have been everywhere in the chunk             ******/

while (num_ends >= 0)
  {
    p_current_atom= p_molecule + atoms_at_ends[0];

    for (ineigh = 0; ineigh <= (p_current_atom->num_neigh); ineigh++)
      {
        this_neigh=  p_current_atom->neighb[ineigh];
        p_flag=  p_flag_list+ this_neigh;

/***************************************************************************/
/*** Shift the neighbours to min_image with respect to this atom ***********/
/***************************************************************************/

        p_neigh= p_molecule+this_neigh;

        join_atoms( p_current_atom, p_neigh, &vec[0]);

/*** take min_image and shift atom2 to it ****************************/

        min_image(&vec[0], &vec[1], &vec[2]);

        p_neigh->x = p_current_atom->x + vec[0];
        p_neigh->y = p_current_atom->y + vec[1];
        p_neigh->z = p_current_atom->z + vec[2];

/***************************************************************************/
/***** If atom1 crops up as a neighbour of someone else in this ************/
/***** chunk then atoms 1 and 2 must be in a ring               ************/
/***************************************************************************/

        if ( this_neigh == atom1 && atoms_at_ends[0] != atom2)
           {
             return FALSE;
           }

/***************************************************************************/
/***** add this neighbour to the ends list if it is not already flagged ****/
/***** dissallow atom1 joining the list to stop the other chunk being   ****/
/***** crossed into                                                     ****/
/***************************************************************************/

        if ( ! *p_flag &&  this_neigh != atom1 )
          { 
            *p_flag= TRUE;
            if (num_ends < MAX_ENDS && (p_molecule+this_neigh)->num_neigh != 0) 
              {
                num_ends++;
                atoms_at_ends[num_ends]= this_neigh;       
//if (DEBUG)
//  {
//                printf("accepted num_ends= %d this_neigh= %d\n",
//                        num_ends,this_neigh);
//  }
              }
            else if (num_ends > MAX_ENDS)
              {
                printf("Run out of ends in find_chunk.c increase MAX_ENDS");
                exit(1);
              }
          }
      } 

/******* shuffle atoms_at_ends list to cover the one we have dealt with ****/

    for (iends= 0; iends < num_ends; iends++)
      {
          atoms_at_ends[iends] =  atoms_at_ends[iends+1];   
      }
    num_ends--;
  }

/************DEBUG*********************/
//if (DEBUG)
//   {
//     printf("\n\nChunk 2 consists of : \n");
//
//     p_flag = p_flag_list;
//     p_atom= p_molecule;
//
//     for (iatom=0; iatom <= num_atoms; iatom++)
//       {
//         if (*p_flag) 
//           {
//             printf(" %d %s\n", iatom, p_atom->label); 
//           }
//         p_flag++;
//         p_atom++;
//       }
//    }
/**************************************/

return TRUE;
}


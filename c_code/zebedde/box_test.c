/************************************************************/
/* bump_check.c                                             */
/* Checks for intramolecular bumps in a molecule            */
/*                                                          */
/* Parameters:                                              */
/*  struct atom                                             */
/* Returns:                                                 */
/*  flag if bumps                                           */
/*                                                          */
/* Started DW/DWL 19/4                                      */
/*                                                          */
/* Updated for multiple molecules March 2013, Dave Willock  */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

int box_test( double *p_box_limits, atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc)

{
#include "header.h"
int i,out_of_your_box, imol;
int x_fail,y_fail,z_fail;
atom *p_atom;

list_partition *p_demarc;

if (symm_set)
  {
     printf("ERROR: Arrived to test box limits with symmetry set???\n");
     printf("ERROR: Currently only periodic systems can use symmetry!\n");
     exit(0);
  }

out_of_your_box= FALSE;

p_demarc= p_guest_demarc;
for (imol=0; imol < num_guests; imol++)
  {
    p_atom= guest_ptrs[imol];

    for (i=0; i <= p_demarc->end; i++)
      {
        x_fail= p_atom->x < p_box_limits[0] || p_atom->x > p_box_limits[1]; 
        y_fail= p_atom->y < p_box_limits[2] || p_atom->y > p_box_limits[3]; 
        z_fail= p_atom->z < p_box_limits[4] || p_atom->z > p_box_limits[5]; 

/* DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG */
        if (x_fail) 
         {
          printf("Out if box in x: %10.6f not within %10.6f < x < %10.6f",
                              p_atom->x, p_box_limits[0], p_box_limits[1]);
         }
        if (y_fail) 
         {
          printf("Out if box in y: %10.6f not within %10.6f < y < %10.6f",
                              p_atom->y, p_box_limits[2], p_box_limits[3]);
         }
        if (z_fail) 
         {
          printf("Out if box in z: %10.6f not within %10.6f < z < %10.6f",
                              p_atom->z, p_box_limits[4], p_box_limits[5]);
         }
/* DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG */

      if ( x_fail || y_fail || z_fail) 
        {
          out_of_your_box= TRUE;
          break;
        }
      p_atom++;
	}
     if (out_of_your_box) break;
     p_demarc++;
   }

return(out_of_your_box);
}

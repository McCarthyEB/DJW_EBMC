/************************************************************/
/* bump_check_with_flags.c                                  */
/* Checks for intramolecular bumps in a molecule            */
/*                                                          */
/* Parameters:                                              */
/*  struct atom                                             */
/* Returns:                                                 */
/*  flag if bumps                                           */
/*                                                          */
/* Started DW/DWL 19/4                                      */
/************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void min_image( double *x, double *y, double *z);

int bump_check_with_flags(atom *guest_ptrs[], int num_guests,
                          list_partition *p_guest_demarc,
                          interaction_indices *p_molmol_list, int num_molmol_list,
                          int *p_flag_list) 
{
#include "header.h"
int i,j,do_this,*p_flag,*p_flag_other;
int imol, jmol, ilist, ind, jnd;
int bump,num_neigh, neigh[10],ineigh;
int jneigh;
double dx,dy,dz,r,v_sum;
atom *p_old_atom, *p_new_atom;

list_partition *p_i_demarc, *p_j_demarc;

bump = FALSE;

for (ilist=0; ilist < num_molmol_list; ilist++)
  {
     imol=p_molmol_list->imol;    
     jmol=p_molmol_list->jmol;    
     p_i_demarc=p_guest_demarc+imol;
     p_j_demarc=p_guest_demarc+jmol;    

     ind=p_molmol_list->ind;
     jnd=p_molmol_list->jnd;

     for (i=0; i <= p_i_demarc->end; i++)
       {
         p_flag= p_flag_list+i; 

         if (*p_flag)
           {  
             p_new_atom= guest_ptrs[ind]+i;

             num_neigh = p_new_atom->num_neigh;
             for (ineigh=0; ineigh <= num_neigh; ineigh++)
                                                neigh[ineigh]= p_new_atom->neighb[ineigh];

	     for (j=0; j< p_j_demarc->end; j++)
	       {
                 do_this= TRUE;
   	         p_old_atom = guest_ptrs[jnd] + j;

                 for (ineigh=0; ineigh <= num_neigh; ineigh++)
                    {
                       if (j == neigh[ineigh]) do_this= FALSE;
                       for (jneigh=0; jneigh <= p_old_atom->num_neigh; jneigh++)
                        {
                           if ( p_old_atom->neighb[jneigh] == neigh[ineigh]) do_this= FALSE;
                        }
                    }
                 p_flag_other= p_flag_list+j;

                 if (do_this && !(*p_flag_other))
                    {

/******* calc separation ******************************************/

                      dx = p_new_atom->x - p_old_atom->x;
                      dy = p_new_atom->y - p_old_atom->y;
                      dz = p_new_atom->z - p_old_atom->z;

/******* for periodic pores use minimum image convention **********/

                      if (pbc) min_image( &dx, &dy, &dz);

/******************************************************************/
                      r = dx*dx + dy*dy + dz*dz;
                      r = sqrt(r);

        /* calc sum vderWaals radii */

                      v_sum = p_new_atom->vdw + p_old_atom->vdw;
                      v_sum = v_sum * vdw_scale; /* mult by scaling factor */

                      if (r < v_sum)
                        {
                           printf("bump true for %d with %d\n", ind, jnd );
		           bump = TRUE;
		           break;
		        }
                  }
              if (bump == TRUE) break;
            }
	}
     if (bump == TRUE) break;
   }
  if (bump == TRUE) break;
 }
return(bump);
}

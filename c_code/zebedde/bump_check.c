/************************************************************/
/* bump_check.c                                             */
/* Checks for inter-molecular bumps in a molecule           */
/*                                                          */
/* Parameters:                                              */
/*  struct atom                                             */
/* Returns:                                                 */
/*  flag if bumps                                           */
/*                                                          */
/* Started DW/DWL 19/4                                      */
/*                                                          */
/* Updated March 2013 to check within molecules and between */
/* them. Dave Willock                                       */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void min_image( double *x, double *y, double *z);

int bump_check(atom *guest_ptrs[], int num_guests,
               list_partition *p_guest_demarc,
               interaction_indices *p_molmol_list,
               int num_molmol_list)
{
#include "header.h"
int i,j,do_this, ineigh, i2neigh, imol, jmol;
int index, iatom, jatom, ilow, ind, jnd; 
int is_neigh;
double dx,dy,dz,r,v_sum;
atom *p_iatom, *p_jatom, *p_neigh;

list_partition *p_demarc, *p_j_dem;

/**************************************************************/
/*** DEBUGGING April 2001, num_old_atoms still has index ******/
/*** DEBUGGING April 2001, including deleted H           ******/
/*** DEBUGGING April 2001, Dave Willock                  ******/
/**************************************************************/


//   printf("Bump checking %d guests\n", num_guests);

    p_demarc=p_guest_demarc;
    for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
      {
        if (symm_set) ind=imol*(num_symm_ops+2);
                                             else ind=imol;

/***** Look for intra-molecule bumps first **********/
        
        p_iatom=guest_ptrs[ind];
        for (iatom=0; iatom<= p_demarc->end; iatom++)   /* loop over first atom */
          {

            for (jatom=iatom+1; jatom<= p_demarc->end; jatom++)   /* loop over second atom */
              {
                 p_jatom=guest_ptrs[ind]+jatom;
/**** check jatom is not a neighbour or second neighbour of iatom ****/
                 is_neigh=FALSE;
                 for (ineigh=0; ineigh <= p_iatom->num_neigh; ineigh++)
                   {
                      if (p_iatom->neighb[ineigh] == jatom) 
                        {
                           is_neigh=TRUE;
                           break;
                        }
                      p_neigh= guest_ptrs[ind]+p_iatom->neighb[ineigh];
  
                      for (i2neigh=0; i2neigh <= p_neigh->num_neigh; i2neigh++)
                        {
                          if (p_neigh->neighb[i2neigh] == jatom) 
                            {
                               is_neigh=TRUE;
                               break;
                            }
                        }
                   }
               
                 if (!is_neigh)
                   { 
/**** Need to test for bumps..... ***/
/******* calc separation ******************************************/

                    dx = p_iatom->x - p_jatom->x;
                    dy = p_iatom->y - p_jatom->y;
                    dz = p_iatom->z - p_jatom->z;

/******* for periodic pores use minimum image convention **********/

                    if (pbc) min_image( &dx, &dy, &dz);

/******************************************************************/
                    r = dx*dx + dy*dy + dz*dz;
                    r = sqrt(r);

        /* calc sum vderWaals radii */

                    v_sum = p_iatom->vdw + p_jatom->vdw;
                    v_sum = v_sum * vdw_scale; /* mult by scaling factor */

                    if (r < v_sum)
                     {
//                       if (DEBUG)
//                         {
//                        printf("Intra Tested %s and %s separation = %10.6f allowed = %10.6f\n",
//                                    p_iatom->label, p_jatom->label,
//                                    r, v_sum); 
//                         }
//     printf("bump_check returning TRUE\n");
//    exit(0);
		       return(TRUE);
                     }
                   }
                  p_jatom++;
                } /* jatom loop ends */
               p_iatom++;
             } /* iatom loop ends */

//        printf("DEBUG>> Clear of intra bumps...checking for guest..guest bumps\n");

/***** Now look for inter-molecule bumps first **********/

        p_j_dem=p_demarc+1;
        for (jmol=imol+1; jmol< num_guests; jmol++)              /* loop over guest molecules */
          {
            if (symm_set) jnd=jmol*(num_symm_ops+2);
                                             else jnd=jmol;

//            printf("molecule %d with %d...\n",ind,jnd);

            p_iatom=guest_ptrs[ind]; 
            for (iatom=0; iatom<= p_demarc->end; iatom++)   /* loop over first guests atoms */
              {
                ilow=0;
                p_jatom=guest_ptrs[jnd]; 
                if (jnd == ind ) ilow=iatom+1;
                for (jatom=ilow; jatom<= p_j_dem->end; jatom++)
                  {
/**** Need to test for bumps..... ***/

/******* calc separation ******************************************/

                    dx = p_iatom->x - p_jatom->x;
                    dy = p_iatom->y - p_jatom->y;
                    dz = p_iatom->z - p_jatom->z;

/******* for periodic pores use minimum image convention **********/

                    if (pbc) min_image( &dx, &dy, &dz);

/******************************************************************/
                    r = dx*dx + dy*dy + dz*dz;
                    r = sqrt(r);

        /* calc sum vderWaals radii */

                    v_sum = p_iatom->vdw + p_jatom->vdw;
                    v_sum = v_sum * vdw_scale; /* mult by scaling factor */

                    if (r < v_sum)
                     {
//                       if (DEBUG)
//                         {
//                          printf("Tested %s and %s separation = %10.6f allowed = %10.6f\n",
//                                      p_iatom->label, p_jatom->label,
//                                      r, v_sum); 
//                         }
		       return(TRUE);
                     }
                  p_jatom++;
                } /* jatom loop ends */
               p_iatom++;
             } /* iatom loop ends */
           p_j_dem++;
         }  /* jmol loop ends */
      p_demarc++;
    } /* imol loop ends */

//  printf("bump_check returning FALSE\n");
//  exit(0);

return(FALSE);
}

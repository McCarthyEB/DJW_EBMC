/*******************************************************/
/***** routine to rotate section of molecule   *********/
/***** assumes first atom in the list is fixed *********/
/***** Dave Willock April 21st 1995            *********/
/***** Adapted for intra-molecular energy terms*********/
/***** Oct. 2006 Dave Willock                  *********/
/***** Twist Action restructured Feb. 2014     *********/
/***** Dave Willock                            *********/
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

/* required functions ------------------------------------*/

double real_random(int done);
 
void rotate_with_flags(atom *p_molecule, double *p_axis, double *p_origin,
                       double theta, int *p_flag_list, int num_template_atoms, 
                       int which_mol);

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm );

/*--------------------------------------------------------*/

void do_the_twist(atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc,
                  int *p_flag_list,
                  symm_ops *p_symm,
                  links *p_link_atoms,
                  double *p_origin, double *p_axis, double *p_theta,
                  int which_mol)
{
#include "header.h"

int isymm, index, index_image;

list_partition *p_demarc;

/****************************************************************/
/**** Start of Executable lines *********************************/
/****************************************************************/
 p_demarc=p_guest_demarc+which_mol;

if (DEBUG) printf("DEBUG>> Arrived in do_the_twist\n");

// printf("DEBUG>> which_test = %d\n", which_test);
// printf("DEBUG>> Gathering molecule consisting of %d atoms\n", p_demarc->num);


  if (symm_set) index=which_mol*(num_symm_ops+2);
                                         else index=which_mol;

 if (DEBUG) printf("twisting around bond by theta= %10.6f\n", *p_theta);

/**** which_mol no longer needed here as only one molecule is sent, Nov. 2013, Dave Willock ****/
          rotate_with_flags(guest_ptrs[index], p_axis, p_origin, *p_theta, 
                            p_flag_list, p_demarc->end, -1);
 
/****** Recalculate symmetry related atoms  ***************/

          if (symm_set)
            {
               for (isymm=0; isymm <= num_symm_ops; isymm++)
                 {
                    index_image= index+isymm+1;
//                  printf("        ....doing image %d with %d atoms to index %d\n",
//                                                    isymm, p_demarc->num, index_image);

                    do_symmetry( guest_ptrs[index], p_demarc->end,
                                 guest_ptrs[index_image], &symm[isymm] );
                 }
            }

return;
}

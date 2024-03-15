/**** Routine takes the basic molecules and uses them to ****/
/**** set up the symmetry images according to the        ****/
/**** defined set of symmetry operations.                ****/
/**** molecules are organised so that the symmetry       ****/
/**** images follow the corresponding basic molecule     ****/
/**** Routine only expects to be called if symmetry set  ****/
/**** Dave Willock February 2013                         ****/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm );

void do_symmetry_all(atom *guest_ptrs[], int num_guests, 
                     list_partition *p_guest_demarc, 
                     symm_ops *p_symm, int num_symm_ops)                      
 {
int index, imol, isymm;
int index_image;
int num_with_symm;

list_partition *p_demarc;

symm_ops *p_this_symm;

    printf("DEBUG>> Symmetrising Guests\n");
    p_demarc=p_guest_demarc;

    printf("Have %d molecules and %d symmetry operations\n", num_guests, num_symm_ops);

    for (imol=0; imol < num_guests; imol++)
      {

/** num_symm_ops is an upper array bond, the base molecule counts for 1 so add 2 **/

        index= imol*(num_symm_ops+2);

        printf("Base molecule.... %d has index %d\n", imol, index);

        p_this_symm= p_symm;
        for (isymm=0; isymm <= num_symm_ops; isymm++)
          {
             index_image= index+isymm+1;
             printf("        ....doing image %d with %d atoms to index %d\n", 
                                                    isymm, p_demarc->num, index_image);

             do_symmetry( guest_ptrs[index], p_demarc->end, 
                          guest_ptrs[index_image], p_this_symm ); 
             p_this_symm++;

          }

         printf("Done all symmetry ops for this molecule\n");

         p_demarc++;
      }
         printf("Done all molecules\n");

    return;
  }


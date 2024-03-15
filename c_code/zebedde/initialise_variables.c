/*******************************************************************/
/********* Routine to initialise variables that are used ***********/
/********* by each template and so need to be reset on   ***********/
/********* each new template startup.                    ***********/
/********* Dave and Dewi began October 1995              ***********/
/*******************************************************************/

#include <stdio.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void initialise_variables( int *frag_weights,
                           int *p_orig_frag_weights,
                           int *p_num_anime_frames,
                           int num_seed_mols,
                           int *guest_hyd_weights_ptrs[],
                           int *p_num_guest_hyds,
                           int *p_sum_hyd_weights,
			   int initial_hydrogen_weight)
{
#include "header.h"
  int i, imol;

*p_num_anime_frames= 0;

/********************************************************************/
/****** Reset fragment weight array (This may have been altered *****/
/****** by builder to abide by fragment conc limitations etc.)  *****/
/********************************************************************/

  for (i=0; i<number_of_fragments; i++)
    {
      frag_weights[i]= *p_orig_frag_weights;
      p_orig_frag_weights++;
    }

/********************************************************************/
/****** Reset template hydrogen weightings **************************/
/********************************************************************/
 
   for (imol=0; imol<num_seed_mols; imol++)
     {
       printf("In initialise_variables processing molecule %d\n",imol);
       printf("In initialise_variables weighting %d hydrogens\n",
                                                    *p_num_guest_hyds);
       printf("In initialise_variables initial_hydrogen_weight %d\n",
                                              initial_hydrogen_weight);
       *p_sum_hyd_weights=0;
       for (i=0; i<= *p_num_guest_hyds; i++)
          {
/***set hydrogen weight*****************************/
            guest_hyd_weights_ptrs[imol][i]= initial_hydrogen_weight;
            *p_sum_hyd_weights += initial_hydrogen_weight;
          }
        p_sum_hyd_weights++;
        p_num_guest_hyds;
     }
  return;
}




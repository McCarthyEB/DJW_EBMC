/*******************************************************************/
/********* Routine to enforce defaults for variables    ************/
/********* that the user couldn't be bothered to set in ************/
/********* the input file.                              ************/
/********* Dave and Dewi began October 1995             ************/
/*******************************************************************/

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void enforce_defaults(double *p_num_rot_steps,
		      int *p_num_actions, int *frag_weights)
{ 

#include "header.h"
int iweight;
double pi;

/*******************************************************************/
/**** DEFAULT values are set in header.h and all variables *********/
/**** that have a default value should be set to something *********/
/**** rediculous in setup_defaults.c (usually -1 for ints!)*********/
/**** so that their set state can be tested here.          *********/
/*******************************************************************/

/*******************************************************************/
/**** Since non_bond and energy are now valid stop_ctf can be ******/
/**** negative!! Dave Willock March 1997                      ******/
/*******************************************************************/
/****  if (stop_ctf < 0) stop_ctf = STOP_CUTOFF_DEFAULT;   *****/

  if (ch_ctf < 0) ch_ctf = CH_CUTOFF_DEFAULT;

  if (nb_ctf < 0) 
    {
       nb_ctf = NB_CUTOFF_DEFAULT;
       nb_ctf_2 = nb_ctf * nb_ctf;
    }

  if (ring_ctf < 0) 
    {
       ring_ctf = RING_CUTOFF_DEFAULT;
       ring_ctf_2 = ring_ctf * ring_ctf;
    }

/**** Catch build weight as a signal this is just MC ****/
  
      fprintf(output_fp, "testing for MC prob_test %10.6f\n", prob_test);
  if (action_weights[0] == 0 && prob_test == 1.0 )
    {
      fprintf(output_fp, "Build weight set to 0, so assume pure MC run.\n");
      fprintf(output_fp, "Setting default action attempts to 1.\n");
      fprintf(output_fp, "Also van der Waals scale to 1.0\n");
      if (vdw_scale < 0) vdw_scale = 1.0;
      if (num_modify_attempts < 0) num_modify_attempts = 1;
      if (num_shake_attempts < 0) num_shake_attempts = 1;
      if (num_rock_attempts < 0) num_rock_attempts = 1;
       
    }


  if (vdw_scale < 0) vdw_scale = VDW_SCALE_DEFAULT;
  if (num_modify_attempts < 0) num_modify_attempts = MODIFY_TRY_DEFAULT;
  if (num_shake_attempts < 0) num_shake_attempts = ATTEMPTS_DEFAULT;
  if (num_rock_attempts < 0) num_rock_attempts = ATTEMPTS_DEFAULT;
  if (max_shake_step < 0) max_shake_step = SHAKE_STEP_DEFAULT;
  if (max_rock_step < 0) max_rock_step = ROCK_STEP_DEFAULT;
  if (box_fraction <= 0.0) box_fraction = BOX_FRACTION_DEFAULT;

/*******************************************************************/
/***** Pick up rediculous values that remain ***********************/
/*******************************************************************/

  if (box_fraction > 1.0)
        {
        fprintf(output_fp, "WARNING : box_fraction > 1, setting to 1\n");
        box_fraction = 1.0;
        }


/*******************************************************************/
/***** Set standard variables **************************************/
/*******************************************************************/

  pi= acos(-1.0);

/*******************************************************************/
/***** Apply scaling factors ***************************************/
/*******************************************************************/

   max_rock_step= max_rock_step/RAD_TO_DEG;
   *p_num_rot_steps= 2.0*pi/max_rock_step;

/*******************************************************************/
/***** If no user weights are given just weight everything *********/
/***** the same                                            *********/
/*******************************************************************/

/*******************************************************************/
/***** Fragment weights first **************************************/
/*******************************************************************/

  if (!frag_weights_given)
    {
      num_frag_weights= number_of_fragments;
      sum_frag_weights= 0;
      for (iweight = 0; iweight < num_frag_weights; iweight++)
         {
            frag_weights[iweight]= 1;
            sum_frag_weights++;
         }
    }
  if (number_of_fragments != num_frag_weights)
    {
    fprintf(output_fp, "Number of Weights does not match number of fragments\n\n");
    fprintf(output_fp, "Number of Weights = %d, Number of fragments = %d\n\n",
                                             num_frag_weights, number_of_fragments);
    fprintf(output_fp, "Exiting....\n");
    exit(1);
    }

/*******************************************************************/
/***** Now action weights ******************************************/
/*******************************************************************/


 if (!action_weights_given)
    {
      fprintf(output_fp, "No user supplied weights so using uniform weighting\n");

      num_action_weights= *p_num_actions;
      sum_action_weights= 0;
      for (iweight = 0; iweight < num_action_weights; iweight++)
         {
            action_weights[iweight]= 1;
            sum_action_weights++;
         }
    }

  if (*p_num_actions != num_action_weights)
    {
    fprintf(output_fp, "Number of Weights does not match number of actions\n\n");
    fprintf(output_fp, "Number of actions = %d, number of weights = %d\n",
                                          *p_num_actions, num_action_weights);
    fprintf(output_fp, "Exiting....\n");
    exit(1);
    }


  return;
}

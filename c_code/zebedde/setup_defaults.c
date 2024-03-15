/*******************************************************************/
/********* Routine to initialise variables that are the ************/
/********* same throughout a run, and to set up DEFAULT ************/
/********* values for variables                         ************/
/********* Dave and Dewi began October 1995             ************/
/*******************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "own_maths.h"
#include "ewald.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

double real_random(int done);

void setup_defaults(time_t *p_start_time, int *p_test_symm,
                    int *p_which_test, int *p_just_position, int *p_num_actions,
                    int *p_seed_type, int *p_want_dock,  int *p_need_monitors, 
				char *p_animation_file)
{ 
#include "header.h"
double rdummy, e_charge, epsilon, avos_number;
double arg;

/*************** time at start of run ***************************/

  *p_start_time = time(NULL); 

/*** Allow trial on symmetry without centre of mass align *******/

  *p_test_symm= FALSE;
  symm_set = FALSE;

/*** Turn off debugging output by default ***********************/

  DEBUG= FALSE;

/*** Turn off insight logfile generation ************************/

  logfile_needed= FALSE;

/*** Set user defined variables defaults ************************/

  prob_test= PROB_TEST_DEFAULT;

  max_templates = 1;
  num_forbidden_bonds= -1;
  num_conc_limits= -1;
  frag_weights_given=0;
  action_weights_given=0;
  num_pore_atoms = 0;

  *p_seed_type = -1; /****no seed supplied yet!*****/


   /* Now done in main since user can control seeding on */
   /* time or "1"                                        */
   /* DWL 21.10.96                                       */
   /* rdummy= real_random(0);                            */

/***** Logical to say we have imposed atom or bond limits *******/

   have_conc_limits= FALSE;

/***** number of allowed torsions, set to -ve to use as flag  ****/
/***** Dave Willock Oct. 2006                                 ****/
   num_allowed_torsions=-1;

/***** Logical to say we have imposed dihedral angles for new ****/
/***** bonds Dave Willock July 98                             ****/

   force_dihedrals= FALSE;

/***** Logical to say we have imposed line restraint          ****/
/***** Dave Willock July 98                                   ****/

   line_restraints= FALSE;
   line_hold= FALSE;
   have_tethers= FALSE;

/***** Logical to say we will try and fix atoms in minisation ****/
/***** Dave Willock July 98                                   ****/

   mini_fix_atoms = FALSE;

/****************************************************************/
/*** Logical to allow just position option, don't normally ******/
/*** do this!! Dave Willock April 1997                     ******/
/****************************************************************/

  *p_just_position= FALSE;
  *p_want_dock = FALSE;

/***** which_test controls what the cost function is ************/
  *p_which_test = STERIC_COST;

/******** Default values fo cost function energy evaluation ****/
  steric = TRUE;
  non_bonded = FALSE;
  charges = FALSE;

/******** default file name for defaults file       ***********/
  strcpy(defaults_file, NOT_SET);

/********default value for output type          ****************/
  verbose = TRUE;

/********default value for centring of seed     ****************/
  centralise_template = TRUE;

/********default value for optimising post mc run *************/
  optimise_post_mc = FALSE; 

/********default do not switch on monitors        *************/
  *p_need_monitors = FALSE;

/******** dummy values for cutoff for default check *************/

  ring_ctf = -1;
  stop_ctf = -1;
  ch_ctf = -1;
  nb_ctf = -1;
  vdw_scale = -1;
  num_modify_attempts= -1;
  num_shake_attempts= -1;
  num_rock_attempts= -1;
  max_shake_step= -1;
  max_rock_step= -1;

/********* Values for Warning reports, added July 09 Dave Willock ***/
/********* told_of_clash controls reporting of atom close calls   ***/

num_angle_warnings = -1;
told_of_clash = FALSE;

/********* Real value for AMBER h-bond cut off *****************/

 hb_ctf=2.5;
 hb_ctf_2= hb_ctf*hb_ctf;

/********* System temperature **********************************/

 temperature = 300.0;

/********* When carrying out optimisation on MC run ***********/

 mc_temp_step = 50.0;
 mc_modi_step = 1000;

/********* Default is to not produce final gulp input - same with car*********/

 want_final_gulp = 0; 
 want_final_car = 0;
/******** Default is not to use condor pool *******************/

num_dock = 50; /*Max number of docked structures to produce*/
finished_dock = 0;
want_condor = 0;
car_with_pore = 0;
exxon = 0; /*Default is not to run on exxons cluster*/

/********* dummy values for box variable ***********************/

  box_fraction = -1.0;
  user_box = FALSE;
  box_limits[0] =  99999.99;
  box_limits[1] = -99999.99;
  box_limits[2] =  99999.99;
  box_limits[3] = -99999.99;
  box_limits[4] =  99999.99;
  box_limits[5] = -99999.99;

/******** dummy file names for intermediate files *************/

  strcpy(template_strategy_file,   NOT_SET) ;
  strcpy(inpore_strategy_file,   NOT_SET) ;

/******** default paths for minimizers ************************/

  strcpy(mopac_path, DEFAULT_MOPAC_PATH);
  strcpy(discover_path, DEFAULT_DISCOVER_PATH);
  strcpy(gulp_path, DEFAULT_GULP_PATH);

/******** default command lines for MOPAC *********************/
  strcpy(mopac_cmdline_molecule, NOT_SET);
  strcpy(mopac_cmdline_inpore, NOT_SET);

/******** default command lines for GULP *********************/
  strcpy(gulp_cmdline_molecule, NOT_SET);
  strcpy(gulp_cmdline_inpore, NOT_SET); 

/******** default file names for intermediate files ***********/

  strcpy(minimizer_name, MINIMIZER_DEFAULT);
  strcpy(discover_forcefield_name, DISCOVER_FORCEFIELD_DEFAULT);
  strcpy(inpore_min_car, INPORE_CAR_DEFAULT);
  strcpy(inpore_min_mdf, INPORE_MDF_DEFAULT);
  strcpy(template_min_car, TEMPLATE_CAR_DEFAULT);
  strcpy(template_min_mdf, TEMPLATE_MDF_DEFAULT);
  strcpy(mopac_root, DEFAULT_MOPAC_OUTPUT);
  strcpy(gulp_root, DEFAULT_GULP_OUTPUT);

  *p_num_actions= NUMBER_OF_ACTIONS;

/********default values for initial minimisation****************/
  initial_minimize_inpore   = TRUE;
  initial_minimize_template = TRUE;

/********default values for initial minimisation****************/
  animate_flag = DEFAULT_ANIMATION_TYPE;
  strcpy(p_animation_file, DEFAULT_ANIMATION_FILE);


/********Mathematical and Physical constants *******************/
 
  one_sixth= 1.0/6.0;
  one_third= 1.0/3.0;
  
  pi = 2.0*acos(0);


  pi_b2 = pi/2;
  two_pi= 2.0*pi;
  four_pi= 4.0*pi;
  four_pi_sqrd= two_pi*two_pi;
  pi_tothehalf = sqrt(pi);

  erf_normalise = 2.0/pi_tothehalf;
  erf_accuracy= 0.0000000001;

/********Ewald sum parameters***********************************/

  ewald_accuracy= 0.00000001;
  arg= log(ewald_accuracy);

  if (arg > 0)
    {
      printf("ERROR : Invalid Ewald sum accuracy parameter: %16.12f\n",ewald_accuracy);
      exit(EXIT_FAILURE);
    }
 
  f_param = sqrt(-arg);

  epsilon = 8.85416E-12;
  e_charge = 1.602192E-19;
  avos_number= 6.022045E23;
  coul_prefactor= avos_number * e_charge * e_charge / (4.0 * pi * epsilon * 1E-10);

/******** To give energies in kcal/mol when working with e charges and angstroms ****/

  coul_prefactor= coul_prefactor / 4184.0;

  printf("Conversion factor to kcal = %10.6f\n",coul_prefactor);

  return;
}

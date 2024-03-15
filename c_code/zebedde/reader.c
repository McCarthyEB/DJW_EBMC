/*****************************************************************************/ 
/******  reader.c                                                *************/ 
/******  Input file reader. started 23/11/94 DWL                 *************/ 
/******  from Tim Bush's reader.c                                *************/ 
/******                                                          *************/ 
/******  22/2/96 DWL new keywords                                  ***********/
/******  initial_minimize                                          ***********/
/******  flag whether to do template min and inpore min before we  ***********/
/******                                                            ***********/
/******  21/10/96 DWL                                              ***********/
/******  Logical flagging now included for most things:            ***********/
/******     keyword [on/off]                                       ***********/
/******        replacing earlier effort with yes/no which was cack ***********/
/******  Keywords implemented in this way are:                     ***********/
/******  steric, non_bonded, charges, initial_minimize             ***********/
/******  centre, check                                             ***********/
/******  NEW KEYWORD : random [on/off]                             ***********/
/******  NEW KEYWORD : determines if real_random is random or not! ***********/
/******  21/10/96 DWL                                              ***********/
/******  NEW KEYWORD : cost_function [steric/non_bond/energy]      ***********/
/******  NEW KEYWORD : determines which cost_function to use       ***********/
/******  7th March 1997 Dave Willock                               ***********/
/******  Decided that cost_function directive can set steric,      ***********/ 
/******  non_bond and charges logicals now so the explicit input   ***********/ 
/******  to set them has been removed!!                            ***********/
/******  Dewi and Dave June 97                                     ***********/
/******  NEW KEYWORDS:                                             ***********/
/******  dihedrals <Num to set>                                    ***********/
/******  elem_A elem_B elem_C elem_D phi                           ***********/
/******  elem_A elem_B elem_C elem_D phi                           ***********/
/****** when a new bond is formed between B and C the dihedral     ***********/
/****** ABCD is set to phi. If the dihedral is not possible, i.e.  ***********/
/****** B(C) has no A(D) element neighbour nothing is done         ***********/ 
/****** 22/7/98 DJW                                                ***********/
/****** NEW KEYWORD structure DJW Aug. 1998                        ***********/
/****** minimizer constraint atoms <num_to_fix>                    ***********/
/****** N ( H D C )                                                ***********/
/****** C ( C C N H )                                              ***********/
/****** C ( C O N )                                                ***********/
/****** N ( C D C )                                                ***********/
/****** C ( C O H )                                                ***********/
/******                                                            ***********/
/******                                                            ***********/
/****** DWL Mods 25/8 Make action_weights require an end as well   ***********/
/******               just like frag_weights does already          ***********/
/******               Leave hydr_weights as is. end will be ignored?**********/
/******                                                             **********/
/****** 2006 edition                                                **********/
/****** NEW KEYWORD:                                                **********/
/****** slab : flag up that the host is a slab                      **********/
/****** poly : flag up that we are building polymers.               **********/
/******                                                             **********/
/****** Alterations May 2009 DJW:                                   **********/
/****** Additional flag to tok_get to allow verbatum reading i.e.   **********/
/****** without tolower being applied.                              **********/
/****** tok_get(input_fp, skip, tolower)                            **********/
/******                                                             **********/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <errno.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "reader.h"

int  parse_line(char *p_animation_file, int *p_seed_type,
                atom *pore, atom *frag_lib,
                int *frag_weights,
                int *seed_weights,
                int *p_num_of_seeds,
                int *p_seeds_to_use, 
                int *p_read_all_frames_as_seeds,
                int *p_random_seed,
                int *p_number_of_frames_in_arc,
		 		 int *p_initial_hydrogen_weight,
		 		 int *p_increment_hydrogen_weight,
                int *p_test_symm,
                char *p_analyse_name,
                int *p_peek_freq,
                int *p_randomisation_method,
                int *p_which_test,
                int *p_just_position, 
                int *p_is_slab, 
                int *p_is_poly,
                int *p_want_dock,
                int *p_is_empty_pore,
                int *p_need_monitors,
                monit *p_monitored,
                int just_count);

int locate_string( char *p_key, char *p_char, int num_of_chars );

int reader(char *file_to_read, char *p_animation_file,
           atom *pore, atom *frag_lib,
           int *frag_weights,
           int *seed_weights,
		    int *p_seed_type, int *p_num_of_seeds, int *p_seeds_to_use,
           int *p_read_all_frames_as_seeds, int *p_random_seed,
           int *p_number_of_frames_in_arc,
		    int *p_initial_hydrogen_weight,
           int *p_increment_hydrogen_weight,
		    int *p_test_symm,
           char *p_analyse_name,
		    int *p_peek_freq,
           int *p_randomisation_method,
           int *p_which_test,
           int *p_just_position,
           int *p_is_slab,
           int *p_is_poly,
           int *p_want_dock,
           int *p_is_empty_pore,
           int *p_need_monitors,
           monit *p_monitored,
           int just_count)
{
#include "header.h"
  int abandon_ship;

  if (!(input_fp = fopen(file_to_read, "r")))
     {
       fprintf(output_fp,"ERROR: Problem opening file %s for input\n", file_to_read);
       exit(1);
     }

  for(;;)   /* get the info from the file */
    {
      abandon_ship = parse_line(p_animation_file, p_seed_type, pore, frag_lib, frag_weights,
                                seed_weights, p_num_of_seeds, p_seeds_to_use,
                                p_read_all_frames_as_seeds,
                                p_random_seed, p_number_of_frames_in_arc,
                                p_initial_hydrogen_weight, 
                                p_increment_hydrogen_weight,
                                p_test_symm, p_analyse_name, p_peek_freq,
                                p_randomisation_method,
                                p_which_test,
                                p_just_position, p_is_slab,
                                p_is_poly, p_want_dock, p_is_empty_pore,
                                p_need_monitors, p_monitored, just_count);

      fclose(input_fp);
      return(abandon_ship);
    }
}

/*****************************************************************************
  Parse_line, chops an input line from the input file and chops it into
  tokens, each of which is tested against the list of allowed directives
  in reader.h.  The correct function to read subsequent program parameters
  is then called..
*****************************************************************************/

int check_element(char *p_element);
int find_kind(char *token, int level);
char * tok_get(FILE *input_fp, int skip_lines, int to_lower_case);
void get_title(void);
int get_pore_file(atom *pore, int just_count);
int get_gulp_pots_file(void);
int get_animation(char *p_animation_file);
int get_max_templates(void);
int get_vdw_scale(void);
int get_prob_test(void);
int get_temperature(void);
void get_inpore_strategy_name(void);
void get_template_strategy_name(void);
/* void get_selection_criteria(void); */
int get_discover_forcefield_name(void);



void get_fragment_library(atom *frag_lib, int just_count);
void get_template_library(void);
int get_seed_file(int *p_seed_type);
int get_seed_archive(int *p_seed_type, int *p_num_archive_frames_as_seeds,
                     int *p_seed_frames_to_read,
                     int *p_read_all_frames_as_seeds,
                     int *p_number_of_frames_in_arc);
int get_seed_fragment(int *p_seed_type, int *p_num_fragments_as_seeds,
                      int *p_fragments_for_seeds, int *p_random_seed);
void parse_error(char *last_tok);
int get_nb_cutoff(void);
int get_ring_cutoff(void);
int get_ch_cutoff(void);
void get_frag_weight(int *frag_weights, int just_count);
void get_seed_weight(int *seed_weight, int just_count);
void get_hydrogen_weight(int *p_initial_hydrogen_weight,
                         int *p_increment_hydrogen_weigot);
void get_action_weight(void);
void get_modify_attempts(void);
void get_rock_attempts(void);
void get_shake_attempts(void);
void get_rock_step(void);
void get_rock_fix_atom(void);
void get_shake_step(void);
int get_stop_cutoff(void);
int get_box_limits(void);
int get_box_fraction(void);
int get_symmetry(void);
void get_forced_dihedrals(void);
void get_allowed_torsions(void);
void get_atoms_to_fix(void);
void get_restraint_elements(void);
int get_tether_labels(void);
int get_minimizer_name(char *tok);
void get_mopac_cmds(int whichone);
void get_mopac_name(void);
void get_concentration_limits(void);
void get_forbidden_bonds(void);
void get_mopac_path(void);
void get_discover_path(void);
void get_forcefield_library(void);
int get_analyse_filename(char *p_analyse_name);
void get_peek_freq(int *p_peek_freq);

void get_mc_opt_actions(void);
void get_mctempstep(void);


void get_defaults_file(void);
void get_gulp_name(void);
void get_gulp_path(void);
void get_gulp_extra(void);
void get_gulp_cmds(int whichone);
int get_dock_energy(void);
int get_num_dock(void);

void print_biosym_header(FILE *file_fp,  int pbc_flag);

int    parse_line(char *p_animation_file, int *p_seed_type,
                  atom *pore, atom *frag_lib,
                  int *frag_weights,
                  int *seed_weights,
                  int *p_num_of_seeds,
                  int *p_seeds_to_use, 
                  int *p_read_all_frames_as_seeds,
                  int *p_random_seed,
                  int *p_number_of_frames_in_arc,
                  int *p_initial_hydrogen_weight,
                  int *p_increment_hydrogen_weight,
                  int *p_test_symm,
                  char *p_analyse_name,
                  int  *p_peek_freq,
                  int *p_randomisation_method,
                  int *p_which_test,
                  int *p_just_position,
                  int *p_is_slab,
                  int *p_is_poly,
                  int *p_want_dock,
                  int *p_is_empty_pore,
                  int *p_need_monitors,
                  monit *p_monitored,
                  int just_count)
{
#include "header.h"
  int   ans,l;
  int   bail_out;
  int   token, second_token, third_token;
  char  *tok, *last_tok;

  l = 0;
  bail_out= FALSE;
  

  for (;;)
    {

      if (!(tok = tok_get(input_fp, TRUE, TRUE)))
    		 {
		     if (++l > 3) break; else continue;
		     }
      l = 0;
/*****************************************************************/
/**** If required print out the token you are about to process ***/
/*****************************************************************/

      token = find_kind(tok, PRIME_DIRECTIVE);

/****DEBUG ****/
/*DEBUG=TRUE;*/
if (DEBUG)
 {
    printf("Token = >>%s<< find_kind assigns %d\n", tok, token);
 }
DEBUG=FALSE;
/****DEBUG ****/

/*****************************************************************/
/**** If the last command was not found let the user know ********/
/**** and then bail out after processing the entire file! ********/
/**** Added 4/3/97 Dave Willock **********************************/
/*****************************************************************/
      if (token == BLANK_DIRECT)
         {
            fprintf(output_fp,"ERROR: The command >>%s<< was not recognised\n",tok);
            bail_out= TRUE;
         }

      last_tok = tok;
      switch (token)
		    {
           case CENTRE_TEMPLATE: 
                tok = tok_get(input_fp, FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);

                last_tok =tok;
                switch (second_token)
                        {
                        case ON2  :centralise_template = TRUE;  break;
                        case OFF2 :centralise_template = FALSE; break;
                        case YES2 :centralise_template = TRUE;  break;
                        case NO2  :centralise_template = FALSE; break;
                        case DOCK : *p_want_dock=TRUE, centralise_template = TRUE; break;

                        case -1: fprintf(output_fp,"\nERROR: centre directive (cent) supplied but unrecognised ");
                                 fprintf(output_fp,"option given: >>%s<<.\n", tok);
                                 bail_out= TRUE;
                                 break;

                        case BLANK_DIRECT : fprintf(output_fp,"\nERROR: centre directive (cent) supplied but no option");
                                            fprintf(output_fp," given. Please add 'on' or 'off' as required.\n"); 
                                            bail_out= TRUE;
                                            break;
		 		 		 }
		 		 break;

		    case TITLE: get_title(); break;

		    case NB_CUTOFF: bail_out= get_nb_cutoff() || bail_out; 
                           break;

		    case RING_CUTOFF: bail_out= get_ring_cutoff() || bail_out; 
                             break;

		    case CH_CUTOFF: bail_out= get_ch_cutoff() || bail_out; 
                           break;

		    case PORE_FILE: ans= get_pore_file(pore, just_count); 
                           
                           if (ans == -10) 
                             {
                               *p_is_empty_pore = TRUE;
                             }
                           else
                             {
                               bail_out = ans || bail_out;
                             }
                           break;

		    case POTS_FILE: bail_out= get_gulp_pots_file() || bail_out;
		 		 		    break;

/*******************************************************************/
/** 12/6/01 convention flag used to set alignment of cell vectors **/
/**         with cartessian axis system                           **/
/*******************************************************************/

           case CONVENTION: 
                tok = tok_get(input_fp, FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);

                last_tok =tok;
                switch (second_token)
                        {
                        case XYZ  : use_xyz =  TRUE; use_zyx = FALSE; break;
                        case ZYX  : use_xyz = FALSE; use_zyx = TRUE ; break;

                        case -1: fprintf(output_fp,"\nERROR: convention directive supplied but unrecognised ");
                                 fprintf(output_fp,"option given: >>%s<<.\n", tok);
                                 fprintf(output_fp,"Please add 'xyz' or 'zyx' as required.\n"); 
                                 bail_out= TRUE;
                                 break;

                        case BLANK_DIRECT : fprintf(output_fp,"\nERROR: convention directive supplied but no option");
                                            fprintf(output_fp," given. Please add 'xyz' or 'zyx' as required.\n"); 
                                            bail_out= TRUE;
                                            break;
		 		 		 }
		 		 break;


/*******************************************************************/
/*******6.6.96 Change to take either biosym or xmol style files ****/
/*******6.6.96 no error checking yet  */ 
/*******18.7.08 AJWL added 'off' to restrict output during dock run*/
/*******************************************************************/

		    case ANIMATION: tok = tok_get(input_fp, FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);

                last_tok =tok;
                switch (second_token)
                  {
                     case BIOSYM: animate_flag = BIOSYM_ANIMATION; break;
                     case PXYZ  : animate_flag = XMOL_ANIMATION_WITHPORE; break;
                     case TXYZ  : animate_flag = XMOL_ANIMATION_NOPORE;   break;
                     case OFF2  : animate_flag = FALSE; break;
		 		      case -1    : bail_out= TRUE;
                                  fprintf(output_fp,"\nERROR: Can not recognise file type >>%s<< for animation.\n",
                                                                                       tok);
                                  fprintf(output_fp,"Current options are: biosym : .arc file output\n");
                                  fprintf(output_fp,"                     txyz   : template only Xmol xyz file\n");
                                  fprintf(output_fp,"                     pxyz   : template/pore Xmol xyz file\n");
                                  fprintf(stdout,"\nERROR: Can not recognise file type >>%s<< for animation.\n",
                                                                                       tok);
                                  fprintf(stdout,"Current options are: biosym : .arc file output\n");
                                  fprintf(stdout,"                     txyz   : template only Xmol xyz file\n");
                                  fprintf(stdout,"                     pxyz   : template/pore Xmol xyz file\n");
                                  break;
                  }
		 		 if (second_token == OFF2) break;
		 		 else if (second_token != -1) bail_out= get_animation(p_animation_file) || bail_out;
		 		 break;

		    case ANALYSE: bail_out= get_analyse_filename(p_analyse_name) || bail_out; 
                         break;

		    case FORCEFIELD: bail_out= get_discover_forcefield_name() || bail_out; 
                            break;

           case SLAB: *p_is_slab=TRUE;
                break;

           case POLY: *p_is_poly=TRUE;
                break;

		    case MC_OPT: optimise_post_mc=TRUE;
		 		 break;

		    case FINAL_GULP: want_final_gulp=TRUE;
		 		 break;

		 case ANNEAL: anneal_run =TRUE;break;

		 case FINAL_CAR: want_final_car=TRUE;
        break;   

		 case MONITOR: *p_need_monitors=TRUE;
                               printf("DEBUG>> Reader found monitor directive\n");
                               if (tok = tok_get(input_fp, TRUE, FALSE))
                                 {
                                   strcpy(p_monitored->host, tok);
                                   printf("Read >>%s<< in monitor section as host atom to monitor.\n", tok);
                                 }
                               else
                                 {
                                   printf("ERROR: No atoms given with monitor directive\n");
                                   exit(0);
                                 }
                               if (tok = tok_get(input_fp, TRUE, FALSE))
                                 {
                                   strcpy(p_monitored->guest, tok);
                                   printf("Read >>%s<< in monitor section as guest atom to monitor.\n", tok);
                                 }
                               else
                                 {
                                   printf("ERROR: No guest atom given with monitor directive\n");
                                   exit(0);
                                 }
		 	       break;

		    case MC_TEMP_STEP:get_mctempstep();break;

           case MC_NO_ACTIONS:get_mc_opt_actions();break; 
case DOCK_ENERGY: bail_out=get_dock_energy() || bail_out; break;         
           case NUM_DOCK: bail_out=get_num_dock() || bail_out; break;
           case WANT_CONDOR: want_condor = TRUE; break;

        case EXXON: exxon = TRUE; break;
        case UCL_: ucl = TRUE; break;

		    case VERBOSE: 
                tok = tok_get(input_fp, FALSE, TRUE);
		 		 second_token = find_kind(tok, SECONDARY_DIRECTIVE);

                if (DEBUG) printf("Second token %s value %d\n", tok, second_token);

		 		 last_tok =tok;
                switch (second_token)
                  {
                    case ON2  : verbose = TRUE;  DEBUG= FALSE; break;
                    case OFF2 : verbose = FALSE; DEBUG= FALSE; break;
                    case DEBUG_IT  : verbose = TRUE; DEBUG=TRUE; printf("Running in DEBUG mode\n");  break;
                    case -1  : bail_out= TRUE;
                               fprintf(output_fp,"\nERROR: Verbose directive (verb) given but with incorrect option.\n");
                               fprintf(output_fp,"       Choices are: on / off / debug \n");
                               break;
                  }
                break;

/**************************************************/
/** Require user to request insight log files *****/
/**************************************************/

		    case LOGFILE: 
                tok = tok_get(input_fp, FALSE, TRUE);
		 		 second_token = find_kind(tok, SECONDARY_DIRECTIVE);

                if (DEBUG) printf("Second token %s value %d\n", tok, second_token);

		 		 last_tok =tok;
                switch (second_token)
                  {
                    case ON2  : logfile_needed = TRUE;  break;
                    case OFF2 : logfile_needed = FALSE; break;
                    case -1  : bail_out= TRUE;
                               fprintf(output_fp,"\nERROR: logfile directive (logf) given but with incorrect option.\n");
                               fprintf(output_fp,"       Choices are: on / off\n");
                               break;
                  }
                break;

		    case CHECK_INPUT: 
                tok = tok_get(input_fp, FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                last_tok =tok;
                switch (second_token)
                  {
                     case ON2  : *p_test_symm = TRUE;  break;
                     case OFF2 : *p_test_symm = FALSE; break;
                     case -1  : bail_out= TRUE;
                                fprintf(output_fp,"\nERROR: Symmetry checking requested (directive chec) ");
                                fprintf(output_fp,"but no option given.\n");
                                fprintf(output_fp,"       Choices are: on / off\n");
                                break;
                  }
		 		 break;

/*******************************************************************************/
/**** cost_function directive to control which cost function to use ************/
/**** added 7th March 1997 Dave Willock                             ************/
/*******************************************************************************/

           case COST_FUNCTION:
                tok = tok_get(input_fp,FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                last_tok =tok;

                switch (second_token)
                  {
                     case  STERIC2   : *p_which_test= STERIC_COST; 
                                       steric = TRUE;
                                       non_bonded = FALSE;
                                       charges = FALSE;  
                                       break;

                     case  NON_BOND2 : *p_which_test= NON_BOND_COST; 
                                       steric = FALSE;
                                       non_bonded = TRUE;
                                       charges = FALSE;  
                                       break;

                     case  ENERGY2   : *p_which_test= ENERGY_COST; 
                                       steric = FALSE;
                                       non_bonded = TRUE;
                                       charges = TRUE;  
                                       break;

                     case -1  : bail_out= TRUE;
                                fprintf(output_fp,"\nERROR: Cost function type defined (directive cost)");
                                fprintf(output_fp,"but with no valid option set.\n");
                                fprintf(output_fp,"       Choices are: ster / non_ / ener\n");
                                break;
                  }
                break;

/*******************************************************************************/
/**** just_position directive to allow just moving molecule around  ************/
/**** added 2nd April 1997 Dave Willock                             ************/
/*******************************************************************************/

           case JUST_POSITION:
                tok = tok_get(input_fp,FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                last_tok =tok;

                switch (second_token)
                  {
                     case  ON2       : *p_just_position = TRUE;  break;
                     case  OFF2      : *p_just_position = FALSE; break;
                     case -1  : bail_out= TRUE;
                                fprintf(output_fp,"\nERROR: just_position directive given with invalid option");
                                fprintf(output_fp,"but with no valid option set.\n");
                                fprintf(output_fp,"       Choices are: on / off\n");
                                break;
                  }
                break;


		    case MAX_TEMPLATES: bail_out= get_max_templates() || bail_out; break;

           case PROB_TEST: bail_out= get_prob_test() || bail_out; break;

           case TEMPERATURE: bail_out= get_temperature() || bail_out; break;

		    case STOP_CUTOFF: bail_out= get_stop_cutoff() || bail_out; break;

		    case VDW_SCALE: bail_out= get_vdw_scale() || bail_out; break;

 		    case MINIMIZER: tok = tok_get(input_fp,FALSE, TRUE);
                           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                           last_tok =tok;
                           switch (second_token)
                             {
                                case CONSTRAIN : tok = tok_get(input_fp,FALSE, TRUE);
                                                 third_token = find_kind(tok, TERTIARY_DIRECTIVE);
                                                 last_tok =tok;
                                            
                                                 switch (third_token)
                                                    {
                                                       case ATOM : get_atoms_to_fix(); break;
                                                    }
                                                 break;
 
                                default : bail_out= get_minimizer_name(tok) || bail_out; break;
                             }
                           break;
                 

		    case INITIAL_MINIMIZE:  tok = tok_get(input_fp,FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                           printf("Reader processing initial minimiser flag\n");
                last_tok =tok;
                switch (second_token)
                  {
                    case TEMPLATE: tok = tok_get(input_fp,FALSE, TRUE);
                                   third_token = find_kind(tok, TERTIARY_DIRECTIVE);
                                   last_tok =tok;
                           printf("..... This is for the template.\n");

                                   switch (third_token)
                                     {
                                       case ON3 : initial_minimize_template = TRUE; 
                                                  printf("In read_input initial_minimize flag set TRUE\n"); break;
                                       case OFF3: initial_minimize_template =FALSE; 
                                                  printf("In read_input initial_minimize flag set FALSE\n"); break;
                                       case YES3 : initial_minimize_template = TRUE; 
                                                   printf("In read_input initial_minimize  yes flag set TRUE\n"); break;
                                       case NO3:   initial_minimize_template =FALSE; 
                                                   printf("In read_input initial_minimize no flag set FALSE\n"); break;
                                       case -1 :  bail_out= TRUE;
                                                  fprintf(output_fp, "\nERROR: Initial minimisation of template directive");
                                                  fprintf(output_fp, " (init temp) given but without an option.\n");
                                                  fprintf(output_fp,"       Choices are: on / off / yes/ no \n");
                                                  break;

                                     }
                                   break;
		 
                case INPORE:  tok = tok_get(input_fp,FALSE, TRUE);
                    third_token = find_kind(tok, TERTIARY_DIRECTIVE);
                    last_tok =tok;
                    switch (third_token)
                      {
                        case ON3 : initial_minimize_inpore = TRUE; break;
                        case OFF3: initial_minimize_inpore =FALSE; break;
                        case YES3 : initial_minimize_inpore = TRUE;  break;
                        case NO3:   initial_minimize_inpore =FALSE;  break;
                        case -1 : bail_out= TRUE;
                                  fprintf(output_fp, "\nERROR: Initial minimisation of template in the pore directive");
                                  fprintf(output_fp, " (init in_p) given but without an option.\n");
                                  fprintf(output_fp,"       Choices are: on / off\n");
                                  break;
                      }
                   break;
                }
                break;

		    case STRATEGY_FILE: tok = tok_get(input_fp,FALSE, TRUE);
                               second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                               last_tok =tok;
                               switch (second_token)
                                 {
		 		 		 		    case TEMPLATE: get_template_strategy_name();break;
		 		 		            case INPORE:   get_inpore_strategy_name();break;
                                 }
		 		 		 		 break;

		    case DISCOVER: case C2DISCOVER: tok = tok_get(input_fp,FALSE, TRUE);
                                           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                                           last_tok =tok;
                                           switch (second_token)
                                             {
                                                case PATH:     get_discover_path();
                                                break;

                                                case FORCEFIELD: bail_out =get_discover_forcefield_name() || bail_out; 
                                                break;
                                             }		 
                                           break;

		    case MOPAC: tok = tok_get(input_fp,FALSE, TRUE);
                       second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                       last_tok =tok;
                       switch (second_token)
                         {
                            case TEMPLATE: get_mopac_cmds(1);break;
                            case INPORE:   get_mopac_cmds(2);break;
                            case NAME:     get_mopac_name();break;
                            case PATH:     get_mopac_path();break;
                         }
                       break;

           case GULP: tok = tok_get(input_fp,FALSE, TRUE);
                       second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                       last_tok =tok;
                       switch (second_token)
                         {
                            case TEMPLATE: get_gulp_cmds(1);break;
                            case INPORE:   get_gulp_cmds(2);break;
                            case NAME:     get_gulp_name();break;
                            case PATH:     get_gulp_path();break;
                         }
                       break;

           case GULP_CMD_LINE: get_gulp_cmds(1); get_gulp_cmds(2);break;

           case GULP_EXTRA: get_gulp_extra();break;


		    case LIBRARY: tok = tok_get(input_fp,FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                last_tok =tok;
                switch (second_token)
                {
                case FORCEFIELD: get_forcefield_library();break;
                case TEMPLATE: get_template_library();break;
                case FRAGMENT: get_fragment_library(frag_lib, just_count);break;
		 		 		 		 }
		 		 		 		 break;

       case DEFAULTS_FILE: get_defaults_file();break;

       case SEED_TEMPLATE: tok = tok_get(input_fp,FALSE, TRUE);
                           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                           last_tok =tok;
                           switch (second_token)
                             {

/***********************************************************************/
/***** only a single seed molecule *************************************/
/***** as previously               *************************************/
/***********************************************************************/
 case OFF2: bail_out = get_seed_file(p_seed_type) || bail_out;
                                              break;
                               case MOLECULE: bail_out = get_seed_file(p_seed_type) || bail_out;
                                              break;

                               case FRAGMENT: bail_out = get_seed_fragment(p_seed_type, 
                                                                           p_num_of_seeds,
                                                                           p_seeds_to_use, 
                                                                           p_random_seed) || bail_out;
                                              break;

/***********************************************************************/
/***** NEW for reading in multiple seeds *******************************/
/***********************************************************************/

                               case ARCHIVE: bail_out = get_seed_archive(p_seed_type,
                                                                         p_num_of_seeds,
                                                                         p_seeds_to_use, 
                      		 		 		 		                          p_read_all_frames_as_seeds,
                                                                         p_number_of_frames_in_arc)
                                                                         || bail_out;
                               break;
                             }
                           break;

       case PEEK_FREQ: get_peek_freq(p_peek_freq);break;

       case RANDOM:    
             tok = tok_get(input_fp,FALSE, TRUE);
 		      second_token = find_kind(tok, SECONDARY_DIRECTIVE);
             last_tok =tok;
             switch (second_token)
		 		 {
                  case ON2  :  *p_randomisation_method =  0; break;
                  case OFF2 :  *p_randomisation_method = -1; break;
		 		 }
             break;

       case SYMMETRY: if (bail_out) printf("bail out set ass symmetry called\n");
                      bail_out=get_symmetry() || bail_out ;
                      if (bail_out) printf("bailout set in get symmetry\n");
                      break;

       case CONCENTRATION_LIMIT:  get_concentration_limits();break;

       case FORBIDDEN_BOND: get_forbidden_bonds();break;

       case ALLOWED_TORSIONS: get_allowed_torsions(); break;

       case WEIGHTS: tok = tok_get(input_fp, FALSE, TRUE);
                second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                printf("Weights secondary token is %s %d\n", tok, second_token);
                last_tok =tok;
                switch (second_token)
                {
		 		 		 		 case FRAGMENT: get_frag_weight(frag_weights, just_count);break;
		 		 		 		 case ACTION: get_action_weight();break;
		 		 		 		 case HYDROGEN: 
		 		 		 		 		 		 get_hydrogen_weight(p_initial_hydrogen_weight,
                         		 		 		 		 		 p_increment_hydrogen_weight);
		 		 		 		 		 		 break;
                                case SEED: get_seed_weight(seed_weights, just_count); break;

                }
		 		 break;

       case BOX: tok = tok_get(input_fp,FALSE, TRUE);
           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
           last_tok =tok;
           switch (second_token)
               {
                 case LIMITS:   bail_out = get_box_limits() || bail_out ;break;
                 case FRACTION: bail_out = get_box_fraction() || bail_out ;break;
               }
               break;

       case MODIFY: tok = tok_get(input_fp,FALSE, TRUE);
           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
           last_tok =tok;
           switch (second_token)
             {
                case ATTEMPTS: get_modify_attempts();break;
             }
             break;

       case ROCK: tok = tok_get(input_fp,FALSE, TRUE);
           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
           last_tok =tok;
           switch (second_token)
             {
               case STEP: get_rock_step();break;
               case FIX: get_rock_fix_atom();break;
               case ATTEMPTS: get_rock_attempts();break;
             }
             break;

       case SHAKE: tok = tok_get(input_fp,FALSE, TRUE);
           second_token = find_kind(tok, SECONDARY_DIRECTIVE);
           last_tok =tok;
           switch (second_token)
             {
                case STEP: get_shake_step();break;
                case ATTEMPTS: get_shake_attempts();break;
             }
             break;
       case FORCE_DIHEDRALS: get_forced_dihedrals();break;

       case LINE_RESTRAINT: tok = tok_get(input_fp,FALSE, TRUE);
                            second_token = find_kind(tok, SECONDARY_DIRECTIVE);
                            last_tok =tok;
                            line_restraints = TRUE;

                            switch (second_token)
                              {
                                case ELEMENTS: tok = tok_get(input_fp,FALSE, TRUE);
                                               third_token = find_kind(tok, TERTIARY_DIRECTIVE);
                                               last_tok =tok;
                                               switch (third_token)
                                                  {
                                                    case HOLD :  line_hold=TRUE;
                                                                 get_restraint_elements(); break;
                                                  }
                                               break;
                              }
                            break;
          case TETHER_RESTRAINT: bail_out= get_tether_labels() || bail_out; break;
       }
  if (DEBUG && bail_out) printf("decided to bail out before here\n");

  }
  
 		 if ((steric == FALSE) && (non_bonded == FALSE) && (charges == FALSE)  )
		  { 
 		   fprintf(output_fp, "No template - pore interactions specified!\n"); 
 		   fprintf(output_fp, "Either steric, non_bonded or charges must be \n"); 
 		   fprintf(output_fp, "set to ON\n\n"); 
 		   fprintf(output_fp, "Aborting.....\n"); 
		   exit(0);
		  }

/*******************************************************************************/
/* errors will have been reported and the code should wind up if bail_out TRUE */  
/*******************************************************************************/

  return (bail_out); 
}

/******************************************************************************
  tok_get.c : parses command lines in a semi intelligent fashion
IMportant change 22/6/95 DWL Now ignores EVERYTHING after a #
******************************************************************************/


char * tok_get(FILE *input_fp, int skip_lines, int to_lower_case)
{
  int          i;
  char         *tok;
  static char  target[BUFFER];
  static char  buf[BUFFER];
  static char  *line, *last_tok;
  extern int   read_new_line, line_no;
  char         *comment;		 

  i = 0;

/**************************************************************************/
/*** If the last token read was a NULL character last_tok will act ********/
/*** like a FALSE logical and force us to read the next line **************/
/**************************************************************************/

  if (!last_tok && skip_lines)
    {

/**************************************************************************/
/**** Read in the next line from input_fp, if it is NULL return NULL ******/
/**************************************************************************/

      if (!(line = fgets(target, BUFFER, input_fp))) return(line);

/**************************************************************************/
/**** Remove commented lines or parts of lines that are comments **********/
/**************************************************************************/

      if ((comment = strchr(line,'#')) != NULL) 
        {

/**************************************************************************/
/****** chop off the bits after the hash **********************************/
/**************************************************************************/

           *(comment) = '\0';

/**************************************************************************/
/****** if this is a simple comment line it will now have length 0 ********/
/****** so simply return a single # character                      ********/
/**************************************************************************/

           if (strlen(line) < 1) 
             {
               *line= '#';
               *(line+1)= '\0';
               return(line);
             }
         }

/**************************************************************************/
/***** read_new_line says we have just read a new line in !! **************/
/**************************************************************************/

       read_new_line = TRUE;
       line_no++;
    }
  else if (!last_tok && !skip_lines)
    {
       return(NULL);
    }

/**************************************************************************/
/****** If we have just read a new line copy it to the buffer *************/
/**************************************************************************/

  if (read_new_line) strncpy(buf, line, BUFFER);

/**************************************************************************/
/***** Now process the latest string of information ***********************/
/**************************************************************************/

  for(;;)
    {

/**************************************************************************/
/***** tok is read from the last line if some tokens are left *************/
/**************************************************************************/

      tok = (read_new_line) ? strtok(buf, " :;,\n\r\t") : strtok(NULL, " :;,\n\r\t");

      read_new_line = (!tok) ? TRUE : FALSE; 
      last_tok = tok;
      if (!tok || !isalpha(*tok)) break;

      if (to_lower_case) for (i = 0; tok[i] != '\0'; i++) tok[i] = tolower(tok[i]);
      if (find_kind(tok, REMAINING_DIRECTIVE) != UNIT) break;
      continue;
    }
  return(tok);
}
  
/******************************************************************************
  Find_kind takes a parsed token, and matches it against the directive list
  contained in reader.h, returns a BLANK_DIRECT value if token is unrecognised,
  and a UNIT value if the token is extraneous clutter..
******************************************************************************/

int find_kind(char *token, int level)
{

list null_namelist[]   = {NULL_DIRECTIVE_LIST};
list first_namelist[]  = {FIRST_DIRECTIVE_LIST};
list second_namelist[] = {SECOND_DIRECTIVE_LIST};
list third_namelist[] = {THIRD_DIRECTIVE_LIST}; 

/*****************************************************************************/
/*** find_kind altered to be context aware by knowing the level from *********/
/*** which the token was read. Dave Willock March 1997 ***********************/
/*****************************************************************************/

  int i, j,k,  l;
  
  i = 0;
  j = 0;
  k = 0;
  l = 0;

  if (!token) return(-1);

if (level == PRIME_DIRECTIVE || level == REMAINING_DIRECTIVE)
  {
    do
      {
        if (!strncmp(token, first_namelist[j].directive,4))
          {
		      return(first_namelist[j].token_index);
		   }
      }
  while(first_namelist[j++].token_index != BLANK_DIRECT);
  }
else if (level == SECONDARY_DIRECTIVE)
  {
    do
      {
        if (!strncmp(token,second_namelist[k].directive,4))
  		   {
		     return(second_namelist[k].token_index);
		   }
      }
    while(second_namelist[k++].token_index != BLANK_DIRECT); 
  }
else if (level == TERTIARY_DIRECTIVE)
  {
    do 
     { 
/**********************************************************/
/**** For 3rd directive compare only to 2nd character *****/
/**********************************************************/
       if (!strncmp(token,third_namelist[l].directive,2))
		  { 
		   return(third_namelist[l].token_index); 
		  } 
     } 
    while(third_namelist[l++].token_index != BLANK_DIRECT);
  }

    do
      {
        if (!strncmp(token, null_namelist[i].directive,2))
		   {
		     return(UNIT);
		   }
      }
    while(null_namelist[i++].token_index != BLANK_DIRECT);

  return(BLANK_DIRECT);
}
  
/******************************************************************************
  Get token routines for individual datatypes, token is taken from line and
  cast to apprpriate type...
******************************************************************************/
 
int get_integer(int *p_error, int skip_lines)
{
#include "header.h"
  int answer;
  char *tok;

  tok= tok_get(input_fp, skip_lines, TRUE);

/*******************************************************************/
/**** Check all is well, tok will give FALSE if it is NULL *********/
/**** If there is a string there check the first character *********/
/**** to see if it is a valid decimal number               *********/
/**** Dave Willock March 1997 **************************************/
/*******************************************************************/

  if (tok)
    {
       if ((*tok >= '0' && *tok <= '9') || *tok == '-')
          {
             answer = atoi(tok);
             *p_error = FALSE;
          }
       else
          {
             answer= 0;
             *p_error = TRUE;
          }
    }
  else
    {
       answer= 0;
       *p_error = TRUE;
    }
  return (answer);

}

char * get_string(void)
{
#include "header.h"
  char *dummy;

  dummy = tok_get(input_fp,FALSE, FALSE);
  return (dummy);
}

double get_double(int *p_error, int skip_lines)
{
#include "header.h"
double answer;
char *tok;
  
  tok= tok_get(input_fp, skip_lines, TRUE);

/*******************************************************************/
/**** Check all is well, tok will give FALSE if it is NULL *********/
/**** If there is a string there check the first character *********/
/**** to see if it is a valid decimal number               *********/
/**** Dave Willock March 1997 **************************************/
/*******************************************************************/

  if (tok)
    {
       if ((*tok >= '0' && *tok <= '9') || *tok == '-' || *tok == '.')
          {
             answer = atof(tok);
             *p_error = FALSE;
          }
       else
          {
             answer= 0;
             *p_error = TRUE;
          }
    }
  else
    {
       answer= 0;
       *p_error = TRUE;
    }
  return (answer);
}
  
/******************************************************************************
  Contains all the files used in parse_line to read specific data into the 
  correct variables.  Also contains error function returning error parsing
  file message, including error line and lasdt token read
******************************************************************************/

void parse_error(char *last_tok)
{
#include "header.h"
  if (!(output_fp =fopen(outputfile, "w")))
    {
      fprintf(output_fp,"\nERROR: opening file \"%s\"\n", outputfile);
      exit(1);
    }
  fprintf(output_fp,"\n\terror parsing input file: %s", outputfile);
  fprintf(output_fp,"\n\t\tlast token read: %s", last_tok);
  fprintf(output_fp,"\t at or near line: %d\n\n", line_no);
  fclose(output_fp);
  exit(1);
} 

void get_title(void)
{ 
#include "header.h"
  fgets(title,BUFFER,input_fp);
  return;
}

int get_pore_file(atom *pore, int just_count)
{
#include "header.h"
  int   iloop,ineigh, i;
  char *p_key;
  int error;

  error=FALSE;

  strcpy(pore_file, get_string());

  if (!(pore_fp = fopen(pore_file, "r")))
		 {
                fprintf(output_fp,"\nERROR: Problem opening Pore file %s\n", pore_file);
                fflush(output_fp);
                printf("\nERROR: Problem opening Pore file %s\n", pore_file);
                error= TRUE;
		 }
  else
		 {

/*****************************************************/
/****** READ IN THE PORE COORDINATES *****************/
/*****************************************************/

		 fgets(dummy_head,BUFFER,pore_fp);
		 fgets(dummy_head,BUFFER,pore_fp);

/*****************************************************/
/****** check if this is a periodic cell pore ********/
/*****************************************************/

        p_key = "PBC=ON";
        pbc = locate_string( p_key, &dummy_head[0], BUFFER );

        if (pbc)
          {
             printf("This is a periodic pore assuming XYZ format!!!\n");
          }
        else
          {
             printf("This is not a periodic pore\n");
          }

		 fgets(pore_title,BUFFER,pore_fp); /* hang on to the title line */
		 fgets(dummy_head,BUFFER,pore_fp);

/*****************************************************************/
/****** if this is a periodic file pick up the cell vectors ******/
/*****************************************************************/

       if (pbc)
         {
            fgets(buffer,BUFFER,pore_fp);

            sscanf(buffer, "%*s%lf%lf%lf%lf%lf%lf", &abc[0],
                     &abc[1], &abc[2], &abc[3], &abc[4],
                     &abc[5]);

            printf("Read Lattce parameters as: ");
            for (iloop = 0; iloop < 6; iloop++)
              {
                printf("%10.6f   ",abc[iloop]);
              }
            printf("\n");
         }

    i=0;
    while(fgets(buffer,BUFFER,pore_fp) != NULL)
		   {
      if (strstr(buffer,"end") != NULL)
       { /*end of car */
       if (i==0) 
         {
            printf("Pore is an empty file\n");
            return(-10);
         }
       break;
       }
      else
       {
     		 if(just_count==TRUE)
        {
          i++;
        }   
        else
        {
          sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf", &(pore[i].label[0]),
                                                  &pore[i].x,
                                                  &pore[i].y, 
                                                  &pore[i].z,
                                                  &(pore[i].group[0]),
                                                  &(pore[i].group_no[0]),
                                                  &(pore[i].pot[0]),
                                                  &(pore[i].elem[0]),
                                                  &pore[i].part_chge);
        

/* set starting values of neighbour indexing bits of the atom structures */

          pore[i].num_neigh=-1;
          for (ineigh=0; ineigh < 4; ineigh++)
                                      pore[i].neighb[ineigh]=-1;
          i++;
        }
       }
      num_pore_atoms = i-1;
      }
     fclose(pore_fp);
    }

  return(error);
}

int get_gulp_pots_file(void)
{
 #include "header.h"
 int error;

 error = FALSE;

 strcpy(gulp_pots_file, get_string());

 if (!(gulp_pots_fp = fopen(gulp_pots_file, "r")))
        {
                fprintf(output_fp,"\nERROR: Problem opening gulp potentials file %s\n",                                                 gulp_pots_file);
                error= TRUE;
        }


 fclose(gulp_pots_fp);
  return(error);
}

void get_gulp_name(void)
{
#include "header.h"
        strcpy(gulp_root, get_string());
        return;
}

void get_gulp_path(void)
{
#include "header.h"
         strcpy(gulp_path, get_string());
         return;
}

void get_gulp_extra(void)
{
#include "header.h"
char extra_file[100];

/***read in next line for commands***/
fgets(buffer,BUFFER,input_fp);
if (buffer == NULL)
        {
        fprintf(output_fp,
        "Error Reading Input File. Expecting GULP commands. End of File Reached\n");
        exit(-1);
        return;
        }

strcpy(extra_file, "extra_gulp.lines");
extra_gulp_fp = fopen(extra_file, "w");

while (strstr(buffer, "end extralines for gulp") == 0)
        {
         fprintf(extra_gulp_fp, "%s", buffer);
         fgets(buffer,BUFFER, input_fp);
        }
fclose(extra_gulp_fp);
return;
}

void get_gulp_cmds(int whichone)
{
#include "header.h"
/***read in next line for commands***/
fgets(buffer,BUFFER,input_fp);
if (buffer == NULL)
        {
        fprintf(output_fp,
        "Error Reading Input File. Expecting GULP commands. End of File Reached\n");
        exit(-1);
        return;
        }

if (whichone ==1)
        {
        strcpy(gulp_cmdline_molecule, buffer);
        }
else if (whichone ==2)
        {
        strcpy(gulp_cmdline_inpore, buffer);
        }
else
        {
        fprintf(output_fp,"SERIOUS ERROR in reader! Contact authors 'cos it has gone haywire.\n");        }
return;
}

int get_dock_energy(void)
{
#include "header.h"
int error;
dock_energy = get_double(&error, FALSE);
return(error);
}

int get_num_dock(void)
{
#include "header.h"
int error;
num_dock = get_integer(&error, FALSE);
return(error);
}

int get_vdw_scale(void)
{
#include "header.h"
int error;

  vdw_scale = get_double(&error, FALSE);

  if (error)
    {
      fprintf(output_fp,"\nERROR: Non-numeric supplied for Van der Waals scaling parameter (directive vdw_).\n");
    }
  else
    {
      if (vdw_scale <= 0)
        {
           fprintf(output_fp,"\nERROR: Van der Waals scaling parameter (directive vdw_) is zero or negative.\n"); 
        }
    }
  return (error);
}

int get_prob_test(void)
{
#include "header.h"
int error;

  prob_test= get_double(&error, FALSE);

  if (error)
    {
      fprintf(output_fp,"\nERROR: Non-numeric supplied for probability of testing  parameter (directive prob).\n");
    }
  else
    {

/***********************************************************************************/
/*** Alteration 25th Sept. 96 DJW: Make input file give prob_test as percentage! ***/
/***********************************************************************************/

       prob_test= prob_test/100.0;
       if (prob_test > 1.0)
         {
            fprintf(output_fp,"WARNING: Test Probability was >100 %% (%9.4f).\n",
                                                                        100.0*prob_test);
            fprintf(output_fp,"       : Setting to 100%% and continuing\n");
		     prob_test = 1.0;
		  }		 
       if (prob_test <= 0.0)
         {
		    fprintf(output_fp,"\nERROR: Test Probability was <=0.0 %% (%9.4f).\n",
                                                                        100.0*prob_test);
           error=TRUE;
		  }
    }
  return (error);
}

int get_stop_cutoff(void)
{
#include "header.h"
int error;

  stop_ctf = get_double(&error, FALSE);
  
  return(error);
}

int get_temperature(void)
{
#include "header.h"
int error;

  temperature = get_double(&error, FALSE);
  
  return(error);
}

int get_max_templates(void)
{
#include "header.h"
int error;

  max_templates = get_integer(&error, FALSE);
  if (error) fprintf(output_fp,"\nERROR: The maximum number of templates directive (maxt) is given but it has no value.\n");
  
  if (max_templates < 0)
     {
        fprintf(output_fp,"\nERROR: The maximum number of templates cannot be negative.\n");
        error=TRUE;
     }
  return (error);
}
int get_nb_cutoff(void)
{
#include "header.h"
int error;

  nb_ctf = get_double(&error, FALSE);
  
  if (error)
    {
       nb_ctf = NB_CUTOFF_DEFAULT;
       fprintf(output_fp,"\nERROR Non-bond cutoff directive (nb_cu) supplied with no value\n");
    } 

  nb_ctf_2 = nb_ctf * nb_ctf;
  return(error);
}

int get_ring_cutoff(void)
{
#include "header.h"
int error;
 
  ring_ctf = get_double(&error, FALSE);

  if (error)
    {
      fprintf(output_fp,"\nERROR: Non-numeric supplied for ring formation cut off parameter (directive ring).\n");
    }
  else
    {
      ring_ctf_2 = ring_ctf * ring_ctf;
    }

  return(error);
}

void get_modify_attempts(void)
{
#include "header.h"
int error;

  num_modify_attempts = get_integer(&error, FALSE);
  return;
}

void get_rock_attempts(void)
{
#include "header.h"
int error;

  num_rock_attempts = get_integer(&error, FALSE);
  return;
}

void get_rock_step(void)
{
#include "header.h"
int error;

  max_rock_step = get_double(&error, FALSE);
  return;
}

void get_rock_fix_atom(void)
{
#include "header.h"

printf("Fixed option not yet available\n");
exit(0);
}

void get_shake_attempts(void)
{
#include "header.h"
int error;

  num_shake_attempts = get_integer(&error, FALSE);
  return;
}

void get_shake_step(void)
{
#include "header.h"
int error;

  max_shake_step = get_double(&error, FALSE);
  return;
}

int get_ch_cutoff(void)
{
#include "header.h"
int error;

  ch_ctf = get_double(&error, FALSE);
  if (error)
     {
        fprintf(output_fp,"\nERROR: charge cutoff directive (ch_cu) supplied with no value\n");
     }
  return(error);
}

void get_peek_freq(int *p_peek_freq)
{
#include "header.h"
int error;

  *p_peek_freq = get_integer(&error, FALSE);
  return;
}

void get_mctempstep(void)
{
#include "header.h"
int error;
                                                                                                                  
  mc_temp_step = get_double(&error, FALSE);
  printf("KJ, mc_temp_step = %10.6f\n", mc_temp_step);
  return;
}

void get_mc_opt_actions(void)
{
#include "header.h"
int error;
                                                                                                                  
  mc_modi_step = get_integer(&error, FALSE);
  printf("KJ, mc_modi_step = %d\n", mc_modi_step);
  return;
}



int  get_box_fraction(void)
{
#include "header.h"
int error;

  box_fraction = get_double(&error, FALSE);
  printf("DD>> box_fraction = %f\n", box_fraction);
  if (error) return TRUE;
  return FALSE;
}

int get_box_limits(void)
{
#include "header.h"
int read_line(FILE *fp, int *p_ichar);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

int next_space( int *p_ichar, int start, int num_of_chars );

int itsanum,place, ichar[BUFFER];
int limit, num_of_chars, failed;

num_of_chars= read_line(input_fp, &ichar[0]);

itsanum= 1;
place = 0;
limit= -1;
failed = FALSE;

while (itsanum)
   {
    limit++;

		 if (limit > 6)
		 		 {
		 		 fprintf(output_fp, "\nERROR: More than 6 box limits suppplied!!\n\n");
		 		 exit(1);
		 		 }

    box_limits[limit]= get_doub(&ichar[0], num_of_chars, &place, &itsanum);
		 printf("limit %d %f\n",limit,box_limits[limit]);

    place= next_space( &ichar[0], place, num_of_chars ); 
   }


if (limit !=6)
		 		 {
    		 fprintf(output_fp, "\nERROR: Not enough box limits suppplied!!\n");
    		 fprintf(output_fp, "       Required as -x,+x,-y,+y,-z,+z\n\n");
        failed= TRUE;
        }

/* Spot stupid limits Dave March 99 */

if ( box_limits[0] > box_limits[1] )
  {
    fprintf(output_fp, "\nERROR: xmin (%10.6f) is greater than xmax (%10.6f) in box definition\n", 
                                    box_limits[0], box_limits[1]);
    failed= TRUE;
  }
if ( box_limits[2] > box_limits[3] )
  {
    fprintf(output_fp, "\nERROR: ymin (%10.6f) is greater than ymax (%10.6f) in box definition\n", 
                                    box_limits[2], box_limits[3]);
    failed= TRUE;
  }
if ( box_limits[4] > box_limits[5] )
  {
    fprintf(output_fp, "\nERROR: zmin (%10.6f) is greater than zmax (%10.6f) in box definition\n", 
                                    box_limits[4], box_limits[5]);
    failed= TRUE;
  }

user_box = TRUE;
return failed;
} 

void get_forbidden_bonds(void)
{
#include "header.h"
		 int i, is_element, error;
		 
		 num_forbidden_bonds = get_integer(&error, FALSE);
		 if (num_forbidden_bonds <1)
		 		 {
		 		 fprintf(output_fp,"Error, number of forbidden bonds not valid\n");
		 		 fprintf(output_fp,"Exiting....\n");
		 		 exit(-1);
		 		 }
		 have_forbidden_bonds = TRUE;

		 for (i=0;i<num_forbidden_bonds;i++)
        {
        fgets(buffer,BUFFER,input_fp);
		 		 printf("i=%i buffer = %s\n", i, buffer);
        sscanf(buffer, "%s %s", forbidden_bond[i].atom1,
                                forbidden_bond[i].atom2);

		 		 /***check valid atoms types***/
		 		 is_element = check_element(forbidden_bond[i].atom1);
		 		 if (is_element == 0)
		 		 		 {
		 		 		 fprintf(output_fp,
		 		 		 		 "Error, Invalid Element %s in Forbidden Bonds Input\n",
		 		 		 		 		 		 		 		 forbidden_bond[i].atom1);
        		 fprintf(output_fp,"Exiting....\n");
		 		 		 exit(-1);

		 		 		 }		 

		 		 is_element = check_element(forbidden_bond[i].atom2);
		 		 if (is_element == 0)
		 		 		 {
		 		 		 fprintf(output_fp,
		 		 		 		 "Error, Invalid Element %s in Forbidden Bonds Input\n",
		 		 		 		 		 		 		 		 forbidden_bond[i].atom2);
        		 fprintf(output_fp,"Exiting....\n");
		 		 		 exit(-1);
		 		 		 }
		 		 
		 		 }

		 return;
}
void get_concentration_limits(void) 
{
#include "header.h"
		 int i,j,error;
		 int found_elem;

		 num_conc_limits = get_integer(&error, FALSE);
		 printf("Getting %i concentration limits\n",num_conc_limits);

		 have_conc_limits = TRUE;   /* flag for testing on concentration */

		 if (num_conc_limits < 1) 
		 		 {
		 		 fprintf(output_fp,"Error, number of concentration limits not valid\n");
		 		 fprintf(output_fp,"Exiting....\n");
		 		 exit(-1);
		 		 }
		 
		 for (i=0;i<num_conc_limits;i++)
		 		 {
        fgets(buffer,BUFFER,input_fp);
        sscanf(buffer, "%s %i", &(atom_limit[i].atom_type[0]), &atom_limit[i].num);

		     /***check valid type and number***/

		 		 if (atom_limit[i].num <0)
		 		 		 {
		 		 		 fprintf(output_fp,"Error, concentration limit %i for atom %s not valid\n", 
                                                               atom_limit[i].num, atom_limit[i].atom_type);
		 		 		 fprintf(output_fp,"Exiting....\n");
		 		 		 exit(-1);
		 		 		 }		 
		 		 
		 		 j=0;
		 		 found_elem = FALSE;
		 		 while ((found_elem == FALSE) && (j<NUM_ELEMENTS))
		 		 		 {
		 		 		 if (strcmp(atom_limit[i].atom_type, period_table[j].elem) == 0)
		 		 		 		 {
		 		 		 		 found_elem = TRUE;
		 		 		 		 }
		 		 		 j++;
		 		 		 }

		 		 if (found_elem == FALSE)
		 		 		 {		 
		 		 		 fprintf(output_fp,"\nERROR: Element %s not recognised in concentration limit input\n", 
                                                                                        atom_limit[i].atom_type);
		 		 		 fprintf(output_fp,"Exiting....\n");
		 		 		 exit(-1);
		 		 		 }
		 		 }
		 
return;
}
int get_symmetry(void)
{
#include "header.h"
  int i,j,ifrac,found_fraction;
  int error, blank_line;
  char translation[20];

  num_symm_ops = get_integer(&error, FALSE);
  printf("Number of symmetry operators: %i\n", num_symm_ops);		 

  if (num_symm_ops > MAX_SYMMOPS)
    {
      fprintf(output_fp,"ERROR: Number of symmetry operators requested, %d, exceeds ",
                                                   num_symm_ops);
      fprintf(output_fp,"maximum allowed, %d.\n", MAX_SYMMOPS);
      printf("ERROR: Number of symmetry operators requested, %d, exceeds ",
                                                   num_symm_ops);
      printf("maximum allowed, %d.\n", MAX_SYMMOPS);
      return TRUE;
    }
  num_symm_ops--;


  if (num_symm_ops <0)
		 {
    		 fprintf(output_fp, "\nERROR: Number of Symmetry Operators < 1\n");
		 		 exit(1);
		 }		 
  symm_set = TRUE;

  /* now read in the 3*n lines of operators */
  /* LITTLE ERROR CHECKING! */

  for(i=0;i<=num_symm_ops; i++) 		 /* for each operator matrix */
		 {
		 for (j=0;j<3;j++) 		 		 		 /* for each line in the matrix */
		 		 {
        blank_line= TRUE;
        while (blank_line)
          {
            fgets(buffer,BUFFER,input_fp);
            blank_line = strlen(buffer) <  5 || buffer[0] == '#';
            if (blank_line) printf("blank line\n"); 
          }

        sscanf(buffer, "%lf%lf%lf%s", 
		 		 		 		 &symm[i].matrix[(j*3)],
		 		 		 		 &symm[i].matrix[(j*3)+1],
		 		 		 		 &symm[i].matrix[(j*3)+2],
		 		 		 		 translation);

		 		 /* scan list of allowed fractions for match */
		 		 found_fraction = FALSE;
        ifrac = 0;
        while (ifrac < num_fractions && found_fraction == FALSE)
		 		 		 {

		 		 		 if (strcmp(translation,fraction_list[ifrac].label) == 0)
		 		 		 		 {
		 		 		 		 symm[i].translation[j] = fraction_list[ifrac].value;
		 		 		 		 found_fraction = 1;
		 		 		 		 }		 
		 		 		 ifrac++;
		 		 		 }
		 		 if (found_fraction == FALSE)
		 		 		 {
		 		 		 if (strstr(translation,"/") != FALSE)
		 		     		 {
		 		 		 		 fprintf(output_fp, "\nERROR: Unknown fraction %s in Symmetry Operator\n",
		 		 		 		 		 		 translation);
		 		 		 		 fprintf(output_fp,"     : Use decimal fraction\n");
		 		 		 		 exit(1);
		 		 		 		 }

		 		 		 else
		 		 		 		 {
		 		 		 		 symm[i].translation[j] = atof(translation);
		 		 		 		 }
		 		 		 }
		 		 printf("symmetry line: %f %f %f %f\n", 
		 		 		 		 symm[i].matrix[(j*3)],
		 		 		 		 symm[i].matrix[(j*3)+1],
		 		 		 		 symm[i].matrix[(j*3)+2],
		 		 		 		 symm[i].translation[j]);
		 		 }
		 }
return FALSE;

}
int get_tether_labels(void)
{
#include "header.h"

int error, num_args;

 have_tethers= TRUE; 
 error= FALSE;

 fgets(buffer,BUFFER,input_fp);

  num_args = sscanf(buffer, "%s%s%lf%lf", &(tether.A[0]),
                                          &(tether.B[0]),
                                          &(tether.r0),
                                          &(tether.k)     );

  if (num_args < 4)
    {
      fprintf(output_fp, "\nERROR>> Too few arguements for tethering directive given\n");
      fprintf(output_fp, "ERROR>> causing incorrect format for tethering directive on input line:\n");
      fprintf(output_fp, "\n%s\n", buffer );
      fprintf(output_fp, "ERROR>> Expecting format : <host atom label> <seed atom label> <r0> <k>\n\n");
      error= TRUE;
    }

       printf("DEBUG>> Read tether request: %s %s %10.6f %10.6f\n",
                       tether.A, tether.B, tether.r0, tether.k);

  return error;
}



void get_restraint_elements(void)
{
#include "header.h"

int i,error;

  num_line_restraints= get_integer(&error, FALSE);

  printf("Number of line restraints= %d\n",num_line_restraints);
  num_line_restraints--;

  for(i=0;i<=num_line_restraints; i++)
    {
       fgets(buffer,BUFFER,input_fp);

       sscanf(buffer, "%s%s%lf", &(line_atom[i].A[0]),
                                 &(line_atom[i].B[0]),
                                 &(line_atom[i].k));

       printf("DEBUG>> Read restraint request: %s %s %10.6f\n",
                       line_atom[i].A, line_atom[i].B, line_atom[i].k);
    }

  return; 
}

void get_atoms_to_fix(void)     
{
#include "constants.h"
#include "header.h"

int i,error, place, done;
int neigh_index;

  mini_fix_atoms=TRUE;
  num_mini_fix= get_integer(&error, FALSE);

  num_mini_fix--;

  if (num_mini_fix <0)
        {
            fprintf(output_fp, "\nERROR: Number of constrained atom types < 1\n");
            exit(1);
        }

/* now read in the constraint sets */
/* NO ERROR CHECKING! */

  for(i=0;i<=num_mini_fix; i++)
      {
        fgets(buffer,BUFFER,input_fp);

        place=0;
        neigh_index=0;
        done = FALSE; 
        
        while (strncmp(&buffer[place]," ",1) == 0 && place < BUFFER ) place++;
        printf("\n place= %d\n", place);
        sscanf(&buffer[place], "%s", &(fix_atoms[i].central_elem[0]));
        place += strlen(&(fix_atoms[i].central_elem[0]))+1;

        while (!done)
          {
             if (buffer[place] != '\0')
               {
                  while (strncmp(&buffer[place]," ",1) == 0 && place < BUFFER ) place++;

                  sscanf(&buffer[place], "%s", &(fix_atoms[i].neigh[neigh_index][0]));
                  place += strlen(fix_atoms[i].neigh[neigh_index])+1;
      
                  if (strcmp(&(fix_atoms[i].neigh[neigh_index][0]),"(") !=0 &&
                      strcmp(&(fix_atoms[i].neigh[neigh_index][0]),")") !=0 &&
                      strcmp(&(fix_atoms[i].neigh[neigh_index][0])," ") !=0 ) neigh_index++; 
               }
             else
               {
                  done = TRUE;
               }
          }
        fix_atoms[i].num_neigh = neigh_index-1;
        printf("\n");
      }

  printf("Will be constraining the following atoms during minimisations:\n");
  for (i=0; i <=num_mini_fix; i++)
    {
      printf("%s with %d neighs: ( ", &(fix_atoms[i].central_elem[0]), fix_atoms[i].num_neigh+1);

      for (neigh_index=0; neigh_index <= fix_atoms[i].num_neigh; neigh_index++) 
                                         printf("%s ",&(fix_atoms[i].neigh[neigh_index][0])); 
      printf(")\n");
    } 

return;
}


void get_forced_dihedrals(void)
{
#include "constants.h"
#include "header.h"

int i,error;

  force_dihedrals=TRUE;
  num_for_dihedrals= get_integer(&error, FALSE);

  printf("Number of fixed dihedrals: %i\n", num_for_dihedrals);
  num_for_dihedrals--;

  if (num_for_dihedrals <0)
		 {
    		     fprintf(output_fp, "\nERROR: Number of set dihedrals < 1\n");
            exit(1);
		 }		 

/* now read in the fixed dihedrals */
/* NO ERROR CHECKING! */

  for(i=0;i<=num_for_dihedrals; i++) 
      {
        printf("Reading buffer\n");
        fgets(buffer,BUFFER,input_fp);
        printf("Scaning buffer\n>>%s<<",buffer);
        sscanf(buffer, "%s%s%s%s%lf", 
		 		 		 		 &(set_diherals[i].A[0]),
		 		 		 		 &(set_diherals[i].B[0]),
		 		 		 		 &(set_diherals[i].C[0]),
		 		 		 		 &(set_diherals[i].D[0]),
                                &(set_diherals[i].phi));

        printf("Printing debug\n");
        set_diherals[i].phi = set_diherals[i].phi/RAD_TO_DEG;
        printf("DEBUG>> Read dihedral to fix >>%s<< >>%s<< >>%s<< >>%s<<  at %10.6f\n",
                set_diherals[i].A, set_diherals[i].B, 
                set_diherals[i].C, set_diherals[i].D, RAD_TO_DEG*set_diherals[i].phi);
      }
return;

}

void get_action_weight(void)
{
#include "header.h"
int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int next_space( int *p_ichar, int start, int num_of_chars );

int iii, itsanum,place, ichar[BUFFER], ndigi, sign;
int iweight, num_of_chars;
int at_end, ifirst;


action_weights_given=1;
itsanum= 1;
place = 0;
num_action_weights= -1;
/*******with end****/

at_end=FALSE;

while (!at_end)
 {
    num_of_chars= read_line(input_fp, &ichar[0]);

/*** Update to allow leading spaces for end statement. DJW Sept 2015 ***/

    ifirst = 0;
    for (iii=0; iii<num_of_chars; iii++) 
      {
         if (ichar[ifirst]==' ') ifirst++;
            else break;
      } 
 
    if (ichar[ifirst]=='e')
      {
        at_end=TRUE;
      }
    else if (num_of_chars == -10 || ifirst == num_of_chars)
      {
        printf("No end statement supplied for action weights input\n");
        exit(0);
      }
    else
      {
        itsanum= TRUE;
        place = 0;
        while (itsanum)
          {
             num_action_weights++;
             action_weights[num_action_weights]= get_int(&ichar[0], &place,
                                             &itsanum, &ndigi, num_of_chars, &sign);

             place= next_space( &ichar[0], place, num_of_chars );
          }
        num_action_weights--;
      }
  }

num_action_weights++;

/****old bit****/
/*
while (itsanum)
   {
    num_action_weights++;
    action_weights[num_action_weights]= get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);

    place= next_space( &ichar[0], place, num_of_chars ); 
   }
*/
/****old bit****/


sum_action_weights= 0;
for (iweight=0; iweight < num_action_weights; iweight++)
   {
     sum_action_weights += action_weights[iweight];
     printf("%d\n",action_weights[iweight]);
   }
printf("Total of weights %d\n",sum_action_weights);
}

void get_hydrogen_weight(int *p_initial_hydrogen_weight,
		 		 		 		 		 		  int *p_increment_hydrogen_weight)
{
#include "header.h"
int read_line(FILE *fp, int *p_ichar);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );
int next_space( int *p_ichar, int start, int num_of_chars );

int itsanum,place, ichar[BUFFER];
int found_num, num_of_chars;

num_of_chars= read_line(input_fp, &ichar[0]);

itsanum= 1;
place = 0;
found_num= -1;

while ((itsanum) && (found_num <1))
   		 {
		 found_num++;
		 if (found_num==0)
		 		 {
		 		 *p_initial_hydrogen_weight = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
		 		 }
		 else
		 		 {
		 		 *p_increment_hydrogen_weight = get_doub(&ichar[0], num_of_chars, 
                                                        &place, &itsanum);
		 		 }
		 place= next_space( &ichar[0], place, num_of_chars );
		 }

return;
}

/***************************************************************/
/**** Get the fragment weights *********************************/
/**** Changed to allow fragment weights to run over     ********/
/**** several lines so that large libraries can be used ********/
/**** 4th July 1998 Dave Willock                        ********/
/***************************************************************/

void get_frag_weight(int *frag_weights, int just_count)
{
#include "header.h"

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int next_space( int *p_ichar, int start, int num_of_chars );

void put_string( int *p_ichar, int length);

int itsanum,place, ichar[BUFFER], ndigi, sign;
int iweight, num_of_chars, at_end, idummy, iii, ifirst;

frag_weights_given=1;
num_frag_weights= -1;
at_end=FALSE;

while (!at_end)
 {
    num_of_chars= read_line(input_fp, &ichar[0]);

/*** Update to allow leading spaces for end statement. DJW Sept 2015 ***/

    ifirst = 0;
    for (iii=0; iii<num_of_chars; iii++) 
      {
         if (ichar[ifirst]==' ') ifirst++;
            else break;
      } 

    if (ichar[ifirst]=='e')
      {
        at_end=TRUE;
      }
    else if (num_of_chars == -10)
      {
        printf("No end statement supplied for fragment weights input\n");
        exit(0);
      }
    else
      {
        itsanum= TRUE;
        place = 0;
        while (itsanum)
          {
             num_frag_weights++;
             if (!just_count)
               {
                  frag_weights[num_frag_weights]= get_int(&ichar[0], &place, 
                                               &itsanum, &ndigi, num_of_chars, &sign);
               }
             else
               {
                  idummy= get_int(&ichar[0], &place, 
                                               &itsanum, &ndigi, num_of_chars, &sign);
               }

             place= next_space( &ichar[0], place, num_of_chars ); 
          }
        num_frag_weights--;
      }
  }

num_frag_weights++;
if (just_count)
  {
     printf("%d Fragment weights counted:\n",num_frag_weights);
  }
else
  {
     printf("%d Fragment weights counted:\n",num_frag_weights);

     sum_frag_weights= 0;
     for (iweight=0; iweight < num_frag_weights; iweight++)
        {
          sum_frag_weights += frag_weights[iweight];
          printf("%d\n",frag_weights[iweight]);
        }
     printf("Total of weights %d\n",sum_frag_weights);
  }
return;
}

/********************************************************************/
/** Routine to read in list of allowed torsions *********************/
/** Added Oct 2006, Dave Willock                *********************/
/********************************************************************/
void get_allowed_torsions(void)
{
#include "constants.h"
#include "header.h"

int i,error;

  num_allowed_torsions = get_integer(&error, FALSE);

  printf("Number of allowed torsions: %i\n", num_allowed_torsions);
  num_allowed_torsions--;

  if (num_allowed_torsions <0)
		 {
    		     fprintf(output_fp, "\nERROR: Number of set allowed torsions < 1\n");
            exit(1);
		 }		 
   else if (num_allowed_torsions > MAX_ALLOWED_TORS)
        {
    		     fprintf(output_fp, "\nERROR: Number of requested allowed torsions\n");
    		     fprintf(output_fp, "         exceeds maximum the program is      \n");
    		     fprintf(output_fp, "         dimensioned for, you asked for %d   \n", num_allowed_torsions);
    		     fprintf(output_fp, "         when the current maximum is %d   \n", MAX_ALLOWED_TORS);
            exit(1);
        }

/* now read in the fixed dihedrals */
/* NO ERROR CHECKING! */

  for(i=0;i<=num_allowed_torsions; i++) 
      {
        printf("Reading buffer\n");
        fgets(buffer,BUFFER,input_fp);
        printf("Scaning buffer\n>>%s<<",buffer);
        sscanf(buffer, "%s%s%s%s", 
		 		 		 		 &(allowed_torsions[i].A[0]),
		 		 		 		 &(allowed_torsions[i].B[0]),
		 		 		 		 &(allowed_torsions[i].C[0]),
		 		 		 		 &(allowed_torsions[i].D[0]));
                               
        allowed_torsions[i].phi=0.0;

        printf("Printing debug\n");
        allowed_torsions[i].phi = allowed_torsions[i].phi/RAD_TO_DEG;
        printf("DEBUG>> Read torsion to allow >>%s<< >>%s<< >>%s<< >>%s<< %10.6f\n",
                allowed_torsions[i].A, allowed_torsions[i].B, 
                allowed_torsions[i].C, allowed_torsions[i].D, RAD_TO_DEG*allowed_torsions[i].phi);
      }
return;

}

/***************************************************************/
/**** Get the seed weights *************************************/
/**** Added to allow seeds in fragment seed option to   ********/
/**** be weighted like fragment choice is               ********/
/**** May 1999 Dave Willock                             ********/
/***************************************************************/

void get_seed_weight(int *seed_weights, int just_count)
{
#include "header.h"

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int next_space( int *p_ichar, int start, int num_of_chars );

void put_string( int *p_ichar, int length);

int itsanum,place, ichar[BUFFER], ndigi, sign;
int iweight, num_of_chars, at_end, idummy;

seed_weights_given=1;
num_seed_weights= -1;
at_end=FALSE;

while (!at_end)
 {
    num_of_chars= read_line(input_fp, &ichar[0]);

    if (ichar[0]=='e')
      {
        at_end=TRUE;
      }
    else if (num_of_chars == -10)
      {
        printf("No end statement supplied for fragment seed weights input\n");
        exit(0);
      }
    else
      {
        itsanum= TRUE;
        place = 0;
        while (itsanum)
          {
             num_seed_weights++;
             if (!just_count)
               {
                  seed_weights[num_seed_weights]= get_int(&ichar[0], &place, 
                                             &itsanum, &ndigi, num_of_chars, &sign);
               }
             else
               {
                  idummy= get_int(&ichar[0], &place, 
                                               &itsanum, &ndigi, num_of_chars, &sign);
               }

             place= next_space( &ichar[0], place, num_of_chars ); 
          }
        num_seed_weights--;
      }
  }

num_seed_weights++;

if (just_count)
  {
    printf("%d Fragment seed weights counted:\n",num_seed_weights);
  }
else
  {
    printf("%d Fragment seed weights read:\n",num_seed_weights);

    sum_seed_weights= 0;
    for (iweight=0; iweight < num_seed_weights; iweight++)
       {
         sum_seed_weights += seed_weights[iweight];
         printf("%d\n",seed_weights[iweight]);
       }
    printf("Total of weights %d\n",sum_seed_weights);
  }

return;
}

/***************************************************************/
/**** Get the name of the minimiser to use *********************/
/***************************************************************/

int get_minimizer_name(char *tok)
{
#include "header.h"
int error;

  error = FALSE;
  if (tok)
    {
       strcpy(minimizer_name, tok);

       if (strcmp(minimizer_name,"discover")   !=0 &&
           strcmp(minimizer_name,"c2discover") !=0 &&
           strcmp(minimizer_name,"mopac")      !=0 &&
           strcmp(minimizer_name,"gulp")       !=0 &&
           strcmp(minimizer_name,"internal")   !=0 )
         {
           fprintf(output_fp,"\nERROR: Unrecognized Miminizer name %s\n", minimizer_name);
           fprintf(output_fp,"Choices are : discover / mopac / internal\n");
           error= TRUE;
         }
       else
         {
           fprintf(output_fp,"Minimiser name : %s\n",minimizer_name);
         }
     }
   else
     {
       fprintf(output_fp,"\nERROR: The minimiser directive (mini) was given but with no arguement\n"); 
       error= TRUE;
     }
   return(error);
}

void get_mopac_name(void)
{
#include "header.h"
		 strcpy(mopac_root, get_string());
		 return;
}

int get_analyse_filename(char *p_analyse_name)
{
#include "header.h"
int error;
  
        error= FALSE;

        strcpy(p_analyse_name, get_string());
		 printf("Analyse run with name : %s\n", p_analyse_name);

/*******************************************************************************/
/***** Check to see if file exists and is readable *****************************/
/*******************************************************************************/

        if (access(p_analyse_name, R_OK) == -1)
                {
                error=TRUE;
                if (errno == 2)
                   {
                       printf("\nERROR: Archive  %s does not exist\n", p_analyse_name);
                   }
                else
                   {
                       fprintf(output_fp,"\nERROR: Archive %s not readable\n", p_analyse_name);
                   }
                }
        return (error);
}

int get_discover_forcefield_name(void)
{
#include "header.h"
int error;

        error= FALSE;
		 strcpy(discover_forcefield_name, get_string());

/**********************************************************************/
/***** check to see if file exists and is readable ********************/
/**********************************************************************/

		 if (access(discover_forcefield_name, R_OK) == -1)
		 		 {
                error= TRUE;
		 		 if (errno == 2)
		 		 		 {
		 		 		 fprintf(output_fp,"\nERROR: Forcefield %s does not exist\n", discover_forcefield_name);
		 		 		 }
		 		 else
		 		 		 {
                        fprintf(output_fp,"\nERROR: Forcefield %s not readable\n", discover_forcefield_name);
		 		 		 }
		 		 }
		 return(error);
}

void get_mopac_path(void)
{
#include "header.h"
		 strcpy(mopac_path, get_string());
    if (access(mopac_path, X_OK) == -1)
       {
       if (errno == 2)
           {
           fprintf(output_fp,"\nERROR: %s does not exist\n", mopac_path);
           }
       else
           {
           fprintf(output_fp,"\nERROR: %s not executable\n", mopac_path);
           }
       fprintf(output_fp,"Exiting....\n");
       fflush(output_fp);
       fflush(stdout);
       exit(EXIT_FAILURE);
       }
		 return;
}

void get_discover_path(void)
{
#include "header.h"
		 strcpy(discover_path, get_string());
    if (access(discover_path, X_OK) == -1)
       {
       if (errno == 2)
           {
           fprintf(output_fp,"\nERROR: %s does not exist\n", discover_path);
           }
       else
           {
           fprintf(output_fp,"\nERROR: %s not executable\n", discover_path);
           }
       fprintf(output_fp,"Exiting....\n");
       fflush(output_fp);
       fflush(stdout);

       exit(EXIT_FAILURE);
		 
       }
   return;
}

int get_animation(char *p_animation_file)

{
#include "header.h"
int error;
/*******16.10.95 Change for multiple animations ********************/

        error= FALSE;
		 strcpy(p_animation_file, get_string());
        if (!*p_animation_file)
          {
             error= TRUE;
             fprintf(output_fp,"\nERROR: No file name given for animation file stem.\n");

          }
		 return(error);
}

void get_mopac_cmds(int whichone)
{
#include "header.h"
/***read in next line for mopac commands***/
fgets(buffer,BUFFER,input_fp);
if (buffer == NULL)
		 {
 		 fprintf(output_fp,
		 "Error Reading Input File. Expecting MOPAC commands. End of File Reached\n");		 
		 exit(-1);
		 return;
		 }

if (whichone ==1)
		 {
		 strcpy(mopac_cmdline_molecule, buffer);
		 }
else if (whichone ==2)
		 {
		 strcpy(mopac_cmdline_inpore, buffer);
		 }
else
		 {
		 fprintf(output_fp,"SERIOUS ERROR in reader! Contact authors 'cos it has gone haywire.\n");
		 }		 		 
return;
}

void get_template_strategy_name(void)
{
#include "header.h"

		 strcpy(template_strategy_file, get_string());
printf("DB>> template_strategy_file %s\n", template_strategy_file);
		 return;

}
void get_inpore_strategy_name(void)
{
#include "header.h"

		 strcpy(inpore_strategy_file, get_string());
printf("DB>> inpore_strategy_file %s\n", inpore_strategy_file);
		 return;

}
		 
void get_defaults_file(void)
{
#include "header.h"

        strcpy(defaults_file, get_string());

        fprintf(output_fp,"DEFAULTS file set to : %s\n", defaults_file);
      
        return;
}


void get_forcefield_library(void)
{
#include "header.h"
		 strcpy(forcefield_library, get_string());
    fprintf(output_fp,"Forcefield library file = %s\n",forcefield_library);
    fprintf(stdout,"Forcefield library file = %s\n",forcefield_library);
    return;
}

void get_template_library(void)
{
#include "header.h"
		 strcpy(gooduns, get_string());
        fprintf(output_fp,"Template Library (Gooduns) file = %s\n",gooduns);
    return;
}

void get_fragment_library(atom *frag_lib, int just_count)
{
#include "header.h"
  int ineigh, i;
  int num_atoms = 0;
  int end_of_frag = 0;

/* referencing fragment library atom array from zero as of 6th April 1995 DW */

  member_start[1] = 0;
  number_of_fragments = 1;
  strcpy(fragment_file, get_string());
  if (!(fragment_fp = fopen(fragment_file, "r")))
        {
                fprintf(output_fp, "ZEBEDDE ERROR: Error Opening Fragment file %s\n", 
                     fragment_file);
                exit(1);
        }
  else
        {
        /* READ IN THE fragment COORDINATES */
        /* get the headers in ; so that I can use again in output */
        fgets(dummy_head,BUFFER,fragment_fp);
        fgets(dummy_head,BUFFER,fragment_fp);
        fgets(fragment_title,BUFFER,fragment_fp); /* hang on to the title */
        fgets(dummy_head,BUFFER,fragment_fp);

/* referencing fragment library atom array from zero as of 6th April 1995 DW */

    i=0;

    while(fgets(buffer,BUFFER,fragment_fp) != NULL)
      {
if (DEBUG) printf("LINE: %s", buffer);
        if (strstr(buffer,"end") != NULL)
          { /*end of car */
             end_of_frag++;

/***********************************************************/
/**** Make fragment library in to an honest to goodness ****/
/**** arc file so that we can put the molecule names in ****/
/**** For zebedde just skip these lines for now         ****/
/**** look at sections with just_one_end variable!!!    ****/
/**** Altered by Dave Willock for SDK 25th Sept. 1996   ****/
/***********************************************************/

             if (end_of_frag == 1)
		 		           {
                   number_of_members[number_of_fragments] = num_atoms-1;
               		    number_of_fragments++;
		 		            total_frag_atoms += num_atoms;
                   num_atoms = 0;
                  }
              else
                  {
                   if (number_of_fragments != 1)
                     {
                      member_start[number_of_fragments] = i;
                     }
                  end_of_frag = 0;
                  fgets(dummy_head,BUFFER,fragment_fp);
if (DEBUG) printf("dumping %s\n",dummy_head);
                  fgets(dummy_head,BUFFER,fragment_fp);
if (DEBUG) printf("dumping %s\n",dummy_head);
                  }
           }
         else
           {
            /*On the counting up call, just count the number of atoms*/
            if(just_count)
             {
               i++;
             }

            /**********************************/
            /*** Read the atom data line ******/
            /**********************************/
           
            else
            {
              sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf", &(frag_lib[i].label[0]), 
      		 		 		                               &frag_lib[i].x,
                                                      &frag_lib[i].y, 
                                                      &frag_lib[i].z,
                                                      &(frag_lib[i].group[0]),
		 		 		                               &(frag_lib[i].group_no[0]),
                                                      &(frag_lib[i].pot[0]),
                                                      &(frag_lib[i].elem[0]),
                                                      &frag_lib[i].part_chge);
  

/* set starting values of neighbour indexing bits of the atom structures */

               frag_lib[i].num_neigh=-1;
               for (ineigh=0; ineigh < 4; ineigh++)
                                      frag_lib[i].neighb[ineigh]=-1;
               num_atoms++;
   
               i++;
             }
             /*Record the number of atoms in the fragments*/
             frag_lib_atno=i;
           }
         }
number_of_fragments--;
  }

if (DEBUG)
 {
  printf("number_of_fragments = %d\n", number_of_fragments);
 for (i=1;i<=number_of_fragments;i++) 
  { 
    printf("Member_start = %d\n",member_start[i]); 
    printf("Number_of_members = %d\n",number_of_members[i]); 
  } 
 }

fclose(fragment_fp);
return;
 }

int get_seed_file(int *p_seed_type)

{
#include "header.h"

int error;

error= FALSE;
strcpy(seed_file, get_string());
if (!(seed_fp = fopen(seed_file, "r")))
    {
    fprintf(output_fp,"\nERROR: Problem opening seed file %s\n", seed_file);
    fflush(output_fp);
    fflush(stdout);
    error= TRUE;
    }
else
		 {
		 fclose(seed_fp);
		 *p_seed_type = MOLE;
		 }

return error;
}
int get_seed_archive(int *p_seed_type, int *p_num_archive_frames_as_seeds,
		 		 		 		 		   int *p_seed_frames_to_read,
		 		 		 		 		   int *p_read_all_frames_as_seeds,
                      int *p_number_of_frames_in_arc) 
{

#include "header.h"

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int string_from_int(int *p_int, char *p_string);

/**********started 25/2/ dwl/Djw**************/

int place, error, ichar[BUFFER];
int itsanum, num_of_chars;
int ndigi, sign, coped_with=FALSE;
char dummy_string[BUFFER];

error= FALSE;
strcpy(seed_file, get_string());

  if (!(seed_fp = fopen(seed_file, "r")))
        {
                fprintf(output_fp, "\nERROR: Problem opening seed archive file %s\n", seed_file);
                fflush(output_fp);
                fflush(stdout);
                error= TRUE;
        }
  else
        {
          while(fgets(buffer,BUFFER,seed_fp) != NULL)
           {
             if (strstr(buffer,"!date") != NULL)
              { 
                (*p_number_of_frames_in_arc)++; 
              }
           }
          (*p_number_of_frames_in_arc)--; 

		 		 fclose(seed_fp);   /**close cos we'll read it elsewhere ***/

        /******* READ IN THE members to be used *******/

		 		 num_of_chars= read_line(input_fp, &ichar[0]);
		 		 place = 0;
		 		 *p_num_archive_frames_as_seeds = 0;

		 		 *p_seed_frames_to_read = get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);

                coped_with= FALSE;
		 		 *p_read_all_frames_as_seeds = FALSE;
		 		 if (itsanum)
		 		 		 {
		 		 		 /*****read in the rest of the frames on the line *****/
                        coped_with= TRUE;
		 		 		 while (itsanum)
		 		 		 		 {
		 		 		 		 p_seed_frames_to_read++;
		 		 		 		 (*p_num_archive_frames_as_seeds)++;

		 		 		 		 *p_seed_frames_to_read= get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);
		 		 		 		 }
		 		 		 }
		 		 else 
		 		 		 {
		 		 		 /*****check for ALL keyword*****/
		 		 		 num_of_chars= string_from_int(&ichar[0], &dummy_string[0]);
		 		 		 if (strstr(dummy_string, "all") || strstr(dummy_string, "ALL"))
                           {
		 		 		 		 *p_read_all_frames_as_seeds = TRUE;
                                 coped_with= TRUE;
                            }
            }

        if (!coped_with)
		 		     {
		 		 		 		 /*******User IQ=0***************/
		 		 		 		 fprintf(output_fp,"\nERROR: No frames specified in seed archive\n");
                                fflush(output_fp);
                                fflush(stdout);
		 		 		 		 error= TRUE;
		 		 		 }
        *p_seed_type = ARCH;
   }
return coped_with;
}
int get_seed_fragment(int *p_seed_type, int *p_num_fragments_as_seeds,
                      int *p_fragments_for_seeds, int *p_random_seed)
{
#include "header.h"

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

int string_from_int(int *p_int, char *p_string);

int read_line(FILE *fp, int *p_ichar);

/**********started 25/2/ dwl/Djw**************/
/*****Re written Oct 99 by DJW to allow ******/
/*****specific Number of Random seeds to be **/
/*****input by the user  *********************/
/*********************************************/

int place,ichar[BUFFER];
int itsanum, num_of_chars;
int ndigi, sign, coped_with, error;
char dummy_string[BUFFER];

/** Read the line under the fragment keyword */

 error=FALSE;
 num_of_chars= read_line(input_fp, &ichar[0]);

/*****  check for RAN keyword            *****/

 num_of_chars= string_from_int(&ichar[0], &dummy_string[0]);

 if (strstr(dummy_string, "ran") || strstr(dummy_string, "RAN"))
   {
      *p_random_seed = TRUE;

/***** Read in the number of seeds ***********/

      *p_num_fragments_as_seeds = get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);
    
      if (itsanum)
        {
          coped_with= TRUE;
        }
      else
        {
          fprintf(output_fp, "ERROR: random key word with fragment seed option must give\n");
          fprintf(output_fp, "       the number of seeds to make for each random choice.\n");
          return TRUE;
        }
   }
 else
   {

/******* Read in the fragments to be used and count them! *******/

      place = 0;
      *p_num_fragments_as_seeds = 0;
      *p_random_seed= FALSE;

      *p_fragments_for_seeds = get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);

      coped_with= FALSE;
      error= FALSE;

      if (itsanum)
        {
/***** read in the rest of the frames on the line *****/
          coped_with= TRUE;
          while (itsanum)
            {
               p_fragments_for_seeds++;
               (*p_num_fragments_as_seeds)++;

               *p_fragments_for_seeds= get_int(&ichar[0], &place, &itsanum,
                                          &ndigi, num_of_chars, &sign);
            }
        }
   }

        if (!coped_with)
          {
/****** User IQ=0 ***********************************************/

            fprintf(output_fp,"\nERROR: No fragments specified as seeds after seed fragment keyword.\n");
             fflush(output_fp);
             fflush(stdout);
             error=TRUE;
          }
 *p_seed_type = FRAG;

return error;
}

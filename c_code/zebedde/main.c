/******************************************************************************/ 
/*                                                                            */
/* Zebedde brought to you by Dewi Lewis and Dave Willock 1995                 */
/* Rejuvenated for Ben Slater 2006                                            */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

#define MAIN 0
#include "header.h"
#include "data.h"
#include "own_maths.h"
#include "ewald.h"
#undef MAIN

/******************************************************************/
/*** Lets not bother about licensing, Dave Willock Nov. 2005 ******/
/******************************************************************/
/* void license(void); */

void analyse_output(char *p_analyse_file, double *p_kvecs, double *p_kvec2,
                     double *p_gvec2, int num_kvecs,
                     double *p_cos_sum, double *p_sin_sum,
                     int *p_need_grad, double *p_grad, atom *pore, 
                     int num_angles_list, int num_torsions_list,
                     int num_vdws_list);

int make_a_template(atom *guest_ptrs[], int *p_num_guests, 
                    list_partition *p_guest_demarc,
                    atom *pore, atom *frag_lib,
                    int *p_frag_weights, 
                    int *guest_hyd_list_ptrs[],
                    int *p_num_guest_hyds,
                    int *guest_hyd_weights_ptrs[],
                    int *p_sum_hyd_weights,
                    atom_number *guest_types_ptrs[], 
                    int *p_num_guest_types,
                    int seed_type,
                    int seed_number,
                    int use_number,
                    list_partition *p_frag_hyd_partition,
                    int *p_frag_hyd_list,
                    atom_number *p_frag_types,
                    list_partition *p_frag_types_list,
                    vec *p_guest_cofm,
                    double *p_total_mass_guests,
                    vec *p_pore_cofm,
                    char *p_animation_file,
                    time_t *p_start_time,
                    int which_test,
                    double *p_average_cc,
                    int test_symm,
                    int increment_hydrogen_weight,
                    int peek_freq,
                    double *p_kvecs, 
                    double *p_kvec2,
                    double *p_gvec2, 
                    int num_kvecs,
                    double *p_cos_sum, 
                    double *p_sin_sum,
                    int *p_need_grad,
                    double *p_grad,
                    double *p_best_energy,
                    int *p_just_position,
                    int have_comb_rules,
                    int is_empty_pore,
                    int *p_have_AB,
                    int want_dock,
                    int is_poly,
                    int num_pore_mols,
                    int *p_num_seed_mols,
                    double cell_volume,
                    double pore_total_mass,
                    int need_monitors,
                    monit *p_monitored,
                    int num_host_mon, 
                    int num_guest_mon);

void print_dashes(int ndashes,FILE *fp);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void setup_defaults(time_t *p_start_time, int *p_test_symm,
                    int *p_which_test, int *p_just_position, int *p_num_actions,
		    int *p_seed_type, int *p_want_dock, int *p_need_monitors,
                    char *p_animation_file );

void enforce_defaults(double *p_num_rot_steps,
                      int *p_num_actions, int *frag_weights);

void initialise_animation(char *p_this_animation_file, char *p_animation_file,
                          int seed_type, int seed_number, int use_number,
                          int build_number);

void initialise_variables( int *frag_weights,
                           int *p_orig_frag_weights,
                           int *p_num_anime_frames,
                           int num_seed_mols,
                           int *guest_hyd_weights_ptrs[],
                           int *p_num_guest_hyds,
                           int *p_sum_hyd_weights,
                           int initial_hydrogen_weight);

void print_biosym_header(FILE *file_fp, int pbc_flag);

void car_end(FILE *fp);

void generate_neighbours( atom *p_molecule, int num_atoms,
                          atom_number *p_types, int *p_num_types,
                          int use_pbc);

void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp);

void print_template_strategy(char *file);

void mdf_end(FILE *fp);

void car_end(FILE *fp);

void nb_print(atom *p_molecule, int num_atoms);

void write_back_input(atom_number *p_frag_types, 
                      atom *pore, atom *frag_lib,
                      int *frag_weights,
                      list_partition *p_frag_types_list,
                      list_partition *p_frag_hyd_partition,
                      int *p_frag_hyd_list, char *p_animation_file,
					  int test_symm,
                      char *p_analyse_name,
                      int randomisation_method, int is_empty_pore);

void banner(FILE *output_fp);

void print_seed_info(int seed_type, int num_of_seeds, int *p_seeds_to_use,
                     int read_all_frames_as_seeds, int random_seed,
                     int number_of_frames_in_arc,
                     int initial_hydrogen_weight, 
                     int increment_hydrogen_weight);

void join_vectors(double *p_A, double *p_B, double *p_A_to_B);

void relabel(atom *p_molecule, int num_to_label, 
                               int *p_num_elems_used, int reset);

void animate(atom *guest_ptrs[], int num_guests,
             list_partition *p_guest_demarc, FILE *an_fp, FILE *lr_fp,
             FILE *ls_fp, char *p_info, char *p_animation_file,
             int *p_frame_counter, atom *pore, int num_seed_mols );

void print_inpore_strategy(char *file);

void regroup(atom *p_molecule, int num_to_group, int group_number);

void cart_latt_vecs( double *p_abc, double *p_latt_vec, double *p_real_latt_sizes, 
                     double *p_recip_latt_vec, double *p_recip_latt_sizes,
                     double *p_cell_volume, int num_atoms, int use_xyz, int use_zyx);

void calculate_box(atom *p_molecule, int num_atoms, double *p_box_limits);

void print_box(double *p_bl, char *box_file);

void sort_out_molecule(atom *p_molecule, int num_atoms, atom_number *p_types,
		       list_partition *p_types_demarc, int *p_hyd_list, int *p_num_hyds, 
                       int use_pbc, vec *p_c_of_m, double *p_total_mass, 
                       int is_pore, int look_for_tethers, int *p_have_AB,
                       int *p_num_mols, list_partition *p_demarc);

void get_seed_molecule(atom *p_original_template,
                       int *p_num_original_template_atoms,
                       int frame_number, int just_count);

void get_pre_anneal_template(atom *p_original_template,
                       int *p_num_original_template_atoms,
                       int frame_number);

double timer(time_t reference_t);

void shift_coords(atom *p_pore, int num_p_atoms);

double real_random(int done);

int box_test( double *p_box_limits, atom *p_molecule, int num_atoms);

int update_partition_list(list_partition *p_list_partition,
                           int index );

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
           int just_count);

void epilogue(char *p_this_animation_file);

void read_potentials(int *p_have_comb_rules);

char *find_zebedde_file(char *file_name, char *executing_file);

void  gen_recip_klist( double *p_kvecs, double *p_kvec2, double *p_gvec2, int *p_num_kvecs,
                       atom *p_host, int num_host_atoms, double *p_cos_sum, double *p_sin_sum);

void order_double_pot(potential_list *p_list_potent, int *p_num_in_list,
                        types *p_pot_types, int *p_num_types, int is_ordered);

int pick_frm_wgt_list( int sum_weights, int *p_weights, int *p_chosen_weight,
                        int num_in_list);

/*******************************************************************/
/*** New Slab processing routines **********************************/
/*******************************************************************/

double neutralise_slab( atom *guest_ptrs[], int num_guests,
                        list_partition *p_guest_demarc, atom *pore);

/*******************************************************************/
/******* Routines stuck in for DEBUG *******************************/
/*******************************************************************/

int neighbour_order(atom *p_molecule, int num_atoms, int atom1, int atom2);

/*******************************************************************/
/******* Main starts here ******************************************/
/*******************************************************************/

  double total_mass_guests[MAX_MOLS], pore_total_mass;
  double num_rot_steps;
  double rdummy;
  double best_energy;
  double cell_volume;

  vec guest_cofm[MAX_MOLS];
  vec *frag_cofm;
  vec pore_cofm;

  int bad_read;
  int success;
  int have_comb_rules;
  int start,ifrag,iloop, jloop, itype, jtype;
  int imol, end;
  int i, ivec, isymm, indx;
  int num_actions;
  int which_test;
  int  num_template_types;
  int num_pore_mols, num_seed_mols, num_frag_mols;
  int keep_order, chosen_seed_weight;
  int have_AB;
  int need_monitors=FALSE;

/**** Variables to do with docking options, Kim Jelfs 2007 ***/

  int temp_store, attempt_store;
  int must_revert_dock, must_revert_centralise;
  int done_pre_annealing;

  FILE *pre_opt_fp;

/*************************************************************/

  double average_cc, tot_charge;

  char discover_root[FILELEN_MAX];
  char animation_file[FILELEN_MAX];
  char *fullstop;
  char info[80];

/******************************************************************/
/*************** Just positioning a molecule variables ************/
/******************************************************************/

int just_position;
int is_slab=FALSE;
int is_poly=FALSE, is_empty_pore=FALSE;

/******************************************************************/
/*************** Original template variables **********************/
/******************************************************************/

  atom *p_original_template;


  /*Make pore a local array that is passed into each function that*/
  /*requires it for memory allocation reasons JL Oct 2007         */
  atom *pore;
  atom dummy_pore[1];

/******************************************************************/
/** New variables for malloced, array of pointers approach to  ****/
/** multiple guests. Dave Willock and Filippo Dec. 2010        ****/
/** guest_ptrs an array of pointers that will be used to hold  ****/
/** the individual molecules that are guests in the host.      ****/
/******************************************************************/

  atom *guest_ptrs[MAX_MOLS];
  atom *p_atom;

  /*Similarly for frag_lib JL Nov 2007*/ 
  atom *frag_lib;
  atom dummy_frag_lib[1];

  /*And frag_weights*/
  int *frag_weights;
  int dummy_frag_weights[1];

  /*And seed_weights*/
  int *seed_weights;
  int dummy_seed_weights[1];

  int num_template_atoms;
  int num_original_template_atoms;
  int *p_original_template_hyd_list, num_original_template_hyds;
  int dummy_list[2], dummy_int;
  int num_original_template_types;
  int num_guest_types[MAX_MOLS];

/******************************************************************/
/****** Original weights data *************************************/
/******************************************************************/

  int orig_frag_weights[MAXTEMPLATE];

  atom_number pore_types[MAX_TYPES],template_types[MAX_TYPES];
  atom_number frag_types[MAX_TYPES];

/*** Pointers for holding seed types, firstly as complete list, then ***/
/*** as individual lists                                             ***/

  atom_number *p_original_template_types;
  atom_number *p_this_orig_type;
  atom_number *guest_types_ptrs[MAX_MOLS];

  list_partition frag_types_list[MAX_TYPES];

/****** symmetry related variables ********************************/

int test_symm; 
int num_temp_atoms_with_symm;

/******************************************************************/

/****** Variables for atom concentration and bond limits **********/

int frag_hyd_list[MAXTEMPLATE], start_hyds,start_types;
int *guest_hyd_list_ptrs[MAX_MOLS], num_guest_hyds[MAX_MOLS];
int *guest_hyd_weights_ptrs[MAX_MOLS], sum_hyd_weights[MAX_MOLS];

list_partition frag_hyd_partition[MAXFRAGMENTS];

/****** Seed  variables *******************************************/
int this_seed, use_number;
int num_of_seeds, read_all_frames_as_seeds;
int seeds_to_use[MAXFRAGMENTS];
int seed_type, number_of_frames_in_arc;
int seed_number, random_seed;

list_partition pore_demarc[MAX_MOLS];
list_partition guest_demarc[MAX_MOLS];
list_partition dummy_demarc[MAX_MOLS];

list_partition *p_seed_types_demarc;
list_partition *p_this_types_demarc;
list_partition dummy_types_demarc;
list_partition pore_types_demarc[MAX_TYPES];

/******************************************************************/

time_t start_time;
time_t end_time;

/****** Hydrogen weighting variables ******************************/

int *p_temp_hyd_weights[MAX_MOLS];
int initial_hydrogen_weight, increment_hydrogen_weight;

/****** variables to allow storing of weighting KJ***************/
int ihyd, iatom, ilink;

/****** Hydrogen weight variables *********************************/

int *p_weight;

int *p_int, *p_num_template_atoms, *p_num_template_hyds;
int *p_sum_temp_hyd_weights, *p_template_hyd_list;

/***end of variable to allow storing of weighting KJ**************/

/****** filename for analyse run ******************************/
char analyse_name[FILELEN_MAX];

/****** Variable for writing out intermediate template ********/
int peek_freq;

/****** Variable controlling method used for randomisation of rnd number generator ****/
        int randomisation_method;

/******Gradient vector for minimiser **************************/
double grad[MAXTEMPLATE3];
int need_grad;

/******************************************************************/
/**** Variables for Ewald Sum *************************************/
/******************************************************************/

double kvecs[MAX_KVECS_ELEM];
double kvec2[MAX_KVECS];
double gvec2[MAX_KVECS];
double cos_sum[MAX_KVECS], sin_sum[MAX_KVECS];

int num_kvecs;

/*** Variables for docking options ********************************/

int want_dock;

/*** Variables for monitoring options *****************************/
monit monitored;
int num_host_mon=-1;
int num_guest_mon=-1;

/******************************************************************/
/*************** Variables for mallocing               ************/
/******************************************************************/
int just_count;

/******************************************************************/
/*********** Exectutable program lines begin here *****************/
/******************************************************************/

/******************************************************************/
/*********** Check we are running on a kosker machine *************/
/*********** FOR SGI ONLY                  DWL 5/2/97 *************/
/******************************************************************/
/* license(); */
/******************************************************************/
/*** Lets not bother about licensing, Dave Willock Nov. 2005 ******/
/******************************************************************/


        /******************************************************************/
        /*********** Check for input file on command line    **************/
        /*********** ready for reading in a while            **************/
        /*********** Was in reader.c but moved 24/20 DWL     **************/
        /*********** To make reading defaults file easier.   **************/
        /******************************************************************/

int main(argc, argv)
  int   argc;
  char  **argv;
{
        if (argc == 1)
           {
              fprintf(stderr,"ZEBEDDE ERROR: Usage ./zebedde input_file\n OR ./zebedde input_file pore_file seed_file\n");
              exit(1);
           }
        else if (argc == 2)
              {
                                strcpy(inputfile,argv[1]);
                                cmd_line_read = FALSE;
                          }
        else if (argc == 4)
                   {
                        cmd_line_read = TRUE;
                        strcpy(inputfile,argv[1]);
                        strcpy(pore_file,argv[2]);
                        strcpy(seed_file,argv[3]);
                        }

                else
                   {
              fprintf(stderr,"ZEBEDDE ERROR: Usage ./zebedde input_file\nOR ./zebedde input_file pore_file seed_file\n");
              exit(1);
           }


                if (!(input_fp = fopen(inputfile, "r")))
           {
             printf("ZEBEDDE ERROR: Error opening file %s\n", inputfile);
             exit(1);
           }

        /********************* Open output file ***************************/

        strcpy(outputfile,inputfile);
        if ((fullstop = strchr(outputfile, '.')) != NULL) *(++fullstop) = '\0';
        strcat(outputfile,"out");

        if ((output_fp = fopen(outputfile,"w")) == NULL)
           {
             fprintf(stderr, "ZEBEDDE ERROR: Unable to open output file %s.\tAborting\n",outputfile);
             exit(1);
           }
        /******************************************************************/
        /*********** Write snazzy banner to output file now  **************/
        /******************************************************************/

        banner(output_fp);

        /******************************************************************/
        /*********** Set up defaults ready to overwrite from **************/
        /*********** input file, do start up things!         **************/
        /******************************************************************/

          setup_defaults( &start_time, &test_symm, &which_test, &just_position,
                                      &num_actions, &seed_type, &want_dock,
                                      &need_monitors,
                                       &animation_file[0]);

          printf("Back from setup defaults\n");

        /******************************************************************/
        /*********** Search for and source zebedde.def       **************/
        /*********** System default file                     **************/
        /*********** So these overwrite wacky ones in setup_defaults ******/
        /******************************************************************/

          strcpy (defaults_file, find_zebedde_file(ZEBEDDE_DEFAULTS, argv[0]));

          if (strcmp(defaults_file, NOT_SET) != 0)
                {
                fprintf(output_fp,
                        "\nDEFAULTS: Default values supplied in %s being applied\n",
                                                                     defaults_file);

                /*** call reader with the defaults file as input file ***/
                /*** First call establishes values for mallocing, second call assigns*/

                just_count=TRUE;
                bad_read= reader(defaults_file, &animation_file[0], dummy_pore, dummy_frag_lib, 
                                 dummy_frag_weights, dummy_seed_weights, &seed_type,
                                 &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                                 &random_seed, &number_of_frames_in_arc,
                                 &initial_hydrogen_weight, &increment_hydrogen_weight,
                                 &test_symm, &analyse_name[0], &peek_freq,
                                 &randomisation_method, &which_test, &just_position,
                                 &is_slab, &is_poly, &want_dock, &is_empty_pore, &need_monitors,
                                 &monitored, just_count);

                if (bad_read)
                  {
                     fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors counting array sizes based on input file.\n");
                     fprintf(output_fp,"Please edit your defaults file before trying again.\n");
                     fflush(stdout);
                     fflush(output_fp);
            /*         exit(0); */
                  }

                 /*Now that we know the number of atom in the host, we can allocate the memory for it*/
                 /*Will generate a symmetry version when we need it now !!!!!! Dave & Filippo 2010   */

                 printf("Malloced pore, for number of atoms, if there is symmetry a new copy will be made.\n");
                 pore=(atom*)malloc((num_pore_atoms+10)*sizeof(atom));

                 frag_lib=(atom*)malloc((frag_lib_atno+10)*sizeof(atom));
                 frag_weights=(int*)malloc((num_frag_weights+10)*sizeof(int));
                 seed_weights=(int*)malloc((num_seed_weights+10)*sizeof(int));

                just_count=FALSE;
                bad_read= reader(defaults_file, &animation_file[0], pore, frag_lib, 
                                 frag_weights, seed_weights, &seed_type,
                                 &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                                 &random_seed, &number_of_frames_in_arc,
                                 &initial_hydrogen_weight, &increment_hydrogen_weight,
                                 &test_symm, &analyse_name[0], &peek_freq,
                                 &randomisation_method, &which_test, &just_position,
                                 &is_slab, &is_poly, &want_dock, &is_empty_pore, 
                                 &need_monitors, &monitored, just_count);

                if (bad_read)
                  {
         fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors in reading section.\n");
         fprintf(output_fp,"Please edit your defaults file before trying again.\n");
         fflush(stdout);
         fflush(output_fp);
            /*         exit(0); */
                  }

                }
          else
                {
                fprintf(output_fp,
                        "\nDEFAULTS: WARNING: System default file (%s) not found\n",
                                                                       ZEBEDDE_DEFAULTS);

                }

/******************************************************************/
/********** Find out what OS we are running on ********************/
/********** Added AJWL 07/08                   ********************/
/******************************************************************/

os_query();

if (os_windows)

fprintf(output_fp,"ZEBEDDE running on Windows\n");

if (os_linux)

fprintf(output_fp,"ZEBEDDE running on LINUX\n");


/******************************************************************/
/*********** read data file (name input by user on the ************/
/*********** command line).                            ************/
/******************************************************************/

  randomisation_method = 0;  /** default is to use real random numbers **/
  peek_freq = 0;
  random_seed= FALSE;
  read_all_frames_as_seeds=FALSE;
  number_of_frames_in_arc= 0;
  strcpy(analyse_name , NO_ANALYSE);

  printf("Going to Count file sizes\n");
  /*** First call establishes values for mallocing, second call assigns*/

  printf("Entering counting phase of reader, expecting ");
  if (is_empty_pore) printf("an empty pore\n");
        else printf("a populated pore\n");

  just_count=TRUE;
  bad_read= reader(inputfile, &animation_file[0], dummy_pore, dummy_frag_lib, dummy_frag_weights, 
                    dummy_seed_weights, &seed_type,
                    &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                    &random_seed, &number_of_frames_in_arc,
                    &initial_hydrogen_weight, &increment_hydrogen_weight,
                    &test_symm, &analyse_name[0], &peek_freq,
                    &randomisation_method, &which_test, &just_position,
                    &is_slab, &is_poly, &want_dock, &is_empty_pore, 
                    &need_monitors, &monitored, just_count);

  printf("Returning from counting phase of reader with ");
  if (is_empty_pore) printf("an empty pore\n");
        else printf("a populated pore\n");

  printf("Found %d pore atoms\n", num_pore_atoms);
  printf("I have counted %d fragment atoms\n", frag_lib_atno);
  printf("Found %d fragment weightings\n", num_frag_weights);
  printf("Found %d seed weightings\n", num_seed_weights);

  if (bad_read)
    {
       fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors counting array sizes based on input file.\n");
       fprintf(output_fp,"Please edit your defaults file before trying again.\n");
       fflush(stdout);
       fflush(output_fp);
       /*  exit(0); */
    }

   /*Now that we know the number of atom in the host, we can allocate the memory for it*/
   printf("Malloced pore for number of atoms found, if there is symmetry this will duplicate the array\n");
   pore=(atom*)malloc((num_pore_atoms+1)*sizeof(atom));

   frag_lib=(atom*)malloc((frag_lib_atno+10)*sizeof(atom));
   frag_weights=(int*)malloc((num_frag_weights+10)*sizeof(int));
   frag_cofm=(vec*)malloc((num_frag_weights+10)*sizeof(vec));
   seed_weights=(int*)malloc((num_seed_weights+10)*sizeof(int));

  just_count=FALSE;

  printf("Entering reader, expecting ");
  if (is_empty_pore) printf("an empty pore\n");
        else printf("a populated pore\n");

  bad_read= reader(inputfile, &animation_file[0], pore, frag_lib, frag_weights, 
                   seed_weights, &seed_type,
                   &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                   &random_seed, &number_of_frames_in_arc,
                   &initial_hydrogen_weight, &increment_hydrogen_weight,
                   &test_symm, &analyse_name[0], &peek_freq,
                   &randomisation_method, &which_test, &just_position,
                   &is_slab, &is_poly, &want_dock, &is_empty_pore, &need_monitors, 
                   &monitored, just_count);

  printf("Returning from reader with ");
  if (is_empty_pore) printf("an empty pore\n");
        else printf("a populated pore\n");

  printf("dat_file read\n");
  fprintf(output_fp, "Cost function will use:\n");

  fprintf(output_fp,"Steric clashes between target and host:                         ");
  if (which_test == STERIC_COST) 
    {
      fprintf(output_fp,"ON\n");
    }
  else
    {
      fprintf(output_fp,"OFF\n");
    }

  fprintf(output_fp,"Just van der Waals interactions between target and host:        ");
  if (which_test == NON_BOND_COST) 
    {
      fprintf(output_fp,"ON\n");
    }
  else
    {
      fprintf(output_fp,"OFF\n");
    }

  fprintf(output_fp,"van der Waals and Coulomb interactions between target and host: ");
  if (which_test == ENERGY_COST) 
    {
      fprintf(output_fp,"ON\n");
    }
  else
    {
      fprintf(output_fp,"OFF\n");
    }

  fprintf(output_fp,"System temperature set as %7.2f K\n\n", temperature);

  if (num_allowed_torsions >=0)
    {
      fprintf(output_fp,"Molecular flexibility allowed for %d torsions as follows:\n",
                                                     num_allowed_torsions+1);
      for (i=0; i<=num_allowed_torsions; i++)
        {
          fprintf(output_fp,"%s - %s - %s - %s\n",allowed_torsions[i].A,
                                                  allowed_torsions[i].B, 
                                                  allowed_torsions[i].C, 
                                                  allowed_torsions[i].D);
        }
    }

 done_pre_annealing = FALSE;
 gulp_input_count =0;

 if (DEBUG)
    {
      fprintf(output_fp,"This is a debugging run, expect lots of output!\n");
      printf("This is a debugging run, expect lots of output!\n");
    }

  if (want_dock) 
    {
      fprintf(output_fp,"This is a docking run\n");
fprintf(output_fp,"The rejection energy is %lf\n\n",dock_energy);
    }
if (exxon) fprintf(output_fp,"Further template modifications will be carried out on Clinton cluster\n");

  must_revert_dock=FALSE;
  must_revert_centralise=FALSE;

if (anneal_run) fprintf(output_fp,"This is a post docking run. Will generate a Discover submission at the end\n");

  if (bad_read)
    {
       fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors in reading section.\n");
       fprintf(output_fp,"Please edit your input data file before trying again.\n");
       fflush(stdout);
       fflush(output_fp);
       exit(0);
    }

       /******************************************************************/
       /*********** If user supplied a defaults_file in the   ************/
       /*********** input file, source it now                 ************/
       /******************************************************************/

  if ((strcmp(defaults_file, NOT_SET) !=0) &&
                             (strcmp(defaults_file, ZEBEDDE_DEFAULTS) !=0))
       {
        if (!(input_fp = fopen(defaults_file, "r")))
           {
             fprintf(output_fp,
                    "\nDEFAULTS: WARNING: Error opening defaults file %s\n",
                                                               defaults_file);
             fprintf(output_fp,"DEFAULTS: WARNING: Skipping defaults file\n");
           }
        else
           {
             fprintf(output_fp,
                     "\nDEFAULTS: Default values supplied in %s being applied\n",
                                                                     defaults_file);
             fprintf(output_fp,
                     "DEFAULTS: Will override values given in input file (%s)\n\n\n",
                                                                         inputfile);

             /*** call reader with the defaults file as input file ***/
             /*** First call establishes values for mallocing, second call assigns*/

             just_count=TRUE;
             bad_read= reader(defaults_file, &animation_file[0], dummy_pore, 
                              dummy_frag_lib, dummy_frag_weights, 
                              dummy_seed_weights, &seed_type,
                              &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                              &random_seed, &number_of_frames_in_arc,
                              &initial_hydrogen_weight, &increment_hydrogen_weight,
                              &test_symm, &analyse_name[0], &peek_freq,
                              &randomisation_method, &which_test, &just_position,
                              &is_slab, &is_poly, &want_dock, &is_empty_pore, &need_monitors, 
                              &monitored, just_count);

             if (bad_read)
               {
                  fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors counting array sizes based on input file.\n");
                  fprintf(output_fp,"Please edit your defaults file before trying again.\n");
                  fflush(stdout);
                  fflush(output_fp);
                  /*  exit(0); */
               }
  
              /*Now that we know the number of atom in the host, we can allocate the memory for it*/
              printf("Malloced pore for number of atoms found, if there is symmetry that will be dealt with\n");
              pore=(atom*)malloc((num_pore_atoms+10)*sizeof(atom));

              frag_lib=(atom*)malloc((frag_lib_atno+10)*sizeof(atom));
              frag_weights=(int*)malloc((num_frag_weights+10)*sizeof(int));
              seed_weights=(int*)malloc((num_seed_weights+10)*sizeof(int));
  
              just_count=FALSE;
              if (DEBUG) printf("Reading defaults\n");
              bad_read= reader(defaults_file, &animation_file[0], pore, frag_lib, 
                               frag_weights, seed_weights, &seed_type,
                               &num_of_seeds, &seeds_to_use[0], &read_all_frames_as_seeds,
                               &random_seed, &number_of_frames_in_arc,
                               &initial_hydrogen_weight, &increment_hydrogen_weight,
                               &test_symm, &analyse_name[0], &peek_freq,
                               &randomisation_method, &which_test, &just_position, &is_slab,
                               &is_poly, &want_dock, &is_empty_pore, &need_monitors, 
                               &monitored, just_count);

              if (bad_read)
                {
                   fprintf(output_fp,"\nZEBEDDE bailing out due to the above errors in reading section.\n");
                   fprintf(output_fp,"Please edit your defaults file before trying again.\n");
                   fflush(stdout);
                   fflush(output_fp);
                   exit(0);
                }
           }
        }

  if (DEBUG)
    {
      printf("DB>> Back in Main.c from reader()\n");
      printf("DB>>initial_hydrogen_weight %d\n",initial_hydrogen_weight);
      printf("DB>>increment_hydrogen_weight %d\n",increment_hydrogen_weight);
    }
  
/******************************************************************/
/*********** Set variables to defaults                 ************/
/******************************************************************/

  enforce_defaults( &num_rot_steps, &num_actions, frag_weights);


/******************************************************************/
/**********  Seed random number generator               ***********/
/**********  0 = seed randomly                          ***********/
/********** -1 = seed at same place each time for check ***********/
/******************************************************************/

  rdummy= real_random(randomisation_method);


/*****************************************************************/
/******* if there is no box yet and its not a periodic system ****/
/******* will need to calculate the box                       ****/
/*****************************************************************/

  if (!pbc)
	{
           if (!user_box) calculate_box(&pore[0], num_pore_atoms, box_limits);
           print_box(box_limits, BOX_OUTPUT_DEFAULT);
           printf("DEBUG>> Back from printing the box\n");
	}

/*****************************************************************/
/******* sort out lattice vectors if required ********************/
/*****************************************************************/

  if (pbc)
    {
       printf("This is periodic\n");

       if ( use_xyz || use_zyx)
          {
            cart_latt_vecs( &abc[0], &latt_vec[0], &real_latt_sizes[0], 
                            &recip_latt_vec[0], &recip_latt_sizes[0], 
                            &cell_volume, num_pore_atoms, use_xyz, use_zyx);
          }
       else
          {
             printf("ERROR: Convention for lattice vectors not supplied.\n");
             printf("       supply convention in input file, either:\n");
             printf("       convention xyz                          \n");
             printf("       or                                      \n");
             printf("       convention zyx                          \n\n");
             printf("       Cerius car files are usually the former. \n");
             exit(0);
          }

/* Added AJWL Dec 08. routine to shift all atoms ********/
/* into the same box ************************************/

shift_coords(&pore[0], num_pore_atoms);


/**** THE BELOW NOW NEEDS TO GO TO THE OUTPUT ROUTINES **/
       fprintf (output_fp,"Cartesian lattice vectors :\n\n");

       ivec=0;
       for (iloop = 0; iloop < 9; iloop++)
         {
            fprintf(output_fp,"%10.6f  ", latt_vec[iloop]);
            if (iloop == 2 || iloop == 5 || iloop == 8) 
		{
                fprintf (output_fp,"  :  size = %10.6f\n",real_latt_sizes[ivec]);
                ivec++;
		}
         }

       fprintf (output_fp,"\nReciprocal lattice vectors :\n\n");

       ivec=0;
       for (iloop = 0; iloop < 9; iloop++)
         {
            fprintf(output_fp,"%10.6f  ", recip_latt_vec[iloop]);
            if (iloop == 2 || iloop == 5 || iloop == 8)
               {
                  fprintf (output_fp,"  :  size = %10.6f\n",recip_latt_sizes[ivec]);
                  ivec++;
               }
         }

  if (is_poly)
      fprintf(output_fp, "Building polymer chains from monomers.\n\n");

  if (is_empty_pore)
      fprintf(output_fp, "Building in an empty periodic cell.\n\n");

  if (is_slab)
      fprintf(output_fp, "This host will be treated as a slab structure\n\n");

/**** THE ABOVE NOW NEEDS TO GO TO THE OUTPUT ROUTINES **/

/********************************************************/
/**** Generate the reciprocal space Ewald sum lists *****/
/********************************************************/
        
       num_kvecs=-1;
       fprintf (output_fp,"\nCell Volume = %10.6f cubic Angstroms\n\n",cell_volume);

       if (!is_empty_pore)
          {
             gen_recip_klist(&kvecs[0], &kvec2[0], &gvec2[0], &num_kvecs,
                                  &pore[0], num_pore_atoms, &cos_sum[0], &sin_sum[0]);

/**** THE BELOW NOW NEEDS TO GO TO THE OUTPUT ROUTINES **/
             fprintf(output_fp,"Maxima for ewald sum: real space= %10.6f, recip space = %10.6f\n\n", 
                                                           real_sum_max,recip_sum_max);
/**** THE ABOVE NOW NEEDS TO GO TO THE OUTPUT ROUTINES **/

          }
    }

/*** Default should be no gradient calculations, Dave Willock May 2013 ***/
need_grad=FALSE;
 
/********************************************************************/
/****** Default intermediate DISCOVER files                    ******/
/****** write out default strategy files if not supplied       ******/
/********************************************************************/

printf("DEBUG>> Off to decide on minimiser\n");

if (strcmp(minimizer_name, DISCOVER_MINIMIZER) ==0 || strcmp(minimizer_name, C2DISCOVER_MINIMIZER) == 0 )
	{
  	if (strcmp(template_strategy_file, NOT_SET) == 0)
			{
			strcpy(template_strategy_file, TEMPLATE_STRATEGY_DEFAULT);
			if (!(strategy_fp = fopen(template_strategy_file,"r")))
				{
				/***no default strategy file so write one****/
				print_template_strategy(template_strategy_file);
                                fprintf(output_fp, "WARNING: No Default template strategy file\n");
                                fprintf(output_fp,"ZEBEDDE has written one to %s\n",
                                                              template_strategy_file);

				}
                        else
    			        { 
                                fclose(strategy_fp);
                                }
			}
else /* template_strategy_file ! = NOT_SET */
         {
        if (!(strategy_fp = fopen(template_strategy_file, "r")))
                {
                fprintf(output_fp,"\nERROR: Could not open Strategy file %s\n", template_strategy_file);
                fflush(output_fp);
                fflush(stdout);
                exit(1);
                }
        /* set up intermediate discover file to the strategy file name root */

        strcpy(template_min_car,template_strategy_file);
        if ((fullstop = strchr(template_min_car, '.')) != NULL) *(++fullstop) = '\0';
        strcpy(template_min_mdf,template_min_car);
        strcat(template_min_car,"car");
        strcat(template_min_mdf,"mdf");
        /* fprintf(output_fp,"template_Minimiser files = %s %s\n", template_min_car, template_min_mdf); */

        fclose(strategy_fp);

         }


  	if (strcmp(inpore_strategy_file , NOT_SET) == 0) 
			{
			strcpy(inpore_strategy_file, INPORE_STRATEGY_DEFAULT);
			if (!(strategy_fp = fopen(inpore_strategy_file,"r")))
                             {
                             /***no default strategy file so write one****/
                             print_inpore_strategy(inpore_strategy_file);
                             fprintf(output_fp, "WARNING: No Default inpore strategy file\n");
                             fprintf(output_fp,"ZEBEDDE has written one to %s\n",
                                                    inpore_strategy_file);

                             }
                        else
                             {
                                fclose(strategy_fp);
                             }

			}

else /* inpore_strategy_file ! = NOT_SET */
         {
        if (!(strategy_fp = fopen(inpore_strategy_file, "r")))
                {
                fprintf(output_fp,"\nERROR: Opening Strategy file %s\n", inpore_strategy_file);
                fflush(output_fp);
                fflush(stdout);
                exit(1);
                }
	/* set up intermediate discover file to the strategy file name root */

        strcpy(inpore_min_car,inpore_strategy_file);
        if ((fullstop = strchr(inpore_min_car, '.')) != NULL) *(++fullstop) = '\0';
        strcpy(inpore_min_mdf,inpore_min_car);
        strcat(inpore_min_car,"car");
        strcat(inpore_min_mdf,"mdf");
        /* fprintf(output_fp,"inpore_Minimiser files = %s %s\n", inpore_min_car, inpore_min_mdf); */

        fclose(strategy_fp);

         }


/*****check that .inp, .car and .mdf are the same : discover crash if not */

	/*get the root of the filename */
	strcpy(discover_root, template_strategy_file);
	if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';

	if ((strncmp(discover_root, template_min_car, strlen(discover_root))!=0)
  	|| (strncmp(discover_root, template_min_mdf, strlen(discover_root))!=0)
  	|| (strncmp(template_min_car, template_min_mdf, strlen(discover_root))!=0))
		{
		fprintf(output_fp,"ZEBEDDE ERROR: Intermediate Discover files must have same root\n");
		fprintf(output_fp,"Files are: %s %s %s\n",template_strategy_file,
			template_min_car,template_min_mdf);
		}

	strcpy(discover_root, inpore_strategy_file);
	if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';

	if ((strncmp(discover_root, inpore_min_car, strlen(discover_root))!=0 )
    	||  (strncmp(discover_root, inpore_min_mdf, strlen(discover_root))!=0)
    	||  (strncmp(inpore_min_car, inpore_min_mdf, strlen(discover_root))!=0))
    	{
    	fprintf(output_fp,
				"ZEBEDDE ERROR: Intermediate Discover files must have same root\n");
    	fprintf(output_fp,"Files are: %s %s %s\n",inpore_strategy_file,
            inpore_min_car,inpore_min_mdf) ;
    	}
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
   /* printf("DEBUG>> Files asscociated with discover: discover_root= %s, inpore_min_car=%s\n",*/
                                                              /* discover_root,inpore_min_car);*/
   /* printf("DEBUG>> Files asscociated with discover: inpore_min_mdf=%s\n", inpore_min_mdf);*/
/**** DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG  DEBUG DEBUG    ***/
}
else if (strcmp (minimizer_name, MOPAC_MINIMIZER) ==0)
	{

        printf("DEBUG>> You asked for mopac minimiser\n");
 	if (strcmp(mopac_cmdline_molecule, NOT_SET) == 0)	
		{
		strcpy(mopac_cmdline_molecule, DEFAULT_MOPAC_COMMANDLINE);
		}

 	if (strcmp(mopac_cmdline_inpore, NOT_SET) == 0)	
		{
		strcpy(mopac_cmdline_inpore, DEFAULT_MOPAC_COMMANDLINE);
		}
	}
else if (strcmp (minimizer_name, GULP_MINIMIZER) ==0)
        {
        printf("DEBUG>> You asked for gulp minimiser, with %s\n", gulp_cmdline_inpore);
        /*printf("gulp_cmdline_inpore before setting is %s, it needs to be %s\n",
                                 gulp_cmdline_molecule, NOT_SET);*/
        if (strcmp(gulp_cmdline_molecule, NOT_SET) == 0)
                {
                printf("Applying default gulp commandline\n");
                strcpy(gulp_cmdline_molecule, DEFAULT_GULP_COMMANDLINE);
                }
        if (strcmp(gulp_cmdline_inpore, NOT_SET) == 0)
                {
                strcpy(gulp_cmdline_inpore, DEFAULT_GULP_COMMANDLINE);
                printf("Applying default gulp commandline\n");
                }
        }

else if (strcmp(minimizer_name,INTERNAL_MINIMIZER) == 0)
        {
             printf("Will use internal minimiser\n");
             /* need_grad=TRUE; */
        }


/******************************************************************/
/**************Set up interaction parameters***********************/
/******************************************************************/
/******************************************************************/
/******* Label pore sequentially and regroup    *******************/
/******************************************************************/
/*** Removed to allow tethering Dave Willock Sept 99  *************/
/******************************************************************/
/*                                                                */
/*  printf("Sending Pore to relabel\n");                          */
/*  relabel(&pore[0], num_pore_atoms, &num_elems_used[0], TRUE);  */
/*                                                                */
/******************************************************************/

  if (non_bonded || charges) 
    {
      if (DEBUG) printf("Reading potentials\n");
      read_potentials(&have_comb_rules);

/**************************************************************************/
/*** reset non-bond index used by atoms according to equivalence list *****/
/**************************************************************************/
   
      if  ( strcmp(pot_info.hbond, H_BOND_AMBER) == 0 )  
         {
             keep_order=FALSE;
             order_double_pot(&h_potent[0], &num_hbonds, &h_pot_types[0], 
                                                     &num_h_pot_types, keep_order);

         }
    }

  if (charges && num_bond_incs >= 0)
   {
/*********************************************************************/
/*** Establish charges using bond increments if set Bond increments  */
/*********************************************************************/
             keep_order= TRUE;
             order_double_pot(&bond_increments[0], &num_bond_incs, 
                                     &bond_inc_types[0], &num_bond_inc_types, keep_order);
   }
	
/******************************************************************/
/************** generate neighbour list for the pore  *************/
/************** (required for writing mdf files)      *************/
/******************************************************************/
  if (!is_empty_pore)
    {
       printf("Sorting out pore file\n");
       sort_out_molecule(&pore[0], num_pore_atoms, &pore_types[0],
                         &pore_types_demarc[0], &dummy_list[0], &dummy_int,
                         pbc, &pore_cofm, &pore_total_mass, TRUE, FALSE,
                         &have_AB, &num_pore_mols, &pore_demarc[0]);

       printf("Number of pore atoms     = %d\n\n",num_pore_atoms);
       printf("Number of pore molecules = %d\n\n",num_pore_mols);

       itype=0;
       for (iloop=0; iloop < num_pore_mols; iloop++)
         {
           printf("    Molecule %d contains %d types\n", iloop+1,
                                                         pore_types_demarc[iloop].num);

           for (jloop=0; jloop < pore_types_demarc[iloop].num; jloop++)
              {
                printf("%s : %d\n",pore_types[itype].atom_type, pore_types[itype].num);
                itype++;
              }
         }
       printf("\nPore total mass          = %10.6f a.m.u.\n\n",pore_total_mass);

       if (is_poly)
         {
           fprintf(output_fp, "Density based on host alone = %10.6f a.m.u. / cubic A\n",
                                                       pore_total_mass/cell_volume);
           fprintf(output_fp, "                           or %10.6f g / cubic cm\n",
                                    1.66*pore_total_mass/cell_volume);
         }


       if (have_tethers) printf("Pore tethered atom is= %s (%d)\n\n", tether.A, tether.index_a);
       print_neighbours(&pore[0], num_pore_atoms, stdout);
    }
  else
    {
       pore_total_mass = 0.0;
       num_pore_mols = 0;
       printf("Pore is an empty periodic cell so no need to generate its neighbours.\n");
    }

/******************************************************************/
/****** Have to do each fragment seperately as neighbours *********/
/****** are referenced independantly within each fragment *********/
/******************************************************************/

  start_types= 0;
  start_hyds=0;

  for (ifrag = 1; ifrag <= number_of_fragments; ifrag++)
     {
       start = member_start[ifrag];

       printf("Sorting out fragment %d, starts at %d\n", ifrag, start);

/*** just send dummy_demarc as fragments should not contain multiple molecules ***/

       sort_out_molecule(&frag_lib[start], number_of_members[ifrag], 
                         &frag_types[start_types], &frag_types_list[ifrag],
                         &frag_hyd_list[start_hyds], &(frag_hyd_partition[ifrag].num),
                         FALSE, &frag_cofm[ifrag], &total_mass_guests[0], FALSE, FALSE, &have_AB,
                         &num_frag_mols, &dummy_demarc[0]);

       printf("Updating hyds partition\n");
       start_hyds=  
             update_partition_list(&frag_hyd_partition[ifrag], start_hyds);

       printf("Updating types partition\n");
       start_types= 
             update_partition_list(&frag_types_list[ifrag],start_types);
    }
/*****************************************************************/
/************Write back the input to the output file**************/
/*****************************************************************/

  printf("Writing back input\n");
  write_back_input(&frag_types[0], pore, frag_lib, frag_weights, &frag_types_list[0],
                   &frag_hyd_partition[0], &frag_hyd_list[0],
                   &animation_file[0], test_symm,
		   &analyse_name[0],
                   randomisation_method, is_empty_pore);
  printf("done - Writing back input\n");

  if (DEBUG)
    {
  printf("DB>> Back in Main.c from write_back_input()\n");
  printf("DB>> mopac path %s\n", mopac_path);
  printf("DB>> discover path %s\n", discover_path);
   } 

/******************************************************************/
/****** print out pore atom non-bonding parameters ****************/
/******************************************************************/

  fprintf(output_fp,"Non-bonding Paramters for Internal Calculations\n");
  fprintf(output_fp,"Parameters read from %s \n", forcefield_library); 	

  printf("done this too!\n");

  if (verbose)
    {

    printf("verbose 1\n");
    fprintf(output_fp,"\nPore Non-bonding parameters\n");
    print_dashes(40,output_fp);

    printf("verbose 1....2\n");
    nb_print(&pore[0], num_pore_atoms);
    printf("verbose 1....3\n");
    fprintf(output_fp,"Read_pot = Potential type given in input file\n");
    fprintf(output_fp,"Use_pot  = Potential type assigned from potentials list\n");
    printf("verbose 2\n");

    fprintf(output_fp,"\n\nFragment Library Non-bonding parameters\n");
    print_dashes(60,output_fp);

    printf("verbose 3\n");
    for (ifrag = 1; ifrag<=number_of_fragments; ifrag++)
        {
        fprintf(output_fp,"\nFragment %i\n",ifrag);
        print_dashes(20,output_fp);
        start =  member_start[ifrag];
        nb_print(&frag_lib[start], number_of_members[ifrag]);
        }
    printf("verbose 4\n");
    fprintf(output_fp,"Read_pot = Potential type given in input file\n");
    fprintf(output_fp,"Use_pot  = Potential type assigned from potentials list\n");
    printf("verbose 5\n");
    }

  printf("and done this too!\n");
/***************************************************************************/
/********* Make a note of the original fragment weights ********************/
/***************************************************************************/
  
   for (iloop=0; iloop < number_of_fragments; iloop++)
      {
         orig_frag_weights[iloop] = frag_weights[iloop];
      }

  printf("done this too as well!\n");

/***************************************************************************/
/***** Checkout that the seeds given exist !! ******************************/
/***************************************************************************/

       if (seed_type == MOLE)
         {
/******************************************************************/
/******** Get the seed if it is from a separate file **************/
/******************************************************************/

            num_of_seeds=1;
            just_count=TRUE;
            get_seed_molecule(p_original_template,
                              &num_original_template_atoms,0, 
                              just_count);

            if (num_original_template_atoms < 0)
              {
                fprintf(output_fp, "ERROR: Seed file contains no atoms, check car and input files\n");
                exit(0);
              }
/**** Now can malloc original_template and go back to fill array ****/

            p_original_template=(atom*)malloc( (num_original_template_atoms+10) * sizeof(atom));
            p_original_template_hyd_list = (int*)malloc( (num_original_template_atoms+10) * sizeof(int));

            if (p_original_template== NULL)
              {
                fprintf(output_fp, "ERROR: Problem assigning memory for seed molecule read\n");
                exit(0);
              }

            just_count=FALSE;
            get_seed_molecule(p_original_template,
                              &num_original_template_atoms,0, 
                              just_count);

            printf("JUST GOT SEED: found %d atoms\n", num_original_template_atoms);
            print_molecule(p_original_template, num_original_template_atoms, stdout, FALSE);
            printf("DEBUG>> Back from print_molecule in main.c\n");

/****************************************************************************/
/**** sort_out_molecule does neighbour generation, finds hydrogen atoms *****/
/**** finds centre of mass and assigns potentials                       *****/
/****************************************************************************/

           printf("Sorting out seed molecule\n");

/*****************************************************************************/
/** Updated Feb 06 to use pbc even with seed molecule as it may be ***********/
/** ripped out of a periodic cell with the host.                   ***********/
/** Dave Willock Feb 06.                                           ***********/
/*****************************************************************************/
          
           printf("Mallocing space %d for types\n",  MAX_TYPES*MAX_MOLS);
           p_original_template_types=(atom_number*)malloc( MAX_TYPES*MAX_MOLS* sizeof(atom_number));
           p_seed_types_demarc= (list_partition*)malloc( MAX_TYPES*MAX_MOLS* sizeof(list_partition));

           if (p_original_template_types == NULL || p_seed_types_demarc == NULL )
             {
               printf("ERROR: Problem with memory assignment for original template types arrays.");
               exit(0);
             }


           sort_out_molecule(p_original_template, num_original_template_atoms,
                             p_original_template_types,
                             p_seed_types_demarc,
                             p_original_template_hyd_list,
                             &num_original_template_hyds, TRUE,
                             &guest_cofm[0], &total_mass_guests[0], FALSE, 
                             TRUE, &have_AB, &num_seed_mols, &guest_demarc[0] );

       printf("\n Summary of demarcation for the %d molecules:\n\n", num_seed_mols);
       printf(" At this point the boundaries still refer to the single list\n");

       p_this_types_demarc= p_seed_types_demarc;
       for (imol=0; imol < num_seed_mols; imol++)
         {
           printf("    Seed molecule %d atoms start at %d end at %d numbering %d\n", imol+1,
                                                                                     guest_demarc[imol].start,
                                                                                     guest_demarc[imol].end,
                                                                                     guest_demarc[imol].num);

           printf("    Seed molecule %d begins %d ends %d contains %d types\n\n", imol+1,
                                                         p_this_types_demarc->start,
                                                         p_this_types_demarc->end,
                                                         p_this_types_demarc->num);
           p_this_types_demarc++;
         }

       printf("\n\nEntering copying loop:\n");

       p_atom= p_original_template;
       p_this_types_demarc= p_seed_types_demarc;

       for (imol=0; imol < num_seed_mols; imol++)
         {
/***** Reserve memory for each guest molecule    ***/
/***** if symmetry is set each molecule will     ***/
/***** be followed by its images so allocate the ***/
/***** space now. Dave Willock March 2013        ***/
/***** These arrays are the same for all symmetry versions so just set up one per ****/
/***** independent molecule                                                       ****/

           guest_hyd_list_ptrs[imol]=(int*)malloc( guest_demarc[imol].num *sizeof(int));
           guest_hyd_weights_ptrs[imol]=(int*)malloc( guest_demarc[imol].num *sizeof(int));

           if (!symm_set) num_symm_ops = -1;

           for (isymm=0; isymm <= num_symm_ops+1; isymm++)
             {
                indx=imol*(num_symm_ops+2) + isymm;
                printf("Mallocing for molecule %d symmetry image %d index %d array with size: %d\n",
                                                        imol, isymm, indx, guest_demarc[imol].num);
                guest_ptrs[indx]=(atom*)malloc( guest_demarc[imol].num *sizeof(atom));

                if (guest_ptrs[indx] == NULL)
                  {
                    printf("ERROR: Problem with memory assignment for guest molecule %d, symmetry image %d,",
                                                                                     imol, isymm);
                    printf(" %d entries long.\n", guest_demarc[imol].num);
                    exit(0);
                  }

/***** Copy guest molecules into new structure ******/

               printf("Copying over atoms to new arrays for molecule %d atoms %d\n",
                                                               imol, guest_demarc[imol].num);

               for (iatom=0; iatom < guest_demarc[imol].num; iatom++)
                 {
                     printf("indx: %d iatom: %d atom: %s\n", indx, iatom, p_atom->label);
                     guest_ptrs[indx][iatom]= *p_atom;
                     p_atom++;
                 }
               guest_demarc[imol].start=0;
               guest_demarc[imol].end=guest_demarc[imol].num-1;
            }

          printf("Copying over types list to new arrays\n");
          guest_types_ptrs[imol]=(atom_number*)malloc( p_this_types_demarc->num * sizeof(atom_number));

          jtype=0;
          p_this_orig_type=p_original_template_types;
          for (itype=p_this_types_demarc->start; itype<=p_this_types_demarc->end; itype++)
            {
               guest_types_ptrs[imol][jtype]= *p_this_orig_type;
               p_this_orig_type++;
               jtype++;
            }
          num_guest_types[imol]=p_seed_types_demarc->num;

           printf("    Reassigned ranges for individual molecule arrays:\n");
           printf("    Seed molecule %d atoms start at %d end at %d numbering %d\n", imol+1,
                                                                                     guest_demarc[imol].start,
                                                                                     guest_demarc[imol].end,
                                                                                     guest_demarc[imol].num);

           printf("    Seed molecule %d begins %d ends %d contains %d types\n", imol+1,
                                                         p_this_types_demarc->start,
                                                         p_this_types_demarc->end,
                                                         p_this_types_demarc->num);

           sort_out_molecule(guest_ptrs[indx], guest_demarc[imol].end,
                             p_original_template_types,
                             &dummy_types_demarc,
                             guest_hyd_list_ptrs[imol],
                             &num_guest_hyds[imol], TRUE,
                             &guest_cofm[imol], &total_mass_guests[imol], FALSE, 
                             TRUE, &have_AB, &dummy_int, &guest_demarc[imol] );

           p_this_types_demarc++;
        }

     printf("Freeing array for original template read from memory.\n");
     free(p_original_template);
     free(p_original_template_types);
     free(p_seed_types_demarc);
     free(p_original_template_hyd_list);

     printf("Original template (i.e. seed) co-ordinates to be used from the separated arrays\n\n");

     printf("Have %d seed molecules\n", num_seed_mols);

      for (imol=0; imol < num_seed_mols; imol++)
        {
           indx=imol*(num_symm_ops+2);
           printf("Molecule %d at index %d\n", imol, indx);
           for (iloop=0; iloop< guest_demarc[imol].num; iloop++)
             {
                printf("%d %s %10.6f  %10.6f  %10.6f\n", iloop,
                                                     guest_ptrs[indx][iloop].label,
                                                     guest_ptrs[indx][iloop].x,
                                                     guest_ptrs[indx][iloop].y,
                                                     guest_ptrs[indx][iloop].z);
            }
          printf("\nCentre of mass : %10.6f  %10.6f  %10.6f  total mass : %10.6f\n\n", 
                                                                     guest_cofm[imol].v[0],
                                                                     guest_cofm[imol].v[1],
                                                                     guest_cofm[imol].v[2],
                                                                     total_mass_guests[imol]);

           printf("Contains %d types: \n", num_guest_types[imol]);

           for (jloop=0; jloop < num_guest_types[imol]; jloop++)
              {
                printf("%s : %d ",guest_types_ptrs[imol][jloop].atom_type, guest_types_ptrs[imol][jloop].num);
              }
           printf("\n\n");
       }

     }
       else if (seed_type == FRAG)
         {
/**** Abandon this option for now ****/

            printf("ERROR: Fragment seeding option not currently available,");
            printf(" please use a separate car file for the seed\n");
            exit(0);

/*************************************/
            if (random_seed)
              {
                for ( this_seed=0; this_seed < num_of_seeds; this_seed++)
                   {
                      seeds_to_use[this_seed]=  1+pick_frm_wgt_list(sum_seed_weights, 
                                                                    &seed_weights[0],
                                                                    &chosen_seed_weight, 
                                                                    num_seed_weights );
                   }
              }
            
            for (this_seed=0; this_seed < num_of_seeds; this_seed++)
               {
                  if (seeds_to_use[this_seed] <0)
                    {
                       printf("ZEBEDDE ERROR: Negative fragment index given as seed\n");
                       exit(EXIT_FAILURE);
                    }
                  else if (seeds_to_use[this_seed] >number_of_fragments)
                    {
                       printf("ZEBEDDE ERROR: Fragment index given as seed exceeds"  
                                               " number of fragments\n");
                       exit(EXIT_FAILURE);
                    }
              }
         }
       else if (seed_type == ARCH)
         {
/* Need to check that none of the indicies exceed the number in arch file*/
           fprintf(output_fp, "ERROR: Archive seed option not yet implemented\n");
           exit(0);
         }

/***************************************************************************/
/***** Print useful info about how the seed is being selected **************/
/***************************************************************************/

      print_seed_info(seed_type, num_of_seeds, &seeds_to_use[0], 
						read_all_frames_as_seeds, random_seed,
						number_of_frames_in_arc,
						initial_hydrogen_weight,
						increment_hydrogen_weight);

      printf("Back from printing seed info...\n");

     if (need_monitors)
       {
          printf("Will monitor distances between adsorbate carbons and host T-sites\n");
          printf("This functionality is under development and only this option is  \n");
          printf("currently allowed.\n");

          printf("Looking to monitor interactions between %s (host) and %s (guest)\n",
                                          monitored.host, monitored.guest);

          printf("ERROR: Monitoring not implemented for multi-guest version\n");
          exit(0);

          if (is_empty_pore)
            {
              printf("ERROR : Cannot monitor host:guest data when host is empty space.\n");
            }
          if (seed_type != MOLE)
            {
              printf("ERROR : Cannot monitor host:guest unless seed is single molecule.\n");
            }

          printf("\nHave total of %d host atoms\n", num_pore_atoms+1);
          printf("Have total of %d guest atoms\n", num_original_template_atoms+1);

/*** Find atoms to monitor in the pore ****/
          num_host_mon = -1;
          for (iloop=0; iloop<= num_pore_atoms; iloop++)
            {
              if (strcmp(monitored.host, pore[iloop].elem) == 0)
                {
                  num_host_mon++;
                  if (num_host_mon < MAX_MONIT)
                    {
                      printf("Monitored pore atom %d is pore atom %d : %s\n", 
                                      num_host_mon, iloop, pore[iloop].label);
  
                      monitored.host_list[num_host_mon]=iloop;
                    }
                  else
                    {
                      printf("ERROR: Cannot cope with more than %d monitored host atoms\n",
                                                MAX_MONIT);
                      exit(0);
                    }
                }
            }

/*** Find atoms to monitor in the guest****/
          p_atom=p_original_template;
          num_guest_mon = -1;
          for (iloop=0; iloop<= num_original_template_atoms; iloop++)
            {
              if (strcmp(monitored.guest, p_atom->elem) == 0)
                {
                  num_guest_mon++;
                  if (num_guest_mon < MAX_MONIT)
                    {
                      printf("Monitored guest atom %d is guest atom %d : %s\n",
                                      num_guest_mon, iloop, p_atom->label);

                      monitored.guest_list[num_guest_mon]=iloop;
                    }
                  else
                    {
                      printf("ERROR: Cannot cope with more than %d monitored guest atoms\n",
                                                MAX_MONIT);
                      exit(0);
                    }
                }
              p_atom++;
            }
       }
     printf("Past monitor set up....\n");
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/*****          CHOOSE Normal run loop or analysis run                ******/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

if (strcmp(analyse_name, NO_ANALYSE) == 0)
{

/***************************************************************************/
/********* Initialise animation file for Good-uns       ********************/
/********* Moved here 4/7/ DWL: previously would ***************************/
/********* overwrite file in analyse mode        ***************************/
/***************************************************************************/


/****** -1 as template number is flag for gooduns! *************************/

   printf("Initialising animation....\n");
   initialise_animation(&these_gooduns[0],  &gooduns[0], seed_type, -1, -1, 0);


/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***** MAIN LOOP Over all the seeds, then over all the maxtemplates   ******/
/***** MAIN LOOP Over all the seeds, then over all the maxtemplates   ******/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

  seed_number=0;
  for (this_seed=0; this_seed < num_of_seeds; this_seed++)
     {
/***************************************************************************/
/********* Loop for the number of templates to build this time *************/
/***************************************************************************/
       printf("Entered seed loop...\n");

       for (iloop=0; iloop < max_templates; iloop++)
          {
                printf("Entered template loop...\n");
                if (want_dock && finished_dock) iloop=max_templates;
                seed_number++;
                use_number=0;
/********************************************************************/
/****** Copy the original template from the use_number fragment *****/
/********************************************************************/

            if ( seed_type != MOLE )
              {
                print_dashes(80, output_fp);
                fprintf(output_fp,"The next %d set of runs will use seed number %d", max_templates, use_number);
                fprintf(output_fp," chosen from the fragment library\n");
                print_dashes(80, output_fp);
              }

            if (have_tethers) fprintf(output_fp,"Seed tethered atom is= %s (%d)\n\n", tether.B, tether.index_b);
            print_dashes(80, output_fp);

            fprintf(output_fp,"Generating template number %d from seed\n", iloop);

            print_dashes(80, output_fp);

            success= FALSE;

/****************************************************/
/** Deal with electrostatics of a slab calculation **/
/****************************************************/

            if (is_slab)
              {
                 printf("Sorting charges for slabs\n");

                 tot_charge= neutralise_slab( guest_ptrs, num_seed_mols,
                                              &guest_demarc[0], pore);

                 printf("Total slab charge now %10.6f\n", tot_charge);
              }

            initialise_variables( frag_weights,
                                  &orig_frag_weights[0], 
                                  &num_anime_frames,
                                  num_seed_mols, 
                                  guest_hyd_weights_ptrs, 
                                  &num_guest_hyds[0],
                                  &sum_hyd_weights[0],
				  initial_hydrogen_weight);

	/*****The way the actions are carried out depends on whether optimising or not ****/
	
	printf("optimise_post_mc is %d\n", optimise_post_mc);
	 
	if(!optimise_post_mc)
	    {
	    printf("A normal run will be used, calling make_a_template at given temp\n");

            printf("number of pore atoms: %d\n", num_pore_atoms+1);
            printf("number of fragment atoms: %d\n",frag_lib_atno);
            printf("number of fragment weights: %d\n",num_frag_weights);
            printf("number of seed atoms: %d\n", num_template_atoms+1);
            printf("number of seed atoms with symmetry (yet to be applied): %d\n", 
                                                             num_temp_atoms_with_symm+1);
            printf("number of seed weights: %d\n",num_seed_weights+1);

            success = make_a_template(guest_ptrs, &num_seed_mols,
                                      &guest_demarc[0],
                                      pore, frag_lib,
                                      frag_weights,
                                      guest_hyd_list_ptrs,
                                      &num_guest_hyds[0],
                                      guest_hyd_weights_ptrs,
                                      &sum_hyd_weights[0],
                                      guest_types_ptrs, 
                                      &num_guest_types[0],
                                      seed_type,
                                      seed_number, 
                                      use_number,
                                      &frag_hyd_partition[0],
                                      &frag_hyd_list[0],
                                      &frag_types[0],
                                      &frag_types_list[0],
                                      &guest_cofm[0],
                                      &total_mass_guests[0],
                                      &pore_cofm,
                                      &animation_file[0],
                                      &start_time, 
                                      which_test,
                                      &average_cc,
                                      test_symm,
                                      increment_hydrogen_weight,
                                      peek_freq,
                                      &kvecs[0], &kvec2[0], &gvec2[0],  
                                      num_kvecs, &cos_sum[0], &sin_sum[0],
                                      &need_grad,
                                      &grad[0],
                                      &best_energy,
                                      &just_position,
                                      have_comb_rules,
                                      is_empty_pore,
                                      &have_AB,
                                      want_dock,
                                      is_poly,
                                      num_pore_mols,
                                      &num_seed_mols,
                                      cell_volume,
                                      pore_total_mass,
                                      need_monitors, 
                                      &monitored,
                                      num_host_mon, 
                                      num_guest_mon);

	    }
	  else
	    {
	      /*******************************************************************************/
	      /**The template will be annealed after the initial make_a_template run**********/
	      /*******************************************************************************/
	      printf("\n\n");
	      printf("You have asked for optimisation post mc, after each normal run\n");
	      printf("the template will be annealed, by gradual reduction of temp\n\n");

	      printf("decrementing by %10.6f K each step\n", mc_temp_step);
	      printf("doing %d actions at each decremented temp\n", mc_modi_step);
		
	      printf("you are starting at %10.6f K\n", temperature);
	      
	      current_temp= temperature;

	      printf(" current_temp is now starting temp %10.6f\n", current_temp);
	      printf("max_templates is %d\n", max_templates);
	      printf("no_annealed_templates so far : %d\n", no_annealed_temps);

	      /********** set so that no docking now - ie dont want zeb to now *************/
	      /********** randomly position seed! and similarly turn of centralising********/
	
	      if(want_dock)
		{
		 printf("You have set docking for each new template - but during annealing docking has defaulted to OFF\n");
		 want_dock = FALSE;
		 must_revert_dock = TRUE;
		}

	      if(centralise_template)
		{
		 printf("You have set template to be centralised for each new template - but during annealing this has defaulted to OFF\n");
                 centralise_template = FALSE;
                 must_revert_centralise = TRUE;
		}
	
		if (done_pre_annealing == FALSE)
		{
		/**** first perform a normal run of make_a_template, at user specified temp ***/
           printf("Off to make a template......\n");
           fflush(stdout);
            success = make_a_template(guest_ptrs, &num_seed_mols,
                                      &guest_demarc[0],
                                      pore, frag_lib,
                                      frag_weights,
                                      guest_hyd_list_ptrs,
                                      &num_guest_hyds[0],
                                      guest_hyd_weights_ptrs,
                                      &sum_hyd_weights[0],
                                      guest_types_ptrs, 
                                      &num_guest_types[0],
                                      seed_type,
                                      seed_number, 
                                      use_number,
                                      &frag_hyd_partition[0],
                                      &frag_hyd_list[0],
                                      &frag_types[0],
                                      &frag_types_list[0],
                                      &guest_cofm[0],
                                      &total_mass_guests[0],
                                      &pore_cofm,
                                      &animation_file[0],
                                      &start_time, 
                                      which_test,
                                      &average_cc,
                                      test_symm,
                                      increment_hydrogen_weight,
                                      peek_freq,
                                      &kvecs[0], &kvec2[0], &gvec2[0],  
                                      num_kvecs, &cos_sum[0], &sin_sum[0],
                                      &need_grad,
                                      &grad[0],
                                      &best_energy,
                                      &just_position,
                                      have_comb_rules,
                                      is_empty_pore,
                                      &have_AB,
                                      want_dock,
                                      is_poly,
                                      num_pore_mols,
                                      &num_seed_mols,
                                      cell_volume,
                                      pore_total_mass,
                                      need_monitors, 
                                      &monitored,
                                      num_host_mon, 
                                      num_guest_mon);

		/**************************************************************************************/	   
		/********* store this structure (will return to this after position optimised) ********/	
		/**************************************************************************************/

/* for(iatom=0; iatom<=num_template_atoms;iatom++) */
/*	{ */
/*	pre_opt_template[iatom]=template[iatom]; */
/*	}	 */
  
/*	if (DEBUG) */
/*	  { */
/*	  printf("KJDB > just stored the template to be used throughout annealing\n"); */
/* 	  printf("KJDB > and no_annealed_temps is %d\n", no_annealed_temps); */
/*	  } */
	
/*	printf("KJ>> template pre annealing\n"); */
 /*             	print_molecule(&template[0], num_template_atoms, stdout, FALSE);	 */
			
/*	done_pre_annealing = TRUE; */

		} /* end if(done_pre_annealing == FALSE) */


		 while(no_annealed_temps <= max_templates)
                {
                if (DEBUG) printf("KJ no_annealed_temps is less than or equal to max_templates %d\n", max_templates);
                anneal_last=FALSE;

		 while(anneal_last==FALSE)
		 {
		  /******** now start decreasing the temperature *********/
		  
		  printf("KJtemp: about to decrease t, current_temp is %10.6f\n", current_temp);

		  current_temp = current_temp - mc_temp_step;

		  printf("KJtemp: have decreased temp, current temp = %10.6f, mc_temp_step = %10.6f\n",
					 current_temp, mc_temp_step);
		  
		  if (current_temp <= 0.0)
		    {
			printf("KJtemp: have found that current_temp is less than or equal to 0 : %10.6f\n",
					 current_temp);
			current_temp = 0.0;
			anneal_last = TRUE;
			no_annealed_temps++;
			if (DEBUG) printf("KJtemp: annealing has reached %10.6f K, annealed_last should be 1 : %d,no_annealed temps increased to %d\n",
				current_temp, anneal_last, no_annealed_temps);
		    }


		/****************************************************************************/
		/*********then running at this decremented temperature **********************/
		/****************************************************************************/
		/*****run make_a_template, but at current_temp and***************************/
		/*****for mc_modi_steps******************************************************/
		/****************************************************************************/

		/****** want to run make_a_template at current_temp, for mc_modi_steps *****/
		/***** so I'll temporarily reset "temperature" and "num_modi_attempts"******/
		
		  temp_store = temperature;
		  attempt_store = num_modify_attempts;
		  num_modify_attempts = mc_modi_step;
		  temperature = current_temp;

            printf("Off to make a template....\n");
            fflush(stdout);

            success = make_a_template(guest_ptrs, &num_seed_mols,
                                      &guest_demarc[0],
                                      pore, frag_lib,
                                      frag_weights,
                                      guest_hyd_list_ptrs,
                                      &num_guest_hyds[0],
                                      guest_hyd_weights_ptrs,
                                      &sum_hyd_weights[0],
                                      guest_types_ptrs, 
                                      &num_guest_types[0],
                                      seed_type,
                                      seed_number, 
                                      use_number,
                                      &frag_hyd_partition[0],
                                      &frag_hyd_list[0],
                                      &frag_types[0],
                                      &frag_types_list[0],
                                      &guest_cofm[0],
                                      &total_mass_guests[0],
                                      &pore_cofm,
                                      &animation_file[0],
                                      &start_time, 
                                      which_test,
                                      &average_cc,
                                      test_symm,
                                      increment_hydrogen_weight,
                                      peek_freq,
                                      &kvecs[0], &kvec2[0], &gvec2[0],  
                                      num_kvecs, &cos_sum[0], &sin_sum[0],
                                      &need_grad,
                                      &grad[0],
                                      &best_energy,
                                      &just_position,
                                      have_comb_rules,
                                      is_empty_pore,
                                      &have_AB,
                                      want_dock,
                                      is_poly,
                                      num_pore_mols,
                                      &num_seed_mols,
                                      cell_volume,
                                      pore_total_mass,
                                      need_monitors, 
                                      &monitored,
                                      num_host_mon, 
                                      num_guest_mon);

		/****** want to return temperature and num_modi_attempts to their original value **/
		temperature = temp_store;
		num_modify_attempts = attempt_store;

		 } /*end of while(anneal_last==FALSE) */

		/*******************************************************************/
		/***Reset to stored structure (in pre_anneal.car)*******************/
		/*******************************************************************/

	
/*for(iatom=0; iatom<=num_template_atoms;iatom++) */
/*	{ */
/*	template[iatom]=pre_opt_template[iatom]; */
/*	}	 */
		
/*printf("KJ>> about to print recovered pre annealed template\n"); */

 /*             print_molecule(&template[0], num_template_atoms, stdout, FALSE); */
  /*            fflush(stdout); */
	

		/****Reset current temp to starting temperature*****************/
		current_temp = temperature;

 		} /*end of while (no_annealed_temps <= max_templates)*/

		if(must_revert_dock)
		 {
		  /***docking was turned off for annealing, turned on again****/
		  printf("Finished annealing template - docking returned to ON\n");
		  want_dock=TRUE;
		 }

                if(must_revert_centralise)
                 {
                  /***centralising was turned off for annealing, turned on again****/
                  printf("Finished annealing template - centralising returned to ON\n");
                  centralise_template=TRUE;
                 }


		if (DEBUG) printf("KJDB, no_annealed_temps: %d, max_templates: %d\n", no_annealed_temps, max_templates);	
		}/*end of if(optimise_post_mc)  */


            printf("Going to epilogue\n");
            epilogue(&animation_file[0]);
            printf("Back from epilogue\n");


           if (success)
             {
                if (which_test == STERIC_COST)
                  {
                     sprintf(info,"Successful Cycle %d, Ave CC %10.6f", 
                                                           iloop, average_cc);
                  }
                else if (which_test == NON_BOND_COST)
                  {
                     sprintf(info,"Successful Cycle %d, Non-bond Energy %10.6f",
                                                               iloop, best_energy);
                  }
                else if (which_test == ENERGY_COST)
                  {
                     sprintf(info,"Successful Cycle %d, Total Energy %10.6f",
                                                               iloop, best_energy);
                  }
              }
           else
             {
                if (which_test == STERIC_COST)
                  {
                     sprintf(info,"Not Successful Cycle %d, Ave CC %10.6f", 
                                                           iloop, average_cc);
                  }
                else if (which_test == NON_BOND_COST)
                  {
                     sprintf(info,"Not Successful Cycle %d, Non-bond Energy %10.6f",
                                                                  iloop, best_energy);
                  }
                else if (which_test == ENERGY_COST)
                  {
                     sprintf(info,"Not Successful Cycle %d, Total Energy %10.6f",
                                                                  iloop, best_energy);
                  }
             }   

        /********* print to gooduns file ****/
        printf("Animating goodun for %d molecules\n", num_seed_mols);
        animate(guest_ptrs, num_seed_mols,
                &guest_demarc[0], gooduns_fp, gooduns_read_fp, gooduns_show_fp, &info[0],
                &gooduns[0], &num_goodun_frames, pore, -1);
         }
      }

} /**** end of normal run loop ******/
else
{
printf("ERROR: Analysis mode not currently implemented\n"); 
/***** analysis run ******/
/* analyse_output(analyse_name, &kvecs[0], &kvec2[0],   */
/*                &gvec2[0], num_kvecs, &cos_sum[0], &sin_sum[0],   */
/*                &need_grad, &grad[0], &pore[0],    */
/*                num_angles_list,  num_torsions_list,   */
/*                num_vdws_list);   */
}

/* print some timing info */
end_time = time(NULL);
fprintf(output_fp, "\n\n");
print_dashes(80,output_fp);
fprintf(output_fp, "Job Finished             : %24.24s\n", ctime(&end_time));
fprintf(output_fp, "Total Job Time           : %f\n", timer(start_time));
print_dashes(80,output_fp);

free(pore);
free(frag_lib);
free(frag_weights);

/*** shut those output files ****/
fclose(output_fp);
/***** and exit gracefully into the night *****/
return 0;

}
	      

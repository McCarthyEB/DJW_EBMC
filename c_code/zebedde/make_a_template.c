#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "own_maths.h"
#include "data.h"

int choose_fragment(atom *p_guest, int num_guest_atoms, int *p_frag_weights,
                    int *p_guest_hyd_list, int num_guest_hyds,
                    int *p_frag_hyd_list, list_partition *p_frag_hyd_partition,
                    atom *frag_lib, atom_number *p_guest_types, int num_guest_types,
                    int *p_num_new_guest_types,
                    atom_number *p_frag_types, list_partition *p_frag_types_list,
                    int have_AB);

void connect_new_fragment(atom *p_guest, list_partition *p_old_demarc, list_partition *p_new_demarc,
                          int *p_place_for_one, int *p_place_for_two,
                          int *p_guest_hyd_list, int *p_num_guest_hyds,
                          int *p_frag_hyd_list, list_partition *p_frag_hyd_partition,
                          atom_number *p_template_types,
                          int *p_num_template_types,
                          int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights,
                          int increment_hydrogen_weight, int *p_have_AB);

void copy_molecule_data( atom *p_B_molecule, atom *p_A_molecule,
                         list_partition *p_B_demarc,  list_partition *p_A_demarc,
                         links *p_B_links, links *p_A_links,
                         int *p_num_B_links, int *p_num_A_links,
                         int *p_B_hyd_list, int *p_A_hyd_list,
                         int *p_B_hyd_weights, int *p_A_hyd_weights,
                         int *p_B_num_hyds, int *p_A_num_hyds,
                         int *p_B_sum_hyd_weights, int *p_A_sum_hyd_weights,
                         atom_number *p_B_types, atom_number *p_A_types,
                         int *p_B_num_types, int *p_A_num_types );

void join_vectors(double *p_A, double *p_B, double *p_A_to_B);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

int  shake_molecule(atom *guest_ptrs[], int num_guests,
                    list_partition *p_guest_demarc,
                    atom *p_pore, double *p_box_limits,
                    interaction_indices *p_molmol_list, 
                    int num_molmol_list,
                    symm_ops *p_symm,
                    int which_test, int num_attempts,
                    double max_move,     
                    double *p_kvecs, double *p_kvec2,
                    double *p_gvec2, int num_kvecs,
                    double *p_cos_sum, double *p_sin_sum,
                    int *p_need_grad, double *p_grad,
                    int have_comb_rules, 
                    int is_empty_pore,
                    int which_mol);

int  rock_molecule(atom *guest_ptrs[], int num_guests,
                   list_partition *p_guest_demarc,
                   atom *p_pore, double *p_box_limits,
                   interaction_indices *p_molmol_list,
                   int num_molmol_list,
                   symm_ops *p_symm, 
                   int which_test, int num_attempts,
                   double max_rot, 
                   double *p_kvecs, double *p_kvec2,
                   double *p_gvec2, int num_kvecs,
                   double *p_cos_sum, double *p_sin_sum,
                   int *p_need_grad, double *p_grad,
                   int have_comb_rules,
                   int is_empty_pore, int which_mol);

void rot_section(atom *guest_ptrs[], int num_guests,
                 list_partition *p_guest_demarc,
                 atom *p_pore, double *p_box_limits,
                 interaction_indices *p_molmol_list, 
                 int num_molmol_list,
                 symm_ops *p_symm,
                 int which_test, links *p_latest_link,
                 int num_to_rot, double max_rot,
                 double *p_kvecs, double *p_kvec2,
                 double *p_gvec2, int num_kvecs,
                 double *p_cos_sum, double *p_sin_sum,
                 int *p_need_grad, double *p_grad,
                 int have_comb_rules, 
                 int is_empty_pore,
                 angle_interact_list *p_angles_list_ptrs[],
                 int *p_num_angles_list, 
                 torsion_interact_list *p_torsions_list_ptrs[],
                 int *p_num_torsions_list, 
                 vdw_interact_list *p_vdw_list_ptrs[],
                 int *p_num_vdws_list,
                 int which_mol);

void set_twist_axis_flags(atom *guest_ptrs[], int num_guests,
                          list_partition *p_guest_demarc,
                          int *p_flag_list,
                          links *p_link_atoms,
                          double *p_origin, double *p_axis, int which_mol);

void do_the_twist(atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc,
                  int *p_flag_list,
                  symm_ops *p_symm,
                  links *p_link_atoms, 
                  double *p_origin, double *p_axis, double *p_theta,
                  int which_mol);

int test_the_twist(atom *p_pore, atom *guest_ptrs[], int num_guests,
                   list_partition *p_guest_demarc,
                   double *p_box_limits,
                   interaction_indices *p_molmol_list,
                   int num_molmol_list,
                   int which_test,
                   double *p_kvecs, double *p_kvec2,
                   double *p_gvec2, int num_kvecs,
                   double *p_cos_sum, double *p_sin_sum,
                   int *p_need_grad, double *p_grad, int have_comb_rules,
                   int is_empty_pore,
                   angle_interact_list *p_angles_list_ptrs[],
                   int *p_num_angles_list,
                   torsion_interact_list *p_torsions_list_ptrs[],
                   int *p_num_torsions_list,
                   vdw_interact_list *p_vdw_list_ptrs[],
                   int *p_num_vdws_list,
                   int which_mol);

	void relabel(atom *p_molecule, int num_to_label,
						    int *p_num_elems_used, int reset);

	void animate(atom *guest_ptrs[], int num_guests,
		     list_partition *p_guest_demarc, FILE *an_fp, FILE *lr_fp, 
		     FILE *ls_fp, char *p_info, char *p_animation_file,
		     int *p_frame_counter, atom *pore);

	void minimize_inpore(atom *p_pore, atom *p_molecule,
			     int num_atoms, char *root_name,
			     double *p_kvecs, double *p_kvec2,
			     double *p_gvec2, int num_kvecs,
			     double *p_cos_sum, double *p_sin_sum,
			     int *p_need_grad, double *p_grad,
			     int have_comb_rules, int is_empty_pore);

	void print_inpore_strategy(char *file);

	void regroup(atom *p_molecule, int num_to_group, int group_number);

	int test_pore_symmetry(FILE *fp, atom *p_pore, int num_atoms,
			       symm_ops *p_symm, int num_symm_ops);

	void do_symmetry_all(atom *guest_ptrs[], int num_guests, 
			     list_partition *p_guest_demarc, 
			     symm_ops *p_symm, int num_sym_ops);

	void copy_types( atom_number *p_to, int *num_to,
			 atom_number *p_from, int num_from );

	void add_types( atom_number *p_types_to, int *p_num_to,
			atom_number *p_types_from, int num_types_from);

	void calculate_energy(atom *p_pore, int num_p_atoms,
			      atom *guest_ptrs[], int num_guests,
			      list_partition *p_guest_demarc,
			      interaction_indices *p_molmol_list,
			      int num_molmol_list,
			      double *p_kvecs, double *p_kvec2,
			      double *p_gvec2, int num_kvecs,
			      double *p_cos_sum, double *p_sin_sum,
			      int *p_need_grad, double *p_grad,
			      int have_comb_rules, int is_empty_pore);

	void calculate_intra_energy( atom *guest_ptrs[], int num_guests,
				     list_partition *p_guest_demarc,
				     angle_interact_list *p_angles_list_ptrs[],
				     int *p_num_angles_list,
				     torsion_interact_list *p_torsions_list_ptrs[],
				     int *p_num_torsions_list,
				     vdw_interact_list *p_vdw_list_ptrs[],
				     int *p_num_vdws_list);

	void print_dashes(int ndashes,FILE *fp);

	void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name,
			       int *p_need_grad, double *p_grad, int is_empty_pore);

	void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp);

	void print_neighbours_heavies( atom *p_molecule, int num_atoms, FILE *fp);

	double timer(time_t reference_t);

	double real_random(int done);

	int box_test( double *p_box_limits, atom *p_molecule, int num_atoms);

	int pick_frm_wgt_list( int sum_weights, int *p_weights, int *p_chosen_weight,
			       int num_in_list );

	int cost_function(atom *guest_ptrs[], int num_guests,
			  list_partition *p_guest_demarc,
			  atom *p_pore, double *p_box_limits, 
			  interaction_indices *p_molmol_list,
			  int num_molmol_list,
			  int which_test, double *p_cc_tot, 
			  int flags_set, int *p_flag_list,
			  double *p_kvecs, double *p_kvec2,
			  double *p_gvec2, int num_kvecs,
			  double *p_cos_sum, double *p_sin_sum,
			  int *p_need_grad, double *p_grad,
			  int have_comb_rules, int is_empty_pore,
			  angle_interact_list *p_angles_list_ptrs[],
			  int *p_num_angles_list, 
			  torsion_interact_list *p_torsions_list_ptrs[],
			  int *p_num_torsions_list, 
			  vdw_interact_list *p_vdw_list_ptrs[],
			  int *p_num_vdws_list,
			  int need_intra);

	void print_energy(FILE *fp, energy *p_energy, internal_energy *p_intra, atom *p_molecule, 
						int num_atoms, int madalung, int have_built,
						int need_mini_ene, int guest_no, int build_no);

	int update_partition_list(list_partition *p_list_partition, int index );

	void initialise_animation(char *p_this_animation_file, char *p_animation_file, 
				  int seed_type, int seed_number, int use_number, 
				  int build_number);

	int ring_maker(atom *p_molecule, int *p_num_atoms, bond *p_forbidden_bond,
		       int *p_template_hyd_list, int *p_num_template_hyds,
		       int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights);

	void print_statistics(stats *p_statistics, int number_of_actions,
				     int *p_action_weights, FILE *output_fp);

	int neighbour_order(atom *p_molecule, int num_atoms, int atom1, int atom2);

	void print_peek_file(atom *guest_ptrs[], int num_guests, list_partition *p_guest_demarc, 
			     int counter, FILE *peek_fp);

	void find_line_atoms(int line, atom *p_mol, int num_atoms, int *p_iA, int *p_iB);

	void unit_vector(double *p_vector, double *p_size);

	void write_pdb(FILE *pdb_fp, atom *guest_ptrs[], int num_guests,
		       list_partition *p_guest_demarc, double *p_abc,
		       int *p_super, double *p_latt_vec);

	int open_file(FILE **p_file, char *p_filename, char *p_status);

	void print_dist_mat( atom *p_molecule, int num_atoms, int use_pbc, FILE *fp);

	int assemble_molmol_list(interaction_indices *p_molmol_list, 
				 int num_guests, int num_expected);

	void assemble_intra_lists(atom *guest_ptrs[], int num_guests, list_partition *p_guest_demarc, 
				  int num_symm_ops, angle_interact_list *p_angles_list[],
				  torsion_interact_list *p_torsions_list[], vdw_interact_list *p_vdw_list[],
				  links *p_links_list_ptrs[],
				  int *p_num_angles_list, int *p_num_torsions_list, int *p_num_vdws_list, 
				  int *p_num_links, int just_count);

	void min_image( double *x, double *y, double *z);

	void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
			    int num_atoms, int which_mol );

void dock_orientate(atom *p_template, int num_template_atoms,
                    double *p_cofm);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc);

void assign_stretch(atom *p_molecule, int num_atoms);

void generate_input(int number_of_actions, int *p_action_weights, 
                    int *p_frag_weights, int num_frag_weights);

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
                    int num_guest_mon)
{
#include "header.h"

  double temp_to_pore[3];
  double cc_tot, rnmol;

  int iii, jjj, kkk, animated, idummy, imol;
  int iloop, ist, ied, mem_prob; 
  int need_test=FALSE, i, iatom, jatom;
  int action=-1,have_built, have_made_ring, have_minimized, animate_this;
  int num_elems_used[NUM_ELEMENTS];
  int reject_template, finished;
  int dummy_int,ilink;
  int last_frag_added;		 
  int weight_of_chosen_action;
  int which_mol, which_built_mol, which_frag, num_action_failed, can_do;
  int *p_int, build_number, last_build_number;
  int *p_num_this_types;
  int place_for_one, place_for_two;   /* indices to retain where build removed H atoms */

  int fragments_in_template;    /* for labelling residues for discover */
  int done_move; /* for flaging that shake, rock etc has been successful */
  int *p_flag_list; /* carries lists of flagged atoms, for example for setting up twist rotation */

  char discover_root[FILELEN_MAX];
  char *fullstop;
  char this_animation_file[FILELEN_MAX];
  char pdb_filename[FILELEN_MAX];
  char xyz_filename[FILELEN_MAX];
  char csv_filename[FILELEN_MAX];
  char debug_filename[FILELEN_MAX];
  char energy_csv_filename[FILELEN_MAX]; /*** csv file for recording energies ***/

  char info[BUFFER];
  char command[FILELEN_MAX];

  int num_old_links;
  links *p_old_links, *p_this_link;

/****** Twist action variables     ********************************/
  int num_to_rot, reject_twist, num_rots;
  int allow, isign;

  double origin[3], axis[3], theta;
/****** symmetry related variables ********************************/

  double min_dist, max_image_dist;
/******************************************************************/
/*** Molecule demarcation *****************************************/
/******************************************************************/
int start, num;
list_partition *p_demarc; 
list_partition *p_old_demarc, *p_which_demarc;

/****** Variables for atom concentration and bond limits      *****/
/****** Keeping enough information to restore after rejection *****/
/****** Dave Willock.                                         *****/

int *p_old_hyd_list, num_old_hyds;
int ihyd, itype;

atom *p_old_guest, *p_old_atom;

/****** Hydrogen weight variables *********************************/

int *p_old_hyd_weights, sum_old_hyd_weights;
int *p_weight;

/****** Variables for Statistics on run ***************************/
 
 stats statistics[NUMBER_OF_ACTIONS];
 int ifrag, composition[MAXFRAGMENTS];

/****** peek  variables *******************************************/
  int peek_now, peek_now_energy;
  int anim_pdb_isopen=FALSE; /***  keep track of anim_pdb file statis ***/
  int energy_file_is_open=FALSE;
  
  FILE *peek_fp; /* fp for peek file; set in defaults.h */
  FILE *anim_pdb_fp; /*** file pointer for pdb version of animation frames ***/
                     /*** added Feb 2006 Dave Willock                      ***/
  FILE *anim_xyz_fp; /*** file pointer for xyz version of animation frames ***/
                     /*** added Nov 2013 Dave Willock                      ***/
  FILE *debug_fp;    /*** file pointer for writing debug files.            ***/
                     /*** added Jan 2014 Dave Willock                      ***/
  FILE *csv_fp; /*file pointer for energies in csv format*/
  FILE *energy_csv_fp; /*file pointer for energies in csv format*/


  time_t loop_time;

  atom *p_atom, *p_guest;
  int *p_this_frag_H_list, start_H_list;
  list_partition *p_this_frag_hyd_partition;
/* atom best_ever_template[MAXTEMPLATE]; */

/*****Added by DJW 28th Aug 1997 to sort out bug in ring maker *****/
/*****Updating its types properly                              *****/

  atom_number *p_type; 

/***** Atom type information to allow reinstallation of guest after an ***/
/***** unsuccessful build.                                             ***/
/***** Dave Willock, Dec. 2013                                         ***/

  atom_number *p_old_types;
  int num_new_types, num_old_types;
 
  vec *p_this_cofm, *p_which_cofm;
  double *p_this_totmass, *p_which_totmass;
/******************************************************************/
/*** Added to cope with energy based cost functions ***************/
/*** Dave Willock March 1997 updates for multiple molecules 2013 **/
/******************************************************************/

  double total_energy, total_inter, total_intra, inter_term, intra_term;
  double best_restraint=0.0;
  double best_intra, tot_non_bond, tot_hbond, tot_charges;
  double best_intra_vdw;
  double best_total;
  double best_ever;
  double boltz_fact;
  double vec[3], vec1[3], shft_vec[3];
  double total_mass_temp, siz;

  interaction_indices *p_molmol_list;
  int num_molmol_list;
                         
/******************************************************************/
/*** Variables for MC, Dave Willock Aug 2009 **********************/
/******************************************************************/

  int entry_no, num_mc_accepted, num_mc_recorded;

  double energy_tot_grand,  energy_tot_grand2; 
  double energy_inter_grand,  energy_inter_grand2; 
  double av_etot, av_einter, stddev_etot, stddev_einter;
/******************************************************************/
/*** Added to cope with restraints DJW July 1998    ***************/
/******************************************************************/

  int num_new;
  int irestr, iA, iB;
  atom *p_A, *p_B;
  double line_ref[3],size;

/******************************************************************/
/*** Variables for pdb files **************************************/
/******************************************************************/
  int super[3];

/******************************************************************/
/*** DEBUG variables **********************************************/
/******************************************************************/

  int debug_loop, isymm;
  atom *p_this_molecule;


/******************************************************************/
/*** Variables for monitoring *************************************/
/******************************************************************/

  int host_mon, guest_mon, host_ind, guest_ind;
  double dist2, dist, dx, dy, dz;

  atom *pore_image;

  int just_count;
  int index, index_image, which_index;
  int mols_written;
  int num_seeds_with_symm;
  int num_angles_list[MAX_MOLS], num_torsions_list[MAX_MOLS], num_vdws_list[MAX_MOLS];
  int num_links_list[MAX_MOLS];
  int num_angles_counted[MAX_MOLS], num_torsions_counted[MAX_MOLS], num_vdws_counted[MAX_MOLS];
  int num_links_counted[MAX_MOLS];

/******************************************************************/
/**** For multiple molecule additions make an array of pointers ***/
/**** so that each molecule in the list can have its own set of ***/
/**** intra-molecular interaction lists.                        ***/
/******************************************************************/

  angle_interact_list *p_angles_list_ptrs[MAX_MOLS];
  torsion_interact_list *p_torsions_list_ptrs[MAX_MOLS]; 
  vdw_interact_list *p_vdw_list_ptrs[MAX_MOLS];
  links *p_links_list_ptrs[MAX_MOLS];

/******************************************************************/
/**** Variables for GCMC ******************************************/
/******************************************************************/

  int have_created, have_removed;

/******************************************************************/
/*********** Exectutable program lines begin here *****************/
/******************************************************************/

if (DEBUG)
  {
    printf("\n\n\nArrived in make_a_template with DEBUG on\n\n\n");
  }
else
  {
    printf("\n\n\nArrived in make_a_template with DEBUG off\n\n\n");
  }

/**** Set the scene for mallocing the molecule arrays *************/
/**** Code combined Feb 2013                          *************/

    printf("\nUsing demarcation to malloc space for guests:\n");
    p_demarc=p_guest_demarc;

    *p_num_guests = *p_num_seed_mols;

    if (symm_set)
      {
        if (*p_num_seed_mols *( num_symm_ops + 1) > MAX_MOLS)
         {
            printf("ERROR: The product of the number of seed molecules (%d) ", 
                                             *p_num_seed_mols);
            printf("and number of symmetry operations (%d) supplied\n", num_symm_ops + 1);
            printf("          is greater than MAX_MOLS (%d).\n", MAX_MOLS);
            exit(0);
         }
      }
    else 
      {
        if (*p_num_seed_mols > MAX_MOLS)
         {
            printf("ERROR: The number of molecules supplied in the seed");
            printf(" (%d) is greater than MAX_MOLS (%d).\n",
                                             *p_num_seed_mols , MAX_MOLS);
            exit(0);
         }
       num_symm_ops=-1;
      }

/*** value for highest index molecule for this set of seeds ***/
       num_seeds_with_symm= *p_num_seed_mols*( num_symm_ops + 2)-1;
       printf("num_seed_mols %d, num_symm_ops %d so num_seeds_with_symm=%d\n", 
                              *p_num_seed_mols, num_symm_ops, num_seeds_with_symm); 

    p_this_cofm=p_guest_cofm;
    p_this_totmass= p_total_mass_guests;
    p_demarc=p_guest_demarc;
    for (imol=0; imol < *p_num_seed_mols; imol++)
      {
         printf("mol %d starts %d ends %d num %d\n", imol, p_demarc->start, 
                                                          p_demarc->end, 
                                                          p_demarc->num);

         printf("\nCentre of mass : %10.6f  %10.6f  %10.6f total mass: %10.6f\n\n",
                                                               p_this_cofm->v[0], 
                                                               p_this_cofm->v[1], 
                                                               p_this_cofm->v[2], 
                                                               *p_this_totmass);

        if (symm_set) index=imol*(num_symm_ops+2);
                                         else index=imol;

        printf("\nNeighbour information for the %d atoms of guest molecule %d\n", p_demarc->num, imol);

        print_neighbours(guest_ptrs[index], p_demarc->end, stdout);
        p_demarc++;
        p_this_cofm++;
        p_this_totmass++;
     }

/**** If needed zero monitoring arrays *****************/
if (need_monitors)
  {
    printf("DEBUGING FOR MONITOR FUNCTION\n");
    printf("For the multi-molecule version\n");
    printf("Monitoring currently disabled.\n");
    printf("DEBUGING FOR MONITOR FUNCTION\n");
    exit(0);
    for (guest_mon=0; guest_mon<=num_guest_mon; guest_mon++)
      {
        guest_ind = p_monitored->guest_list[guest_mon];

      /*  printf("%d guest to monitor: %d %s \n", guest_mon, guest_ind, (p_template+guest_ind)->label); */
      }

    for (host_mon=0; host_mon<=num_host_mon; host_mon++)
      {
        host_ind = p_monitored->host_list[host_mon];

        printf("%d host to monitor: %d %s \n", host_mon, guest_ind, pore[host_ind].label);

        for (guest_mon=0; guest_mon<=num_guest_mon; guest_mon++)
          {
            guest_ind = p_monitored->guest_list[guest_mon];

            p_monitored->entry[host_ind][guest_ind]=0;
          }
      }
   }

/****************************************************************************/
/*** Zero composition counters **********************************************/
/****************************************************************************/
for (ifrag= 1; ifrag <= number_of_fragments; ifrag++) composition[ifrag]=0;

/****************************************************************************/
/*********************** The meat of the program ****************************/
/*********************** The meat of the program ****************************/
/****************************************************************************/

if (which_test == STERIC_COST) printf("will do steric cost_function\n");
if (which_test == NON_BOND_COST) printf("will do non-bond cost_function\n");
if (which_test == ENERGY_COST) printf("will do energy cost_function\n");

/****************************************************************************/
/*********** Initialise variable for peeking at the template ****************/
/****************************************************************************/
peek_now = peek_freq;
peek_now_energy = peek_freq/10;
if (peek_now_energy == 0 ) peek_now_energy = 1;

/***WHY????*****/
cc_tot=500;
best_total = 100000.0;
best_intra = 100000.0;
best_intra_vdw = 100000.0;
best_ever = 100000.0;

animated=FALSE;
/****************************************************************************/
/*********** Initialise statistics vars for this template *******************/
/****************************************************************************/

for (i=0; i<  NUMBER_OF_ACTIONS; i++)
  {
    statistics[i].tries    = 0;
    statistics[i].accepted = 0;
    statistics[i].tries_after_build    = 0;
    statistics[i].accepted_after_build = 0;
  }
 
printf("Initialised run statistics in make_a_template\n");


/****************************************************************************/
/*********** Initialise animation files for this template                 ***/
/****************************************************************************/

build_number=0;
if (animate_flag != 0)
  {
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
    if (use_number > 0)
      {
        printf("ERROR: This call to make_a_template is not possible, use_number = %d\n", use_number);
        exit(0);
      }
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
    initialise_animation(&this_animation_file[0], p_animation_file, 
                         seed_type, seed_number,use_number,
                         build_number);
/*************************************/
/** Additions for pdb file writer ****/
/*************************************/

    if (seed_type != MOLE)
      {
         sprintf(pdb_filename,"%s%d_%d_%d.pdb",p_animation_file,
                                       use_number, seed_number, build_number);

         sprintf(xyz_filename,"%s%d_%d_%d.xyz",p_animation_file,
                                       use_number, seed_number, build_number);
      }
    else
      {
         sprintf(pdb_filename,"%s%d_%d.pdb",p_animation_file,
                                                seed_number, build_number);

         sprintf(xyz_filename,"%s%d_%d.xyz",p_animation_file,
                                                seed_number, build_number);
      }

    if (open_file(&anim_pdb_fp, pdb_filename, "w") == EXIT_FAILURE)
      {
         printf("ERROR>> Cannot open pdb file %s for writing\n", pdb_filename);
         exit(0);
      }

    if (open_file(&anim_xyz_fp, xyz_filename, "w") == EXIT_FAILURE)
      {
         printf("ERROR>> Cannot open xyz file %s for writing\n", xyz_filename);
         exit(0);
      }

    anim_pdb_isopen=TRUE;
  
    super[0]=1;
    super[1]=1;
    super[2]=1;
  }
    /************************************************/
    /***Create a csv file for this attempt JL Ot 07**/
    /************************************************/
    if (animate_flag !=0)
      {
        strcpy(csv_filename, pdb_filename);
        *strrchr(csv_filename, '_')='\0';
        strcat(csv_filename, ".csv");

        if (open_file(&csv_fp, csv_filename, "w")== EXIT_FAILURE)
          {
             printf("ERROR>> Cannot open csv file %s for writing\n", csv_filename);
             exit(0);
          }

        strcpy(energy_csv_filename, "total_energy.csv");

        printf("Will put Energies in csv format into file: %s\n", energy_csv_filename);  

        if (open_file(&energy_csv_fp, energy_csv_filename, "w")== EXIT_FAILURE)
          {
             printf("ERROR>> Cannot open energy csv file %s for writing\n", energy_csv_filename);
             exit(0);
          }
        energy_file_is_open=TRUE;
        fprintf(energy_csv_fp, "iteration, intra, inter, total (kcal /mol)\n");
      }


  if (DEBUG) printf("DEBUG >> Will put Energies in csv format into file: %s\n", csv_filename);  

  if (DEBUG) printf("DEBUG >> Have %d pore mols and %d seed mols\n",
                                     num_pore_mols, *p_num_seed_mols);
/***************************************************************************/
/**** Set reference for constraints which hold template in initial *********/
/**** orientation DJW July 1998                                    *********/
/***************************************************************************/

 if (line_restraints && line_hold)
   {

      printf("ERROR: line_restraints not currently available in multi-guest version\n");
      exit(0);

      if (DEBUG) 
       {
         printf("Setting up line restraints\n");
         printf("DEBUG>> line_restraints = %d, line_hold= %d\n", 
                                        line_restraints, line_hold);
       }
      for (irestr=0; irestr <= num_line_restraints; irestr++)
        {
/*         find_line_atoms(irestr, p_template, *p_num_template_atoms, &iA, &iB); */

           if (iA == iB)
             {
               printf("ERROR : Could not find line atoms in seed check input\n");
               exit(0);
             }

/*         p_A = p_template+iA; */
/*         p_B = p_template+iB; */

           printf("DEBUG>> Setting line restraint between atoms %d : %s and %d : %s\n",
                          iA, p_A->label, iB, p_B->label);

           line_ref[0] = p_B->x-p_A->x;
           line_ref[1] = p_B->y-p_A->y;
           line_ref[2] = p_B->z-p_A->z;

           unit_vector( &line_ref[0], &size);

           line_atom[irestr].ref_x= line_ref[0];
           line_atom[irestr].ref_y= line_ref[1];
           line_atom[irestr].ref_z= line_ref[2];

        } 
   }

/***************************************************************************/
/*** Prepare for MC ********************************************************/
/***************************************************************************/

  entry_no=0;
  if (prob_test >= 1.0)
    {
       printf("In make_a_template found pure MC (prob_test = 100 %%)\n");
       energy_tot_grand = 0.0;
       energy_inter_grand = 0.0;
       num_mc_accepted=0;
    }
  else if (prob_test < 1.0 && action_weights[0] > 0)
    {
       printf("In make_a_template found Build job.\n");
       energy_tot_grand = 0.0;
       energy_inter_grand = 0.0;
       num_mc_accepted=0;
       num_mc_recorded=0;
    }
/***************************************************************************/
/*********** Move template to pore centre of mass **************************/
/***************************************************************************/

    if (centralise_template)
      {
        fprintf(output_fp,"Centralising template\n");

/***** If the pore is empty use the centre of the periodic box as the point to start from ****/
        if (is_empty_pore)
          {
            fprintf(output_fp,"Empty pore so placing template at centre of periodic box\n");

            p_pore_cofm->v[0] = 0.5* (latt_vec[0] +  latt_vec[3] +  latt_vec[6] ); 
            p_pore_cofm->v[1] = 0.5* (latt_vec[1] +  latt_vec[4] +  latt_vec[7] ); 
            p_pore_cofm->v[2] = 0.5* (latt_vec[2] +  latt_vec[5] +  latt_vec[8] ); 
          }
  
        if (want_dock)
          {
         p_demarc=p_guest_demarc;
         p_this_cofm=p_guest_cofm;
         p_this_totmass= p_total_mass_guests;
         for (imol=0; imol < *p_num_seed_mols; imol++)
           {
              if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

//              printf("Before centre of mass of molecule %d: %10.6f %10.6f %10.6f\n", imol,
//                                                     p_this_cofm->v[0],
//                                                     p_this_cofm->v[1],
//                                                     p_this_cofm->v[2]);         
/***********************************************/
/** Place molecule at random position and ******/
/** orientation for dock calcs            ******/
/***********************************************/
/*** First move molecule to origin *************/
/***********************************************/
              vec[0]= - p_this_cofm->v[0];
              vec[1]= - p_this_cofm->v[1];
              vec[2]= - p_this_cofm->v[2];

              move_molecule( guest_ptrs[index], p_demarc->end, &vec[0], -1);

/***********************************************/
/*** Now do random move inside unit cell *******/
/*** to set docking position             *******/
/***********************************************/

              vec[0] = real_random(1) * (latt_vec[0] +  latt_vec[3] +  latt_vec[6] ); 
              vec[1] = real_random(1) * (latt_vec[1] +  latt_vec[4] +  latt_vec[7] ); 
              vec[2] = real_random(1) * (latt_vec[2] +  latt_vec[5] +  latt_vec[8] ); 

              move_molecule( guest_ptrs[index], p_demarc->end, &vec[0], -1);

/*** test DEBUG ***/
             centre_of_mass(&(p_this_cofm->v[0]), p_this_totmass, 
                            guest_ptrs[index], p_demarc->end, -1 );

              printf("Randomised centre of mass of molecule %d: %10.6f %10.6f %10.6f\n", imol,
                                                     p_this_cofm->v[0],
                                                     p_this_cofm->v[1],
                                                     p_this_cofm->v[2]);         

              p_demarc++;
              p_this_cofm++;
              p_this_totmass++;
            }
/***********************************************/
/*** Now do random rotation inside cell  *******/
/*** to set docked orientation           *******/
/***********************************************/

         p_demarc=p_guest_demarc;
         p_this_cofm=p_guest_cofm;
         for (imol=0; imol < *p_num_seed_mols; imol++)
           {
             if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

             dock_orientate(guest_ptrs[index], p_demarc->end,
                                               &(p_this_cofm->v[0]));
             p_demarc++;
             p_this_cofm++;
           }
            
          }
        else
          {
/****** Need to move all the atoms to be near cofm of pore NOT each molecule ****/

            printf("Moving the %d seed molecules to the pore centre of mass\n",*p_num_seed_mols);

/***** find average cofm for all guests ****/
            p_this_cofm=p_guest_cofm;
            vec[0]=0.0;vec[1]=0.0;vec[2]=0.0;

            for (imol=0; imol < *p_num_seed_mols; imol++)
              {
                if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;
                vec[0]+= p_this_cofm->v[0];
                vec[1]+= p_this_cofm->v[1];
                vec[2]+= p_this_cofm->v[2];
                p_this_cofm++;
              }

            rnmol= (double) *p_num_seed_mols;

/**** Make up shift vector ***/
            vec[0]=vec[0]/rnmol; vec[1]=vec[1]/rnmol; vec[2]=vec[2]/rnmol; 
            join_vectors(&vec[0], &(p_pore_cofm->v[0]), &temp_to_pore[0]);

/**** Do shift ***/
            p_demarc=p_guest_demarc;
            for (imol=0; imol < *p_num_seed_mols; imol++)
              {
                if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

                move_molecule( guest_ptrs[index], p_demarc->end, &temp_to_pore[0], -1);

/**** Lets see structure now...DEBUG***/
                for (iii=0; iii<p_demarc->num; iii++)
                  {
                     printf("%s   %10.6f  %10.6f  %10.6f\n", guest_ptrs[index][iii].label,
                                                             guest_ptrs[index][iii].x,
                                                             guest_ptrs[index][iii].y,
                                                             guest_ptrs[index][iii].z);
                  }
/**** Lets see structure now...DEBUG***/
                p_demarc++;
              }
          }
        }
/***********************************************************************/
/************* Apply symmetry if required ******************************/
/***********************************************************************/

       if (symm_set)
          {

/**** First test that the symmetry operations are valid for the pore ****/

             printf("Testing symmetry on the pore. which has %d atoms..\n", num_pore_atoms);

             if (test_pore_symmetry(output_fp, &pore[0], num_pore_atoms,
                                    &symm[0], num_symm_ops) ) printf("Host follows symmetry operations\n");

/**** Need to apply symmetry operations to guest structures now *********/

             do_symmetry_all(guest_ptrs, *p_num_seed_mols, 
                             p_guest_demarc, &symm[0], num_symm_ops);
          }

/******************************************************************************/
/*** Note the energy if required **********************************************/
/******************************************************************************/
     
          if (which_test == NON_BOND_COST)
            {
               for (imol=0; imol < *p_num_seed_mols; imol++) interaction_energy[imol].charges = 0.0;

               just_count=TRUE;
               assemble_intra_lists(guest_ptrs, *p_num_seed_mols, p_guest_demarc, 
                                    num_symm_ops, p_angles_list_ptrs,
                                    p_torsions_list_ptrs, p_vdw_list_ptrs,
                                    p_links_list_ptrs,
                                    &num_angles_counted[0], &num_torsions_counted[0], 
                                    &num_vdws_counted[0], &num_links_counted[0], just_count);

/**** Need to malloc angles_, torsions_ and vdw_ lists and the link_atoms ***/
/**** Note we only need think about the symmetry unique molecules         ***/

               printf("Memory allocation info for intra-molecular potentials: just_counting\n");
               printf("molecule    angles    torsions    links    vdws\n");

               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                   if (num_angles_counted[imol] > 0 )
                                 p_angles_list_ptrs[imol]= (angle_interact_list*) malloc((10+num_angles_counted[imol])*sizeof(angle_interact_list));
                   else
                                 p_angles_list_ptrs[imol]= (angle_interact_list*) malloc(sizeof(angle_interact_list));

                   if (num_torsions_counted[imol] > 0 )
                                 p_torsions_list_ptrs[imol]= (torsion_interact_list*) malloc((10+num_torsions_counted[imol])*sizeof(torsion_interact_list));
                   else
                                 p_torsions_list_ptrs[imol]= (torsion_interact_list*) malloc(sizeof(torsion_interact_list));

                   if (num_vdws_counted[imol] > 0 )
                                 p_vdw_list_ptrs[imol]= (vdw_interact_list*) malloc((10+num_vdws_counted[imol])*sizeof(vdw_interact_list));
                   else
                                 p_vdw_list_ptrs[imol]= (vdw_interact_list*) malloc(sizeof(vdw_interact_list));

                   if (num_links_counted[imol] > 0 )
                                 p_links_list_ptrs[imol]= (links*) malloc((10+num_links_counted[imol])*sizeof(links));
                   else
                                 p_links_list_ptrs[imol]= (links*) malloc(sizeof(links));

                   if (     p_angles_list_ptrs[imol] == NULL
                       || p_torsions_list_ptrs[imol] == NULL
                       ||      p_vdw_list_ptrs[imol] == NULL
                       ||    p_links_list_ptrs[imol] == NULL )
                     {
                       printf("ERROR: Problem assigning memory for intra molecular potentials or links\n");
                       exit(0);
                     }
                   printf(" %d     %d     %d       %d      %d\n", imol, num_angles_counted[imol]+1,
                                                                        num_torsions_counted[imol]+1,
                                                                        num_links_counted[imol]+1,
                                                                        num_vdws_counted[imol]+1);
                                                              
                   num_angles_list[imol]= num_angles_counted[imol];
                   num_torsions_list[imol]= num_torsions_counted[imol];
                   num_links_list[imol]= num_links_counted[imol];
                   num_vdws_list[imol]= num_vdws_counted[imol];
                 }

               just_count=FALSE;
               printf("DEBUG: Now assembling intra lists\n");

               assemble_intra_lists(guest_ptrs, *p_num_seed_mols, p_guest_demarc, 
                                    num_symm_ops, p_angles_list_ptrs,
                                    p_torsions_list_ptrs, p_vdw_list_ptrs,
                                    p_links_list_ptrs,
                                    &num_angles_list[0], &num_torsions_list[0], 
                                    &num_vdws_list[0], &num_links_list[0], just_count);

               printf("DEBUG>> back after filling lists\n");

               printf("Array filling info for intra-molecular potentials:\n");
               printf("molecule    angles    torsions    links    vdws\n");
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                   printf(" %d     %d     %d       %d      %d\n", imol, num_angles_list[imol]+1,
                                                                        num_torsions_list[imol]+1,
                                                                        num_links_list[imol]+1,
                                                                        num_vdws_list[imol]+1);
                 }

               mem_prob=FALSE;
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                    if (num_angles_list[imol] > num_angles_counted[imol]) mem_prob=TRUE;
                    if (num_torsions_list[imol] > num_torsions_counted[imol]) mem_prob=TRUE;
                    if (num_links_list[imol] > num_links_counted[imol]) mem_prob=TRUE;
                    if (num_vdws_list[imol] > num_vdws_counted[imol]) mem_prob=TRUE;
                 }

              if ( mem_prob )
                {
                   printf("ERROR: Poor matching of counted and assigned array sizes will cause memory problems\n");
                   printf("ERROR: Could be that a potential is not available for one of the atoms sets found.\n");
                   printf("ERROR: Update frc file.\n");
                   exit(0);
                }

//          for (imol=0; imol < *p_num_seed_mols; imol++)
//             {
//               printf("Current links list for molecule %d:\n", imol);
//               for (iii=0; iii <= num_links_list[imol]; iii++)
//                  {
//                     ist=p_links_list_ptrs[imol][iii].start;
//                     ied=p_links_list_ptrs[imol][iii].end;
//
//                     printf("%d %d %d\n", iii, ist, ied);
//
//                     printf("%d %d %s %s\n", ist, ied,
//                                            guest_ptrs[imol][ist].label, 
//                                            guest_ptrs[imol][ied].label);
//                  }
//             }

/**** Assemble molmol interaction indices list ******/
               num_molmol_list = (num_seeds_with_symm+1)*(num_seeds_with_symm+2)/2;

//             printf("DEBUG>> mallocing %d space for molmol list\n", num_molmol_list);

               p_molmol_list= (interaction_indices*) malloc(num_molmol_list*sizeof(interaction_indices));
               
//             printf("Have %d molecules so expect %d in the pair interaction list\n",
//                                     num_seeds_with_symm+1, num_molmol_list);

               num_molmol_list= assemble_molmol_list(p_molmol_list, *p_num_seed_mols,
                                                     num_molmol_list);

               if (num_molmol_list<0)
                 {
                   printf("ERROR: Found too few interactions in molmol list\n");
                   exit(0);
                 }
               else printf("Found the expected number of entries in the molmol list\n");

// printf("Molmol list in make_a_template:\n");
// for (iii=0; iii< num_molmol_list; iii++)
//   {
//     printf("ind: %d jnd: %d\n",(p_molmol_list+iii)->ind, (p_molmol_list+iii)->jnd);
//   }

  
//                if (DEBUG)
//                 {
//                    printf("Current links list:\n");
//                    for (imol=0; imol < *p_num_seed_mols; imol++)
//                      {
//                         printf("For molecule %d contains %d entries\n", imol,  num_links_list[imol]);
//                         if (symm_set) index=imol*(num_symm_ops+2);
//                                                       else index=imol;
// 
//                         for (iii=0; iii <= num_links_list[imol]; iii++)
//                           {
//                              ist=p_links_list_ptrs[imol][iii].start;
//                              ied=p_links_list_ptrs[imol][iii].end;
// 
//                              printf("%d %d %s %s\n", ist, ied,
//                                                      guest_ptrs[index][ist].label, 
//                                                      guest_ptrs[index][ied].label);
//                          }
//                      }
//                  }

/*******************************************/
/*** Calculate energy **********************/
/*******************************************/

               printf("DEBUG>> Off to calculate energy....\n");

               printf("\nfirst atom in lists:\n");
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                   printf("%d, %s (elem %s) %10.6f %10.6f %10.6f\n",
                                imol, guest_ptrs[imol][0].label, 
                                      guest_ptrs[imol][0].elem,
                                      guest_ptrs[imol][0].x,
                                      guest_ptrs[imol][0].y,
                                      guest_ptrs[imol][0].z);
                 }

               calculate_energy(&pore[0],num_pore_atoms, 
                                guest_ptrs, *p_num_seed_mols,
                                p_guest_demarc,
                                p_molmol_list, num_molmol_list,
                                p_kvecs, p_kvec2, p_gvec2,
                                num_kvecs, p_cos_sum, p_sin_sum, 
                                p_need_grad, p_grad,
                                have_comb_rules, is_empty_pore);

               printf("DEBUG>> Back from inter-energy calculation, now for intra...\n");

               calculate_intra_energy( guest_ptrs, *p_num_seed_mols,
                                       p_guest_demarc,
                                       p_angles_list_ptrs,
                                       &num_angles_list[0], 
                                       p_torsions_list_ptrs, 
                                       &num_torsions_list[0],
                                       p_vdw_list_ptrs,
                                       &num_vdws_list[0]);

               *p_best_energy= 0.0;
               best_intra=0.0;
               tot_non_bond=0.0; tot_hbond=0.0;
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                   tot_non_bond += interaction_energy[imol].non_bonded;
                   tot_hbond    += interaction_energy[imol].hbond;
                   best_restraint  +=  interaction_energy[imol].restraint;

                   printf("At Start intra energy terms for molecule %d are:\n", imol+1);
                   printf("Bond Stretch : %10.6f\n", intra_energy[imol].stretch);
                   printf("Angle bend   : %10.6f\n", intra_energy[imol].angle  );
                   printf("Torsion      : %10.6f\n", intra_energy[imol].torsion);
                   printf("Van der Waals: %10.6f\n", intra_energy[imol].vdw);

                   printf("\nContribution to inter-molecular non-bond : %10.6f\n\n",
                                                   interaction_energy[imol].non_bonded);

                   best_intra += intra_energy[imol].total;
                 }

               *p_best_energy =  tot_non_bond + tot_hbond + best_restraint;

               printf("At start non_bond= %10.6f hbond=%10.6f restraint =  %10.6f total = %10.6f\n",
                           tot_non_bond, tot_hbond, best_restraint, *p_best_energy);


/*** Don't bother with madelung energies for polymer builder *****/

               if (is_poly)
                 {
                   fprintf(output_fp, "Reporting from Starting values, poly\n");
                 }
               else
                 {
                   fprintf(output_fp, "Reporting from Starting values\n");
                 }
                   
               p_demarc=p_guest_demarc;
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                    if (symm_set) index=imol*(num_symm_ops+2);
                                                  else index=imol;

                    fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

                    print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                                 guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                                 imol, build_number);

                    p_demarc++;
                 }
            }
          else if (which_test == ENERGY_COST)
            {
               calculate_energy(&pore[0], num_pore_atoms, 
                                guest_ptrs, *p_num_seed_mols,
                                p_guest_demarc,
                                p_molmol_list, num_molmol_list,
                                p_kvecs, p_kvec2, p_gvec2,
                                num_kvecs, p_cos_sum, p_sin_sum, 
                                p_need_grad, p_grad,
                                have_comb_rules, is_empty_pore);

               *p_best_energy= 0.0;
               best_intra=0.0;
               tot_non_bond=0.0; tot_hbond=0.0;
               best_restraint = 0.0;
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                   tot_non_bond  += interaction_energy[imol].non_bonded;
                   tot_hbond     += interaction_energy[imol].hbond;
                   tot_charges   += interaction_energy[imol].charges;
                   best_restraint+=  interaction_energy[imol].restraint;
                }
               *p_best_energy=  tot_non_bond + tot_charges
                               +tot_hbond + best_restraint;

               printf("Running full energy cost function\n");
               printf("At start non_bond= %10.6f charge= %10.6f hbond=%10.6f restraint =  %10.6f total = %10.6f\n",
                           tot_non_bond, tot_charges, tot_hbond, best_restraint, *p_best_energy);
            }

/****************************************************************************/
/************** minimize the seed *******************************************/
/************** label template sequentially Starting with fresh numbers *****/
/****************************************************************************/
 
     
     fragments_in_template=1;
     p_demarc=p_guest_demarc;
     for (imol=0; imol < *p_num_seed_mols; imol++)
       {
          if (symm_set) index=imol*(num_symm_ops+2);
                                               else index=imol;

          printf("make_a_template 1: Calling relabel with num_to_label = %d\n", p_demarc->end);
          relabel( guest_ptrs[index], p_demarc->end, &num_elems_used[0], TRUE);
          printf("make_a_template 1: Calling regroup with fragments_in_template = %d\n", fragments_in_template);
          regroup( guest_ptrs[index], p_demarc->end, fragments_in_template);
          printf("Back from regrouping...\n");

          p_demarc++;
       }


     if (initial_minimize_template)
       {
         if ((strcmp(minimizer_name,DISCOVER_MINIMIZER) == 0) || 
             (strcmp(minimizer_name,C2DISCOVER_MINIMIZER) == 0))
           {
             strcpy(discover_root, template_strategy_file);
             if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';
           }

         printf("ERROR: Cannot do minimise in multi-molecule version. Doing minimise_molecule initially\n");
         exit(0);
/*       minimize_molecule(p_template, *p_num_template_atoms, discover_root, */
/*                            p_need_grad, p_grad, is_empty_pore); */

         calculate_energy(&pore[0],num_pore_atoms, 
                          guest_ptrs, *p_num_seed_mols,
                          p_guest_demarc,
                          p_molmol_list, num_molmol_list,
                          p_kvecs, p_kvec2, p_gvec2,
                          num_kvecs, p_cos_sum, p_sin_sum, 
                          p_need_grad, p_grad,
                          have_comb_rules, is_empty_pore);
       }

	       if (animate_flag != 0) 
          {
          if (initial_minimize_template)
            {
              sprintf(&info[0],"ANIMATE: Initial seed Minimize ");
            }
          else
            {
              sprintf(&info[0],"ANIMATE: Initial seed NOT Minimized ");
            }


/***********************************************************************/
/************* Apply symmetry if required ******************************/
/***********************************************************************/

if (symm_set)
   {
       do_symmetry_all(guest_ptrs, *p_num_seed_mols, 
                             p_guest_demarc, &symm[0], num_symm_ops);
   }

          animate( guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                   anim_file_fp, anim_read_fp, anim_show_fp, 
                   &info[0], &this_animation_file[0], &num_anime_frames, 
                   pore);
          }

/**************************************************************************/
/******* Exit if just testing symmetry ************************************/
/**************************************************************************/

        if (test_symm) 
          {
            fclose(anim_file_fp);
            if (anim_pdb_isopen)
              {
                 fclose(anim_pdb_fp);
                 fclose(anim_xyz_fp);
              }
            printf("Exiting [test_symm] in make_a_template\n");
            print_dashes(80,output_fp);
            fprintf(output_fp,"\t\t\t\tSYMMETRY TEST RUN - STOPPING\n");
            print_dashes(80,output_fp);
            fclose(output_fp);
            exit (0);
         }

/**************************************************************************/
/*******ONLY do the initial shaking of flagged to do so********************/
/*******  NEW! on 22/2/96 DWL                          ********************/
/**************************************************************************/

if (initial_minimize_inpore)
		 {

    /**************************************************************************/
    /******** Rigid body move and rotate                   ********************/
    /******** on the seed to try and orientate it sensibly ********************/
    /**************************************************************************/

         calculate_energy(&pore[0],num_pore_atoms, 
                          guest_ptrs, *p_num_seed_mols,
                          p_guest_demarc,
                          p_molmol_list, num_molmol_list,
                          p_kvecs, p_kvec2, p_gvec2,
                          num_kvecs, p_cos_sum, p_sin_sum, 
                          p_need_grad, p_grad,
                          have_comb_rules, is_empty_pore);

       if (is_poly)
        {
          fprintf(output_fp, "Reporting from  initial minimise in pore, poly\n");
          p_demarc=p_guest_demarc;
          for (imol=0; imol < *p_num_seed_mols; imol++)
            {
              if (symm_set) index=imol*(num_symm_ops+2);
                                                  else index=imol;

              fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

              print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                           guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                           imol, build_number);

              p_demarc++;
           }
        }
       else
        {
          fprintf(output_fp, "Reporting from initial minimise in pore\n");

          p_demarc=p_guest_demarc;
          for (imol=0; imol < *p_num_seed_mols; imol++)
            {
               if (symm_set) index=imol*(num_symm_ops+2);
                                                  else index=imol;

               fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

               print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                            guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                            imol, build_number);

               p_demarc++;
            }
        }

       if (*p_num_seed_mols > 1)
        {
          for ( which_mol = 1; which_mol < *p_num_seed_mols + 1; which_mol++)
           {
             idummy= shake_molecule(guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0],
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_shake_attempts, max_shake_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules,
                                    is_empty_pore,
                                    which_mol);

             idummy= rock_molecule (guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0], 
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_rock_attempts, max_rock_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules, 
                                    is_empty_pore, which_mol);
           }
        }
     else
        {
             idummy= shake_molecule(guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0],
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_shake_attempts, max_shake_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules,
                                    is_empty_pore,
                                    which_mol);

             idummy= rock_molecule (guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0], 
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_rock_attempts, max_rock_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules, 
                                    is_empty_pore,
                                    which_mol);
        }

/*********************************************************************/
/*** keep template so that its centre of mass is within min_image ****/
/*** of the origin.              Added Oct 2006 Dave Willock      ****/
/*********************************************************************/
         
//       centre_of_mass(p_c_of_m_temp, &total_mass_temp, p_template,
//                                   *p_num_template_atoms, -1 );

/*** vec is the current vector from the cofm template to that of pore ***/

//       vec[0]= *p_c_of_m_pore     - *p_c_of_m_temp;
//       vec[1]= *(p_c_of_m_pore+1) - *(p_c_of_m_temp+1);
//       vec[2]= *(p_c_of_m_pore+2) - *(p_c_of_m_temp+2);

/*** vec1 is the minimum image version of vec ***************************/

//       vec1[0]=vec[0];
//       vec1[1]=vec[1];
//       vec1[2]=vec[2];

//       min_image(&vec1[0], &vec1[1], &vec1[2]);

/*** Will now shift by the difference between the two *******************/

//       shft_vec[0]= vec[0] - vec1[0];
//       shft_vec[1]= vec[1] - vec1[1];
//       shft_vec[2]= vec[2] - vec1[2];

//       siz = sqrt(  shft_vec[0]*shft_vec[0] 
//                   +shft_vec[1]*shft_vec[1] 
//                   +shft_vec[2]*shft_vec[2]);

//       if (siz > 1E-4)
//          {
//            if (DEBUG) printf("Moving image by %10.6f %10.6f %10.6f\n", shft_vec[0],
//                                                                        shft_vec[1],
//                                                                        shft_vec[2]);

//             move_molecule(p_template, *p_num_template_atoms, &shft_vec[0], -1);
               
//          }
/***********************************************************************/
/************* Apply symmetry if required ******************************/
/***********************************************************************/

       if (symm_set)
          {
              do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
                                    p_guest_demarc, &symm[0], num_symm_ops);
          }


        } /*end if (initial_minimize_inpore)*/

/**************************************************************************/

      if (animate_flag  != 0)  
       {
         if (initial_minimize_inpore)
            {
            sprintf(&info[0],"ANIMATE: after initial rock/shake ");
            }
        else
            {
            sprintf(&info[0],"ANIMATE: NO initial rock/shake");
            }

        animate( guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                 anim_file_fp, anim_read_fp, anim_show_fp, 
                 &info[0], &this_animation_file[0], &num_anime_frames, 
                 pore);
       }

/**************************************************************************/
/****************** Main loop for building the template *******************/
/**************************************************************************/

  have_built= FALSE;
  have_made_ring= FALSE;
  have_minimized = FALSE;
  num_to_rot = 0;
  
//  for (imol=0; imol < *p_num_seed_mols; imol++)
//    {
//      num_old_template_hyds[imol]= *p_num_template_hyds;
//    }

//  sum_old_temp_hyd_weights= *p_sum_temp_hyd_weights;

  iloop =0;
  finished= FALSE;
  need_test= FALSE;
  last_build_number=0;

  while (iloop <= num_modify_attempts && !finished) 
    {
/* reset peek checker if it zeroed out last time ***/
 
      done_move= FALSE;
      animate_this = FALSE;
      loop_time = time(NULL);  /* time at start of this loop */

/***************************************************************************/
/******* Decide what action to take next  **********************************/
/***************************************************************************/

// if (DEBUG) 
//   {
//     printf("SI DEBUG>> Choosing action using pick_frm_wgt_list\n");
//     printf("SI DEBUG>> sending %d %d \n",
//                         sum_action_weights, num_action_weights);
//
//   for (debug_loop=0; debug_loop <= num_action_weights; debug_loop++)
//     {
//       printf("Action %d weight %d\n", 
//                  debug_loop, action_weights[debug_loop]);
//     }
//   }

/****Pick the molecule now, moved out of individual actions, March 2013 Dave Willock ***/

   if (*p_num_seed_mols > 1)
     {
        which_mol = real_random(1) * *p_num_seed_mols;

/*** Assign "which" pointers appropriately ***/
        p_which_demarc = p_guest_demarc+which_mol;

        if (symm_set) which_index=which_mol*(num_symm_ops+2);
                                                  else which_index=which_mol;
        p_which_cofm=p_guest_cofm+which_mol;
        p_which_totmass= p_total_mass_guests+which_mol;
     }
   else
     {
        which_mol = 0;
        p_which_demarc = p_guest_demarc;
        p_which_cofm=p_guest_cofm;
        p_which_totmass= p_total_mass_guests;
     }
   if (DEBUG) printf("Will work with molecule %d\n", which_mol);

   if (which_mol > *p_num_seed_mols)
     {
        printf("ERROR: chosen seed molecule index to shake that is out of range\n");
        printf("       selected : %d when there are only %d molecules.\n", 
                                                       which_mol, *p_num_seed_mols);
        exit(0);
     }

/**** Update Aug 09: make sure action is possible before proceeding, Dave Willock ****/

      can_do=FALSE;
      num_action_failed = 0;
      while(!can_do && num_action_failed < 100)
        {
           action = pick_frm_wgt_list( sum_action_weights, &action_weights[0], 
                                  &weight_of_chosen_action, num_action_weights);

           can_do = TRUE;
           if (action == BUILD_ACTION && (have_built || have_made_ring))
              {
                 can_do = FALSE;
                 num_action_failed++;
              }
           else if (action == RING_MAKER_ACTION  && (have_built || have_made_ring))
              {
                 can_do = FALSE;
                 num_action_failed++;
              }
           else if (action == TWIST_ACTION && num_links_list[which_mol] < 0)
              {
                 can_do = FALSE;
                 num_action_failed++;
              }
        }

      if (num_action_failed >= 1000)
        {
           printf("ERROR: Failed to find a viable action after 1000 attempts....\n");
           printf("       Check action weightings.\n");
           exit(0);
        }

//    if (DEBUG) 
//      {
//        printf("DEBUG>> Trying action %d this time with molecule %d\n", action, which_mol);
//        if (action == 4)
//          {
//             printf("Molecule %d, num_links= %d\n", 
//                                                  which_mol, num_links_list[which_mol]);
//          }
//      }
 
/************************ Update Statistics ********************************/
      statistics[action].tries++;
      statistics[action].tries_after_build += have_built;

      if (action==-1)
        {
           printf("ERROR >> Arrived at action list without a valid action option\n");
           printf("         check that at least some action weights are non-zero\n");
           exit(0);
        }

/***************************************************************************/
/******* select thing to do according to the value of action ***************/
/***************************************************************************/

/**********************************************************************/
/******** action 1 => Join a new fragment to a selected      **********/
/********             guest.                                 **********/
/**********************************************************************/

      if (action == BUILD_ACTION && !have_built && !have_made_ring)
        {
//          printf("\n------------------------\n --Trying to Build on molecule %d-- \n--------------------------\n",
//                           which_mol);
/***************************************************************************/
/***** Remember information for old template for rejection later ***********/
/***************************************************************************/
/***** Need to remember which molecule was used in the latest build ********/

          which_built_mol= which_mol;

/***************************************************************************/
/**** Remember the guest itself ********************************************/
/***************************************************************************/

          p_demarc=p_guest_demarc+which_built_mol;

          p_old_guest= (atom*) malloc((p_demarc->num)*sizeof(atom));
          p_old_demarc= (list_partition*) malloc(sizeof(list_partition));


          p_old_links= (links*) malloc((1+num_links_list[which_built_mol])*sizeof(links));

          num_old_hyds= *(p_num_guest_hyds+which_built_mol);
          p_old_hyd_list= (int*) malloc((1+num_old_hyds)*sizeof(int));
          p_old_hyd_weights= (int*) malloc((1+num_old_hyds)*sizeof(int));

          p_num_this_types = p_num_guest_types+which_built_mol;
          num_old_types = *p_num_this_types;
          p_old_types = (atom_number*) malloc((num_old_types+1)*sizeof(atom_number));

          copy_molecule_data( p_old_guest, &guest_ptrs[which_built_mol][0],
                              p_old_demarc, p_demarc,
                              p_old_links,  &p_links_list_ptrs[which_built_mol][0],
                              &num_old_links, &num_links_list[which_built_mol],
                              p_old_hyd_list, &guest_hyd_list_ptrs[which_built_mol][0],
                              p_old_hyd_weights, &guest_hyd_weights_ptrs[which_built_mol][0], 
                              &num_old_hyds, p_num_guest_hyds+which_built_mol,
                              &sum_old_hyd_weights, p_sum_hyd_weights+which_built_mol,
                              p_old_types, &guest_types_ptrs[which_built_mol][0],
                              &num_old_types, p_num_guest_types+which_built_mol );

          printf("TESTING MOLECULE COPYING ROUTINE===>\n\n");

          printf("DEBUG>>Old demarc: start:%d  num:%d end:%d\n", p_demarc->start,  p_demarc->num,  p_demarc->end);

          p_atom= p_old_guest;
          for (iatom=0; iatom <=p_old_demarc->end; iatom++)
            {
              if (iatom < 10) printf("DEBUG>> Recorded old atom %s %10.6f  %10.6f  %10.6f \n", p_atom->label,
                                                                                               p_atom->x,
                                                                                               p_atom->y,
                                                                                               p_atom->z);
              p_atom++;
            }
          printf("DEBUG>>1\n");

/***************************************************************************/
/*** Link list is now more than inter-fragment so need num_links saved too */
/*** added Aug. 09 Dave Willock                                            */
/*** Now need to malloc space for the old_links list.                      */
/***************************************************************************/

          printf("Recorded %d entries in old links list:\n", num_old_links+1);
          p_this_link=p_old_links; 
          for (ilink=0; ilink <= num_old_links; ilink++)            
             {                                                   
                printf("Start: %d End: %d\n", p_this_link->start, p_this_link->end);
                p_this_link++;
             }                                                     

/******************************************************************************/
/******** Remember old hydrogen                         ***********************/
/******** and its weights for rejection later           ***********************/
/******************************************************************************/

          printf("Recorded %d entries in old hyds list:\n", num_old_hyds+1);

          p_int= p_old_hyd_list;                       
          p_weight= p_old_hyd_weights;                    

          for (ihyd=0; ihyd <= num_old_hyds; ihyd++) 
             {                                        
                printf("%d: Hyd index: %d label %s (elem %s) weight %d\n", ihyd, *p_int,
                                                                            guest_ptrs[which_built_mol][*p_int].label, 
                                                                            guest_ptrs[which_built_mol][*p_int].elem, 
                                                                           *p_weight); 
                p_int++;                             
                p_weight++;                           
             }                                         

          printf("Recorded sum H weights as: %d\n", sum_old_hyd_weights);

/***************************************************************************/
/*** Remember old types list ***********************************************/
/***************************************************************************/

        p_type= p_old_types;
        for (itype=0; itype<num_old_types; itype++)
           {
              printf("copying: %d %s %d\n", itype, guest_types_ptrs[which_built_mol][itype].atom_type,
                                                   guest_types_ptrs[which_built_mol][itype].num);
              p_type++; 
           }

        printf("Recorded %d old types:\n", num_old_types);

        p_type= p_old_types;
        for (itype=0; itype<num_old_types; itype++)
          {
             printf("%d %s %d\n", itype, p_type->atom_type, p_type->num);
             p_type++;
          }

/***************************************************************************/
/******* Build new template section from fragment library                  */
/******* Even though this is a multi-molecule version build only need      */
/******* see one molecule at a time. Now need to split build into a        */
/******* fragment decision section and the join_frag_to_template function  */
/******* This allows a suitalbe realloc call to be made for space to be    */
/******* allocated to the new part of the array.                           */
/***************************************************************************/

      printf("DEBUG>> In make_a_template before choosing fragment sum_hyd_weights %d old %d\n",  
                                                                         *(p_sum_hyd_weights+which_built_mol),
                                                                         sum_old_hyd_weights );

        which_frag= choose_fragment(&guest_ptrs[which_built_mol][0], p_demarc->end, p_frag_weights,
                                    guest_hyd_list_ptrs[which_built_mol], *(p_num_guest_hyds+which_built_mol),
                                    p_frag_hyd_list, p_frag_hyd_partition,
                                    frag_lib, &guest_types_ptrs[which_built_mol][0], 
                                    *p_num_this_types, &num_new_types, 
                                    p_frag_types, p_frag_types_list, *p_have_AB);

        printf("DEBUG>> In make_a_template Chosen fragment %d for a build....\n", which_frag);
        printf("DEBUG>> 1. sum_hyd_weights now %d old %d\n", *(p_sum_hyd_weights+which_built_mol), sum_old_hyd_weights );

/*** realloc space for new guest atoms ****/
   
         num_new = p_demarc->num + number_of_members[which_frag]+1;
         printf("DEBUG>> Need to allocate space for %d atoms to do build\n", num_new); 

         guest_ptrs[which_built_mol]= (atom*) realloc( guest_ptrs[which_built_mol], 2*num_new * sizeof(atom));

/**** Copy fragment into new part of the array ****/

         p_atom= &guest_ptrs[which_built_mol][p_demarc->num];

         for (iatom=0; iatom <= number_of_members[which_frag]; iatom++)
           {
             *p_atom = frag_lib[ member_start[which_frag] + iatom ];
             p_atom++;
           }

/**** update the demarcation structure, note that the start will remain the same ****/
 
         p_demarc->num = num_new; p_demarc->end=num_new-1;
           
       p_atom= &guest_ptrs[which_built_mol][0];
       for (iatom=0; iatom < num_new; iatom++)
         {
            printf("%d %s %10.6f %10.6f %10.6f %s\n", iatom,
                                                      p_atom->label,
                                                      p_atom->x,
                                                      p_atom->y,
                                                      p_atom->z,
                                                      p_atom->elem);
            p_atom++;
         }

//       printf("So now before connecting have %d atoms, max index %d\n", p_demarc->num, p_demarc->end); 

/*** realloc space for new guest atoms ****/

//       printf("Will need %d types for new guest\n",num_new_types);

         guest_types_ptrs[which_built_mol] =
              (atom_number*) realloc( guest_types_ptrs[which_built_mol], (num_new_types+1)* sizeof(atom_number) );

//       printf("Made Space....\n");
/*** Update types array ****/

         copy_types( &guest_types_ptrs[which_built_mol][0], p_num_this_types,
                     p_old_types, num_old_types);

//     printf("Copied types\n");
//   p_type= p_old_types;
//    for (itype=0; itype<*p_num_this_types; itype++)
//      {
//         printf("%d %s %d\n", itype, guest_types_ptrs[which_built_mol][itype].atom_type,
//                                     guest_types_ptrs[which_built_mol][itype].num );
//
//         if (itype < num_old_types) 
//           {
//              printf("Old: %d %s %d\n", itype, p_type->atom_type,
//                                                       p_type->num );
//              p_type++;
//           }
//      }

          p_type = p_frag_types+(p_frag_types_list+which_frag)->start;

//      printf("DEBUG>> 2. sum_hyd_weights now %d old %d\n", *(p_sum_hyd_weights+which_built_mol), sum_old_hyd_weights );

//      printf("Adding fragment types...\n");

          add_types( &guest_types_ptrs[which_built_mol][0], p_num_this_types,
                     p_type, (p_frag_types_list+which_frag)->num);

//    printf("Updated types array has %d members:\n", *p_num_this_types);
//    for (itype=0; itype<*p_num_this_types; itype++)
//      {
//         printf("%d %s %d\n", itype, guest_types_ptrs[which_built_mol][itype].atom_type,
//                                     guest_types_ptrs[which_built_mol][itype].num );
//      }
 
/**** Deal with hydrogen indexing lists *****/

         p_this_frag_hyd_partition= p_frag_hyd_partition + which_frag;
         start_H_list = p_this_frag_hyd_partition->start;
         p_this_frag_H_list = p_frag_hyd_list + start_H_list;

//     printf("Sending fragment %d to connect_new_fragment\n", which_frag);
//     printf("This has %d hydrogens:\n", p_this_frag_hyd_partition->num);
//     for (ihyd=0; ihyd<=p_this_frag_hyd_partition->num; ihyd++)
//       {
//         iatom = p_old_demarc->num+*(p_this_frag_H_list+ihyd);
//         printf("%d -> %d index %d %s (elem %s)\n", ihyd, *(p_this_frag_H_list+ihyd), iatom,
//                                          guest_ptrs[which_built_mol][iatom].label,  
//                                          guest_ptrs[which_built_mol][iatom].elem );
//       }

/** Reallocate hydrogen list and weight arrays to deal with new guest *******/

         guest_hyd_list_ptrs[which_built_mol] = (int*) realloc( guest_hyd_list_ptrs[which_built_mol],
                                                              ( *(p_num_guest_hyds+which_built_mol) + p_this_frag_hyd_partition->num + 200)
                                                                                                                            * sizeof(int));
 
         guest_hyd_weights_ptrs[which_built_mol] = (int*) realloc(  guest_hyd_weights_ptrs[which_built_mol],
                                                                 ( *(p_num_guest_hyds+which_built_mol) + p_this_frag_hyd_partition->num + 200)
                                                                                                                               * sizeof(int));
 
//
// Need to sort out weights for fragment H atoms.
//
//  printf("DEBUG>> 3. sum_hyd_weights now %d old %d\n", *(p_sum_hyd_weights+which_built_mol), sum_old_hyd_weights );
//  printf("DEBUG>> 3. guest has  %d hydrogens available\n", *(p_num_guest_hyds+which_built_mol) );

         connect_new_fragment(&guest_ptrs[which_built_mol][0], p_old_demarc, p_demarc,
                              &place_for_one, &place_for_two,
                              &guest_hyd_list_ptrs[which_built_mol][0], (p_num_guest_hyds+which_built_mol),
                              p_this_frag_H_list, p_this_frag_hyd_partition,
                              &guest_types_ptrs[which_built_mol][0], p_num_this_types,
                              &guest_hyd_weights_ptrs[which_built_mol][0], p_sum_hyd_weights+which_built_mol,
                              increment_hydrogen_weight, p_have_AB);
                              
//   printf("So now after connecting have %d atoms, max index %d\n", p_demarc->num, p_demarc->end); 
//   printf("Updated types array after new connection has %d members:\n", *p_num_this_types);

//   for (itype=0; itype<*p_num_this_types; itype++)
//     {
//        printf("%d %s %d\n", itype, guest_types_ptrs[which_built_mol][itype].atom_type,
//                                    guest_types_ptrs[which_built_mol][itype].num );
//     }

//   printf("Sum H-weights for the %d hyds after connecting is %d old %d\n", 
//                                             *(p_num_guest_hyds+which_built_mol),
//                                             *(p_sum_hyd_weights+which_built_mol), 
//                                             sum_old_hyd_weights);

/*** Update intra-molecular interactions and links list for just this molecule *****/

         assign_stretch(guest_ptrs[which_built_mol], p_demarc->end);

         just_count=TRUE;
         assemble_intra_lists(&guest_ptrs[which_built_mol], 1, p_demarc, 
                              num_symm_ops, &p_angles_list_ptrs[which_built_mol],
                              &p_torsions_list_ptrs[which_built_mol], &p_vdw_list_ptrs[which_built_mol],
                              &p_links_list_ptrs[which_built_mol],
                              &num_angles_counted[which_built_mol], &num_torsions_counted[which_built_mol], 
                              &num_vdws_counted[which_built_mol], &num_links_counted[which_built_mol], just_count);

/**** Need to malloc angles_, torsions_ and vdw_ lists and the link_atoms ***/
/**** Note we only need think about the symmetry unique molecules         ***/

//     printf("Memory re-allocation info for intra-molecular potentials: just_counting for the newly built molecule %d\n", 
//                                                                                                       which_built_mol);
//     printf("molecule    angles    torsions    links    vdws\n");

               if (num_angles_counted[which_built_mol] > 0 )
                  p_angles_list_ptrs[which_built_mol]= 
                       (angle_interact_list*) realloc( p_angles_list_ptrs[which_built_mol],
                                                      (10+num_angles_counted[which_built_mol])*sizeof(angle_interact_list));
               else
                   p_angles_list_ptrs[which_built_mol]= (angle_interact_list*) realloc( p_angles_list_ptrs[which_built_mol],
                                                                                        sizeof(angle_interact_list));

               if (num_torsions_counted[which_built_mol] > 0 )
                  p_torsions_list_ptrs[which_built_mol]= 
                      (torsion_interact_list*) realloc( p_torsions_list_ptrs[which_built_mol],
                                                (10+num_torsions_counted[which_built_mol])*sizeof(torsion_interact_list));
               else
                  p_torsions_list_ptrs[which_built_mol]= 
                      (torsion_interact_list*) realloc( p_torsions_list_ptrs[which_built_mol], sizeof(torsion_interact_list));

               if (num_vdws_counted[which_built_mol] > 0 )
                 p_vdw_list_ptrs[which_built_mol]= 
                       (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_built_mol],
                                                    (10+num_vdws_counted[which_built_mol])*sizeof(vdw_interact_list));
               else
                 p_vdw_list_ptrs[which_built_mol]= 
                       (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_built_mol], sizeof(vdw_interact_list));

               if (num_links_counted[which_built_mol] > 0 )
                             p_links_list_ptrs[which_built_mol]= 
                                  (links*) realloc(p_links_list_ptrs[which_built_mol], 
                                                         (10+num_links_counted[which_built_mol])*sizeof(links));
               else
                             p_links_list_ptrs[which_built_mol]= 
                                  (links*) realloc(p_links_list_ptrs[which_built_mol], sizeof(links));

               if (     p_angles_list_ptrs[which_built_mol] == NULL
                   || p_torsions_list_ptrs[which_built_mol] == NULL
                   ||      p_vdw_list_ptrs[which_built_mol] == NULL
                   ||    p_links_list_ptrs[which_built_mol] == NULL )
                 {
                   printf("ERROR: Problem assigning memory for intra molecular potentials or links\n");
                   exit(0);
                 }
//           printf(" %d     %d     %d       %d      %d\n", which_built_mol, num_angles_counted[which_built_mol]+1,
//                                                                num_torsions_counted[which_built_mol]+1,
//                                                                num_links_counted[which_built_mol]+1,
//                                                                num_vdws_counted[which_built_mol]+1);
                                                          
               num_angles_list[which_built_mol]= num_angles_counted[which_built_mol];
               num_torsions_list[which_built_mol]= num_torsions_counted[which_built_mol];
               num_links_list[which_built_mol]= num_links_counted[which_built_mol];
               num_vdws_list[which_built_mol]= num_vdws_counted[which_built_mol];

               just_count=FALSE;
//           printf("DEBUG: Now assembling intra lists for the just built molecule %d\n",which_built_mol);

               assemble_intra_lists(&guest_ptrs[which_built_mol], 1, p_demarc, 
                                    num_symm_ops, &p_angles_list_ptrs[which_built_mol],
                                    &p_torsions_list_ptrs[which_built_mol], &p_vdw_list_ptrs[which_built_mol],
                                    &p_links_list_ptrs[which_built_mol],
                                    &num_angles_list[which_built_mol], &num_torsions_list[which_built_mol], 
                                    &num_vdws_list[which_built_mol], &num_links_list[which_built_mol], just_count);

//           printf("DEBUG>> back after filling lists\n");

//           printf("Array filling info for intra-molecular potentials:\n");
//           printf("molecule    angles    torsions    links    vdws\n");
//           printf(" %d     %d     %d       %d      %d\n", which_built_mol, num_angles_list[which_built_mol]+1,
//                                                                num_torsions_list[which_built_mol]+1,
//                                                                num_links_list[which_built_mol]+1,
//                                                                num_vdws_list[which_built_mol]+1);

               mem_prob=FALSE;
               if (num_angles_list[which_built_mol] > num_angles_counted[which_built_mol]) mem_prob=TRUE;
               if (num_torsions_list[which_built_mol] > num_torsions_counted[which_built_mol]) mem_prob=TRUE;
               if (num_links_list[which_built_mol] > num_links_counted[which_built_mol]) mem_prob=TRUE;
               if (num_vdws_list[which_built_mol] > num_vdws_counted[which_built_mol]) mem_prob=TRUE;

              if ( mem_prob )
                {
                   printf("ERROR: Poor matching of counted and assigned array sizes will cause memory problems\n");
                   printf("ERROR: Could be that a potential is not available for one of the atoms sets found.\n");
                   printf("ERROR: Update frc file.\n");
                   exit(0);
                }

//           printf("Current links list for molecule %d:\n", imol);
//           for (iii=0; iii <= num_links_list[which_built_mol]; iii++)
//              {
//                 ist=p_links_list_ptrs[which_built_mol][iii].start;
//                 ied=p_links_list_ptrs[which_built_mol][iii].end;
//
//                 printf("%d %d %s %s\n", ist, ied,
//                                        guest_ptrs[which_built_mol][ist].label, 
//                                        guest_ptrs[which_built_mol][ied].label);
//              }

/**** Assemble molmol interaction indices list ******/
               num_molmol_list = (num_seeds_with_symm+1)*(num_seeds_with_symm+2)/2;

//           printf("DEBUG>> mallocing %d space for molmol list\n", num_molmol_list);

               p_molmol_list= (interaction_indices*) malloc(num_molmol_list*sizeof(interaction_indices));
               
//           printf("Have %d molecules so expect %d in the pair interaction list\n", num_seeds_with_symm, num_molmol_list);

               num_molmol_list= assemble_molmol_list(p_molmol_list, *p_num_seed_mols,
                                                     num_molmol_list);

               if (num_molmol_list<0)
                 {
                   printf("ERROR: Found too few interactions in molmol list\n");
                   exit(0);
                 }
//           else printf("Found the expected number of entries in the molmol list\n");

//
// Now write a file to see what we have built
//
//
//     printf("....................Writing a debug pdb file........\n");
//     sprintf(debug_filename,"debug.pdb");
//
//     if (open_file(&debug_fp, debug_filename, "w") == EXIT_FAILURE)
//       {
//          printf("ERROR>> Cannot open pdb file %s for writing\n", debug_filename);
//          exit(0);
//       }
//
//     write_pdb( debug_fp, guest_ptrs, *p_num_seed_mols, p_guest_demarc,
//                &abc[0], &super[0], &latt_vec[0]);
//     
//     fclose(debug_fp);
//
//     printf("....................Debug pdb file written........\n");

         have_built= TRUE;

/*************************************************************************/
/**** Calculate the energy immediately after the new build ***************/
/*************************************************************************/

         idummy = cost_function(guest_ptrs, *p_num_guests, p_guest_demarc,
                                pore, &box_limits[0],
                                p_molmol_list, num_molmol_list,
                                which_test, &cc_tot,
                                FALSE, &dummy_int,
                                p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                p_cos_sum, p_sin_sum,
                                p_need_grad, p_grad,
                                have_comb_rules,
                                is_empty_pore, 
                                p_angles_list_ptrs, &num_angles_list[0],
                                p_torsions_list_ptrs,&num_torsions_list[0], 
                                p_vdw_list_ptrs,&num_vdws_list[0],
                                TRUE);

//        printf("DEBUG>> cost_function calculated\n");
       }
 
/**********************************************************************/
/******** action 2 => Rotate latest addition to molecule     **********/
/********             accept if the cost function is         **********/
/********             improved                               **********/
/**********************************************************************/

     else if (action == ROTATE_LAST_FRAG_ACTION && (have_built && !have_made_ring))
       {


          printf("ERROR: Rotate last fragment action attempted when new multi-molecule code\n");
          printf("ERROR: is not ready to do this.\n");		 		   
          exit(0);

//        if (DEBUG) printf("DEBUG>> Trying Rotating last segment\n");
// 
//        rot_section( p_template, *p_num_template_atoms, &pore[0], 
//                     &box_limits[0], 
//                     interaction_indices *p_molmol_list, 
//                     int num_molmol_list,
//                     &symm[0],
//                     which_test, &link_atoms[num_links], 
//                     num_to_rot, max_rock_step, 
//                     p_kvecs, p_kvec2, p_gvec2, num_kvecs,
//                     p_cos_sum, p_sin_sum, p_need_grad, p_grad,
//                     have_comb_rules, 
//                     is_empty_pore, 
//                     p_angles_list_ptrs, &num_angles_list[0],
//                     p_torsions_list_ptrs,&num_torsions_list[0], 
//                     p_vdw_list_ptrs,&num_vdws_list[0],
//                     which_mol);

       }

/**********************************************************************/
/******** action 3 => Move whole molecule randomly and       **********/
/********             accept if the cost function is         **********/
/********             improved                               **********/
/********             Added multiple molecule seed options   **********/
/********             March 07 Dave Willock                  **********/
/**********************************************************************/

     else if (action == SHAKE_ACTION)
       {
//          printf("\n-----------------------------------------------------");
//          printf("\niloop= %d Attempting SHAKE for molecule %d max_shake_step = %10.6f\n", iloop, which_mol, max_shake_step);
//          if (DEBUG)
//           {
//             printf("Distance matrix before shake:\n");
//             if (symm_set)
//              {
//                 printf("symmetry set case, %d operations\n", num_symm_ops);
//                 for (isymm= 0; isymm <= num_symm_ops; isymm++)   
//                            print_dist_mat( guest_ptrs[which_index+isymm], 
//                                            p_which_demarc->end, TRUE, stdout);
//              }
//             else
//              {
//                print_dist_mat( guest_ptrs[which_index], p_which_demarc->end, TRUE, stdout);
//              }
//           }

         done_move = shake_molecule(guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0],
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_shake_attempts, max_shake_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules,
                                    is_empty_pore,
                                    which_mol);

         if (DEBUG)
          {
            if (done_move) fprintf(output_fp, "successful shake\n");
                                         else fprintf(output_fp,"failed shake\n");
          }

/*********************************************************************/
/*** keep template so that its centre of mass is within min_image ****/
/*** of the origin.              Added Oct 2006 Dave Willock      ****/
/*********************************************************************/
         
        if (done_move)
          {
            centre_of_mass(&(p_which_cofm->v[0]), p_which_totmass, 
                            guest_ptrs[which_index], p_which_demarc->end, -1 );

/*** vec is the current vector from the cofm template to that of pore ***/

            vec[0]= p_pore_cofm->v[0] - p_which_cofm->v[0];
            vec[1]= p_pore_cofm->v[1] - p_which_cofm->v[1];
            vec[2]= p_pore_cofm->v[2] - p_which_cofm->v[2];

/*** vec1 is the minimum image version of vec ***************************/

            vec1[0]=vec[0];
            vec1[1]=vec[1];
            vec1[2]=vec[2];

            min_image(&vec1[0], &vec1[1], &vec1[2]);

/*** Will now shift by the difference between the two *******************/

            shft_vec[0]= vec[0] - vec1[0];
            shft_vec[1]= vec[1] - vec1[1];
            shft_vec[2]= vec[2] - vec1[2];

            siz = sqrt(  shft_vec[0]*shft_vec[0] 
                        +shft_vec[1]*shft_vec[1] 
                        +shft_vec[2]*shft_vec[2]);

            if (siz > 1E-4)
               {
//                if (DEBUG)
//                  {
//                     printf("Moving image by %10.6f %10.6f %10.6f\n", shft_vec[0],
//                                                                         shft_vec[1],
//                                                                         shft_vec[2]);
//                     printf("Distance matrix before cofm move:\n");
//                     print_dist_mat( guest_ptrs[which_index], p_which_demarc->end, TRUE, stdout);
//                  }

                  move_molecule(guest_ptrs[which_index], p_which_demarc->end, &shft_vec[0], -1);
             
//                if (DEBUG)
//                  {
//                     printf("Distance matrix after cofm move:\n");
//                     print_dist_mat( guest_ptrs[which_index], p_which_demarc->end, TRUE, stdout);
//                  }
               }


/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/

             if (symm_set)
               {
                  do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
                                        p_guest_demarc, &symm[0], num_symm_ops);
               }
           }

       }


/**********************************************************************/
/******** action 4 => Rotate whole molecule randomly and     **********/
/********             accept if the cost function is         **********/
/********             improved                               **********/
/**********************************************************************/

     else if (action == ROCK_ACTION)
       {
//          printf("\n-----------------------------------------------------");
//          printf("\niloop= %d Attempting ROCK for molecule %d\n", iloop, which_mol);

         done_move = rock_molecule (guest_ptrs, *p_num_seed_mols,
                                    p_guest_demarc,
                                    &pore[0], &box_limits[0], 
                                    p_molmol_list, num_molmol_list,
                                    &symm[0], which_test,
                                    num_rock_attempts, max_rock_step,
                                    p_kvecs, p_kvec2, p_gvec2,
                                    num_kvecs, p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules, 
                                    is_empty_pore,
                                    which_mol);

         if (DEBUG)
          {
            if (done_move) fprintf(output_fp, "successful rock\n");
                                         else fprintf(output_fp,"failed rock\n");
          }

        if (DEBUG && done_move)
          {
            printf("DEBUG>> We have rocked Distances after rocking \n");

            print_dist_mat( guest_ptrs[which_index], p_which_demarc->end, TRUE, stdout);

           if (have_built)
            { 
               printf("\n DEBUG animating a rock after build  as check\n");
            }
          }
       else if (DEBUG)
         {
           printf("Failed rock\n");
         }

/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/

          if (symm_set)
            {
               do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
                                 p_guest_demarc, &symm[0], num_symm_ops);
            }
       }

/**********************************************************************/
/******** action 5 => Rotate a random section of the template**********/
/********             choosen from the links list at random  **********/
/**********************************************************************/

     else if (action == TWIST_ACTION && num_links_list[which_mol] >= 0)
       {
//          printf("\n-----------------------------------------------------");
//          printf("\niloop= %d Attempting twist for molecule %d , num_attempts %d\n", iloop, which_mol, num_rock_attempts);
/************ choose a link to rotate round *********************/

          ilink=  real_random(1)*(num_links_list[which_mol]+1); 
//          printf("Number of possible links: %d choosing %d\n", num_links_list[which_mol]+1, ilink);

          if (ilink > num_links_list[which_mol]) 
            {
               printf("ERROR: During twist action ilink = %d but num_links is %d\n", 
                                                                   ilink, num_links_list[which_mol]);
               exit(0);
            }

         if (which_mol > *p_num_seed_mols)
           {
             printf("ERROR: chosen seed molecule index to twist that is out of range\n");
             printf("       selected : %d when there are only %d molecules.\n",
                                                               which_mol, *p_num_seed_mols);
             exit(0);
           }

/**------------------------------------------**/
/**-- Need to re-do vdw list after a twist --**/
/**-- to avoid steric clashes ---------------**/
/**---Dave Willock, Feb. 2014 ---------------**/
/**------------------------------------------**/
//

        allow=TRUE; num_rots=0;
        p_flag_list=malloc(p_which_demarc->num*sizeof(int));

        set_twist_axis_flags(guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                             p_flag_list, &(p_links_list_ptrs[which_mol][ilink]),
                             &origin[0], &axis[0], which_mol);

        while (allow && num_rots < num_rock_attempts)
          {
              num_rots++;
              theta = real_random(1)*max_rock_step;
              isign =  real_random(1) < 0.5;
              if (isign) theta= -1.0*theta;

//              printf("Set twist angle %10.6f degrees\n", theta*RAD_TO_DEG);
//              printf(" Before anything the twist, vdw energy now: %10.6f\n", intra_energy[which_mol].vdw );

              do_the_twist(guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                           p_flag_list, &symm[0], &(p_links_list_ptrs[which_mol][ilink]), 
                           &origin[0], &axis[0], &theta, which_mol);

              just_count=TRUE;
              num_vdws_counted[which_mol]=assemble_vdw_list(guest_ptrs[which_mol], p_which_demarc->end, 
                                                            p_vdw_list_ptrs[which_mol], just_count);

//              printf("After twist have %d vdw interactions counted\n", num_vdws_counted[which_mol]);

              if (num_vdws_counted[which_mol] > 0 )
                 p_vdw_list_ptrs[which_mol]=
                      (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_mol],
                                                  (10+num_vdws_counted[which_mol])*sizeof(vdw_interact_list));
               else
                 p_vdw_list_ptrs[which_mol]=
                      (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_mol], sizeof(vdw_interact_list));

              just_count=FALSE; mem_prob=FALSE;
              num_vdws_list[which_mol]=assemble_vdw_list(guest_ptrs[which_mol], p_which_demarc->end,
                                                         p_vdw_list_ptrs[which_mol], just_count);

              if (num_vdws_list[which_mol] > num_vdws_counted[which_mol]) mem_prob=TRUE;

              if ( mem_prob )
                {
                   printf("ERROR: Poor matching of counted and assigned array sizes will cause memory problems\n");
                   printf("ERROR: When assigning vdw intra-potentials during a twist action.                 \n");
                   printf("ERROR: Could be that a potential is not available for one of the atoms sets found.\n");
                   printf("ERROR: Update frc file.\n");
                   exit(0);
                }

//              printf("In twist option have now assigned %d vdw interactions\n", num_vdws_list[which_mol]);

/**** Now need to get energy for new configuration and test if we should accept it. ****/
/**** If not need to undo the rotation and re-set the vdw_list for the molecule     ****/
//        printf("Before testing the twist, vdw energy: %10.6f\n", intra_energy[which_mol].vdw );
//        printf("Before testing the twist, non_bond: %10.6f intra: %10.6f\n", interaction_energy[which_mol].non_bonded,
//                                                                             intra_energy[which_mol].total );

        allow= test_the_twist(pore, guest_ptrs, *p_num_guests, p_guest_demarc,
                              &box_limits[0], p_molmol_list, num_molmol_list,
                              which_test,
                              p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                              p_cos_sum, p_sin_sum,
                              p_need_grad, p_grad,
                              have_comb_rules,
                              is_empty_pore,
                              p_angles_list_ptrs, &num_angles_list[0],
                              p_torsions_list_ptrs,&num_torsions_list[0],
                              p_vdw_list_ptrs,&num_vdws_list[0],
                              which_mol);

//         printf(" tested the twist, vdw energy now: %10.6f\n", intra_energy[which_mol].vdw );
//         if (allow) printf("twist test says OK\n"); else printf("twist test says forget this one...\n");

/*** If we allow once we should flag that the move has been done, Dave Willock Feb. 2016 ***/
         if (allow) done_move=TRUE;

/***** Act on twist result ********/

         if (!allow)
           {
//               printf("Reverting to previous structure before the twist\n");

/**** Fix structure back to former configuration ****/
/**** But note that this only undoes the final step in a multi-step twist event so all others may be OK ***/
             theta= -theta;

             do_the_twist(guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                          p_flag_list, &symm[0], &(p_links_list_ptrs[which_mol][ilink]), 
                          &origin[0], &axis[0], &theta, which_mol);

/**** re-do the vdw_list ***/
             just_count=TRUE;
             num_vdws_counted[which_mol]=assemble_vdw_list(guest_ptrs[which_mol], p_which_demarc->end, 
                                                           p_vdw_list_ptrs[which_mol], just_count);

//             printf("After twist reject have %d vdw interactions counted\n", num_vdws_counted[which_mol]);

             if (num_vdws_counted[which_mol] > 0 )
                 p_vdw_list_ptrs[which_mol]=
                      (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_mol],
                                                  (10+num_vdws_counted[which_mol])*sizeof(vdw_interact_list));
             else
                 p_vdw_list_ptrs[which_mol]=
                      (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_mol], sizeof(vdw_interact_list));

             just_count=FALSE; mem_prob=FALSE;
             num_vdws_list[which_mol]=assemble_vdw_list(guest_ptrs[which_mol], p_which_demarc->end,
                                                         p_vdw_list_ptrs[which_mol], just_count);

             if (num_vdws_list[which_mol] > num_vdws_counted[which_mol]) mem_prob=TRUE;

             if ( mem_prob )
                {
                   printf("ERROR: Poor matching of counted and assigned array sizes will cause memory problems\n");
                   printf("ERROR: When assigning vdw intra-potentials after a rejection of a twist action.    \n");
                   printf("ERROR: Could be that a potential is not available for one of the atoms sets found.\n");
                   printf("ERROR: Update frc file.\n");
                   exit(0);
                }

//              printf("After twist rejection have now assigned %d vdw interactions\n", num_vdws_list[which_mol]);

/**** Now need to get energy for new configuration and test if we should accept it. ****/
/**** If not need to undo the rotation and re-set the vdw_list for the molecule     ****/
//          printf("Before testing the twist, vdw energy: %10.6f\n", intra_energy[which_mol].vdw );

/*** Use the test to reset the energies ***/
          allow= test_the_twist(pore, guest_ptrs, *p_num_guests, p_guest_demarc,
                                &box_limits[0], p_molmol_list, num_molmol_list,
                                which_test,
                                p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                p_cos_sum, p_sin_sum,
                                p_need_grad, p_grad,
                                have_comb_rules,
                                is_empty_pore,
                                p_angles_list_ptrs, &num_angles_list[0],
                                p_torsions_list_ptrs,&num_torsions_list[0],
                                p_vdw_list_ptrs,&num_vdws_list[0],
                                which_mol);

//          printf(" tested the twist restored, vdw energy now: %10.6f\n", intra_energy[which_mol].vdw );
//          if (allow) printf("twist test says OK\n"); else printf("twist test says forget this one...\n");

           }
//         else
//           {
//             printf("Keeping the new configuration after test\n");
//           }

        }
      free(p_flag_list);

/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/

//        if (symm_set)
//          {
//             do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
//                                   p_guest_demarc, &symm[0], num_symm_ops);
//          }
//         printf("\n");
       }

/**********************************************************************/
/******** action 6 => Attempt to find rings to form in the   **********/
/********             template                               **********/
/**********************************************************************/

         else if (action == RING_MAKER_ACTION  && !have_built && !have_made_ring)
           {
             if (DEBUG) printf("DEBUG>> Trying to make a ring\n");

             printf("ERROR: Attempt to use ring_maker action when this is not available\n");
             printf("ERROR: in multi-molecule version.\n");
             exit(0);
             
//           have_made_ring= ring_maker(p_template, p_num_template_atoms, 
//                                      &forbidden_bond[0], p_template_hyd_list,
//                                      p_num_template_hyds, p_temp_hyd_weights,
//                                      p_sum_temp_hyd_weights);  
//
//           if ( have_made_ring )
//             {
//                printf("Ring_maker sucess!! returned %d atoms:\n", *p_num_template_atoms);
//
//                relabel(p_template, *p_num_template_atoms, &num_elems_used[0], TRUE);

/****** Minimize the template w. r. to itself (no pore) *************/

//           if ((strcmp(minimizer_name,DISCOVER_MINIMIZER) == 0) || 
//                      (strcmp(minimizer_name,C2DISCOVER_MINIMIZER) == 0))
//                {
//                 strcpy(discover_root, template_strategy_file);
//                  if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';
//                }

//                minimize_molecule(p_template, *p_num_template_atoms,
//                                  discover_root, p_need_grad, p_grad, is_empty_pore);

//                calculate_energy(&pore[0],num_pore_atoms, 
//                                 guest_ptrs, *p_num_seed_mols,
//                                 p_guest_demarc,
//                                 p_molmol_list, num_molmol_list,
//                                 p_kvecs, p_kvec2, p_gvec2,
//                                 num_kvecs, p_cos_sum, p_sin_sum, 
//                                 p_need_grad, p_grad,
//                                 symm_set, num_symm_ops,
//                                 have_comb_rules, is_empty_pore);

/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/
//               if (symm_set)
//                  {
//                      do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
//                                            p_guest_demarc, &symm[0], num_symm_ops);
//                  }

/****************************************************************************/
/**** Re-do intra molecular interaction lists *******************************/
/****************************************************************************/

//          num_links=-1;

//          assemble_intra_lists(p_template, *p_num_template_atoms, &num_angles_list,
//                                  &num_torsions_list, &num_vdws_list);

//          }
         }

/**********************************************************************/
/******** action 7 => Energy minimise the template           **********/
/********             Must be last option until PBC minimise **********/
/********             option is allowed                      **********/
/**********************************************************************/

		  else if (action == MINIMIZE_ACTION)
		    {

             if (DEBUG) printf("DEBUG>> Trying Minimising in gas phase\n");
             printf("ERROR: Attempt to use minimiser action when this is not available\n");
             printf("ERROR: in multi-molecule version.\n");
             exit(0);

//      relabel(p_template, *p_num_template_atoms, &num_elems_used[0], TRUE);

/******************************************************************/
/**************minimize the template w. r. to itself (no pore) ****/
/******************************************************************/


//           if ((strcmp(minimizer_name,DISCOVER_MINIMIZER) == 0) ||
//                      (strcmp(minimizer_name,C2DISCOVER_MINIMIZER) == 0))
//                {
//                 strcpy(discover_root, template_strategy_file);
//                  if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';
//                }


// 		 minimize_molecule(p_template, *p_num_template_atoms,
//                                discover_root, p_need_grad, p_grad, is_empty_pore); 

//                calculate_energy(&pore[0],num_pore_atoms, 
//                                 guest_ptrs, *p_num_seed_mols,
//                                 p_guest_demarc,
//                                 p_molmol_list, num_molmol_list,
//                                 p_kvecs, p_kvec2, p_gvec2,
//                                 num_kvecs, p_cos_sum, p_sin_sum, 
//                                 p_need_grad, p_grad,
//                                 symm_set, num_symm_ops,
//                                 have_comb_rules, is_empty_pore);

// 		 have_minimized = TRUE;

/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/
//               if (symm_set)
//                   {
//                       do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
//                                             p_guest_demarc, &symm[0], num_symm_ops);
//                   }

             }

/**********************************************************************/
/******** action 8 => Energy minimise the template in the pore ********/
/********             Must be last option until PBC minimise **********/
/********             option is allowed                      **********/
/**********************************************************************/

     else if (action == MINIMIZE_INPORE_ACTION)
       {
         if (DEBUG) printf("DEBUG>> Trying Minimising in environment of pore\n");
             printf("ERROR: Attempt to use minimiser inpore action when this is not available\n");
             printf("ERROR: in multi-molecule version.\n");
             exit(0);

//     relabel(p_template, *p_num_template_atoms, &num_elems_used[0], TRUE);

/*******************************************/
/*** minimize the template in the pore) ****/
/*******************************************/

//      strcpy(discover_root, inpore_strategy_file);
//      if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';

//      minimize_inpore(&pore[0], p_template, 
//                      *p_num_template_atoms, discover_root,
//                      p_kvecs, p_kvec2, p_gvec2, num_kvecs,
//                      p_cos_sum, p_sin_sum, p_need_grad, p_grad,
//                      have_comb_rules, is_empty_pore);

/*************************************************************************/
/******** up date symmetry related atoms *********************************/
/*************************************************************************/

//      if (symm_set)
//         {
//             do_symmetry_all(guest_ptrs, *p_num_seed_mols,      
//                                   p_guest_demarc, &symm[0], num_symm_ops);
//         }

//      calculate_energy(&pore[0],num_pore_atoms, 
//                       guest_ptrs, *p_num_seed_mols,
//                       p_guest_demarc,
//                       p_molmol_list, num_molmol_list,
//                       p_kvecs, p_kvec2, p_gvec2,
//                       num_kvecs, p_cos_sum, p_sin_sum, 
//                       p_need_grad, p_grad,
//                       symm_set, num_symm_ops,
//                       have_comb_rules, is_empty_pore);

//     if (is_poly)
//       {
//        fprintf(output_fp, "Reporting from minimise in pore action, poly\n");
//        p_demarc=p_guest_demarc;
//        for (imol=0; imol < *p_num_seed_mols; imol++)
//          {
//             if (symm_set) index=imol*(num_symm_ops+2);
//                                           else index=imol;
//
//             fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

//             print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
//                          guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
//                          seed_number, build_number);

//             p_demarc++;
//          }
//       }
//     else
//       {
/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/
//         fprintf(output_fp, "Reporting from minimise in pore action\n");
//         p_demarc=p_guest_demarc;
//         for (imol=0; imol < *p_num_seed_mols; imol++)
//           {
//              if (symm_set) index=imol*(num_symm_ops+2);
//                                                else index=imol;
//
//              fprintf(output_fp, "Energy information for molecule %d\n", imol);
//
/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

//              print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
//                           guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
//                           seed_number, build_number);

//              p_demarc++;
//           }
//       }


//      have_minimized = TRUE;
       }
/**********************************************************************/
/******** action 9 => Create a new guest molecule, for GCMC     *******/
/********             Added August 2014, Dave Willock           *******/
/**********************************************************************/
     else if (action == CREATE_GUEST_ACTION)
       {
         
         have_created=TRUE;
       }
  
/**********************************************************************/
/******** action 10=> Remove an entire guest, for GCMC         ********/
/********             Added August 2014, Dave Willock           *******/
/**********************************************************************/
     else if (action == REMOVE_GUEST_ACTION)
       {

         have_removed=TRUE;
       }

/****************************************************************************/
/********** End of action list **********************************************/
/*************************************************************************/

//     printf("----oooo----oooo---- Dropped out of end of action list-----oooo----oooo----\n");

/****************************************************************************/
/**** Decide if we need a test  *********************************************/
/****************************************************************************/
         need_test= real_random(1) < prob_test;

/****************************************************************************/
/*********** Peek at current template if time to do so State energies too****/
/****************************************************************************/
                peek_now--;
/****** If peek has zeroed out wait for the next accepted move to reset *****/
                if (peek_now <= 0) peek_now=0;

/****** Now need to see if we should test the current structure ************/
/******                                                         ************/
/****** Two main cases:                                         ************/
/******                                                         ************/
/****** Template design                                         ************/
/****** Should only test new additions if need_test is set.     ************/
/****** otherwise just report energy.                           ************/
/******                                                         ************/
/****** pure MC without modifications to template               ************/
/****** should have prob_test == 1.0, build weight = 0          ************/
/****** MC test will have been done in routines so just need    ************/
/****** to report energies and do any required monitoring.      ************/

/****** Always test last pass !! Bug fix Feb. 99 DJW ***********************/
/****** Always test annealing run if user has set prob_test to 100%. 07 KJ**/
/****** Note that prob_test is 0-1 as reader re-normalises the % entered  **/
/****** by the user.                                               07 DJW **/

     need_test= need_test || (have_built && iloop == num_modify_attempts);

//     if (need_test) printf("DEBUG>> iloop %d needs testing\n", iloop);

/*** Template design case **************************************************/

     if (need_test && have_built && prob_test < 1.0)
       {
//         printf("DEBUG>> TESTING Current template\n");
/***************************************************************************/
/******* calculate cost function returns false if the new template *********/
/******* better than the old                                       *********/
/******* Using prob_test to allow the user to specify pure MC      *********/
/******* in this case need_test should be true if the move         *********/
/******* were accepted in shake, rock or twist and so no further   *********/
/******* MC test should be applied.                                *********/
/******* Dave Willock, Dec. 06                                     *********/
/***************************************************************************/

/**** cost_function will do intra_energy so make sure lists are up to date ****/
/**** and check for bumps.                                                 ****/
/**** Dave Willock March 07                                                ****/

         if (which_test == NON_BOND_COST || which_test == ENERGY_COST )
//             assemble_intra_lists(guest_ptrs, *p_num_seed_mols, p_guest_demarc,
//                                  num_symm_ops, p_angles_list_ptrs,
//                                  p_torsions_list_ptrs, p_vdw_list_ptrs,
//                                  p_links_list_ptrs,
//                                  &num_angles_list[0], &num_torsions_list[0],
//                                  &num_vdws_list[0], &num_links_list[0], just_count);


/**** Call cost_function, this will update the energy terms and return a TRUE if there are serious ***/
/**** steric problems.                                                                             ***/

//         printf("Calling cost function with:\n");
//         printf("%d guests\n", *p_num_guests);
//         printf("%d in molmol list\n", num_molmol_list);
//         if (is_empty_pore) printf("Using an empty pore\n");
//                      else  printf("Pore contains a host\n");

         reject_template = cost_function(guest_ptrs, *p_num_guests, p_guest_demarc,
                                         pore, &box_limits[0],
                                         p_molmol_list, num_molmol_list,
                                         which_test, &cc_tot,
                                         FALSE, &dummy_int,
                                         p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                         p_cos_sum, p_sin_sum,
                                         p_need_grad, p_grad,
                                         have_comb_rules,
                                         is_empty_pore, 
                                         p_angles_list_ptrs, &num_angles_list[0],
                                         p_torsions_list_ptrs,&num_torsions_list[0], 
                                         p_vdw_list_ptrs,&num_vdws_list[0],
                                         TRUE);


         *p_average_cc = cc_tot/ p_which_demarc->num;

         if (reject_template) 
           {
             printf("DEBUG>> This has intra_bumps or is out of the box!!\n");
           }
//        else
//         {
//           printf("DEBUG>> This has passed intra_bumps check\n");
//         }
       }

/*** Pure MC case or incidental move during template design *****/
     else if (done_move)
       {
         if (DEBUG) fprintf(output_fp,"done_move is accepted so processing energy\n");

//         printf("\nDEBUG>> iloop=%d, done_move is accepted so processing energy\n\n", iloop);

         if (which_test == NON_BOND_COST || which_test == ENERGY_COST)
           {

/*** If have_built we are moving a new build around before it has been assessed ****/

               if (have_built)
                 {
                 }
               else
                 {

//                   printf("Previous total energy: %10.6f\n", total_energy);
//                   printf(".....with non_bond = %10.6f, charges = %10.6f, hbond = %10.6f, restraint = %10.6f ",
//                                      interaction_energy[0].non_bonded,
//                                      interaction_energy[0].charges,
//                                      interaction_energy[0].hbond,
//                                      interaction_energy[0].restraint);
//                   printf(" and intra = %10.6f\n", intra_energy[0].total);


                   total_inter=0.0; total_intra=0.0; total_energy=0.0; *p_best_energy=0.0;
                   for (imol=0; imol < *p_num_seed_mols; imol++)
                     {
                       inter_term=    interaction_energy[imol].non_bonded
                                     +interaction_energy[imol].charges
                                     +interaction_energy[imol].hbond 
                                     +interaction_energy[imol].restraint;

                       intra_term=  intra_energy[imol].total;

                       total_inter += inter_term;
                       total_intra += intra_term;

                       total_energy+= inter_term + intra_term;

                       *p_best_energy+=   inter_term
                                         +intra_energy[imol].vdw;
                    }

//                   printf("After move total energy: %10.6f\n", total_energy);
//                   printf(".....with non_bond = %10.6f, charges = %10.6f, hbond = %10.6f, restraint = %10.6f ",
//                                     interaction_energy[0].non_bonded,
//                                      interaction_energy[0].charges,
//                                      interaction_energy[0].hbond,
//                                      interaction_energy[0].restraint);
//                   printf(" and intra = %10.6f\n", intra_energy[0].total);

/*** For MC report running averages ***/

                  if (last_build_number != build_number)
                    {
                      num_mc_accepted=0;
                      energy_tot_grand=0.0;
                      energy_tot_grand2=0.0;
                      energy_inter_grand=0.0;
                      energy_inter_grand2=0.0;
                      last_build_number=build_number;
                    }

                  if (!peek_now)
                    {
                  
                  if (which_test == STERIC_COST)
                     {
                       sprintf(&info[0],"ANIMATE: pass %d action = %d Av. cl. cont= %7.3f",
                                                                   iloop, action, *p_average_cc );
                     }
                   else if (which_test == NON_BOND_COST || which_test == ENERGY_COST)
                     {
                       sprintf(&info[0],"ANIMATE: accepted move %d act= %d Int.E = %7.3f kcal, tot.=%7.3f",
                                     num_mc_accepted, action, total_inter, total_energy );

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/
                      printf("Pure MC energy output(iloop=%d) : inter %10.6f intra %10.6f tot : %10.6f\n",
                                     iloop, total_inter, total_intra, total_energy );
//                      fprintf(output_fp, "Pure MC energy output : inter %10.6f intra %10.6f tot : %10.6f\n",
//                                     total_inter, total_intra, total_energy );

                      p_demarc=p_guest_demarc;
                      for (imol=0; imol < *p_num_seed_mols; imol++)
                        {
                           if (symm_set) index=imol*(num_symm_ops+2);
                                                         else index=imol;

                           fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

                           print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                                        guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                                        imol, build_number);

                           p_demarc++;
                        }
                     }

                    animate_this =TRUE;

                  }
/****************************************************************************/
/*** Carry out any monitoring operations required ***************************/
/****************************************************************************/

                       if (need_monitors)
                         {
                            printf("ERROR: Monitoring options not implemented for multi-molecule version\n");
                            exit(0);

//                          for (host_mon=0; host_mon<=num_host_mon; host_mon++)
//                            {
//                              host_ind = p_monitored->host_list[host_mon];
//
//                              for (guest_mon=0; guest_mon<=num_guest_mon; guest_mon++)
//                                {
//                                  guest_ind = p_monitored->guest_list[guest_mon];
//
//                                  p_guest = p_template+guest_ind;
//
//                                  dist2= atom_separation_squared(&pore[host_ind], p_guest, pbc);

/*** hardwired to contacts less than 4.0A ***/
//                                  if ( dist2 <= 16.0 )
//                                    {
//                                       p_monitored->entry[host_mon][guest_mon]++;
//                                    }
//                                }
//                            }
                         }
                    }
           }
      }

/*** Pure MC case must record the energy of the existing structure again ****/
/*** if the move has been rejected. Always record new moves too.         ****/
/*** Only accept reasonable structures into average, this avoids large   ****/
/*** Clashes at the beginning of the run messing up the calculation      ****/
/*** Dave Willock, May 2013.                                             ****/

if (prob_test >= 1.0 && total_inter < 10.0)
   {
                  num_mc_accepted++;
                  energy_tot_grand+=total_energy;
                  energy_inter_grand+=total_inter;

                  energy_tot_grand2+=total_energy*total_energy;
                  energy_inter_grand2+=total_inter*total_inter;

                  av_etot  = energy_tot_grand/num_mc_accepted;
                  av_einter= energy_inter_grand/num_mc_accepted;

                  stddev_etot = sqrt(energy_tot_grand2/num_mc_accepted - av_etot*av_etot);
                  stddev_einter = sqrt(energy_inter_grand2/num_mc_accepted - av_einter*av_einter);

                  entry_no++;
                  if (total_inter > 100) fprintf(output_fp,"SEEMS A BIT BIGGGGGGGG.......\n");
                  fprintf(output_fp,"%d seed: %d bld %d trials %d Accptd %d Einter %5.2f kcal/mol",
                                                                           entry_no, seed_number, 
                                                                           build_number, iloop+1, num_mc_accepted,
                                                                           total_inter);

                  fprintf(output_fp," (av %5.2f stddev: %5.2f), tot: %5.2f (av %5.2f std dev %5.2f) success/trials: %6.4f\n",
                                   av_einter, stddev_einter, total_energy, av_etot, stddev_etot, 
                                   (double) num_mc_accepted/ (double) (iloop+1)); 
/*** Write accepted energies to csv file ****/

                  if ( peek_now_energy == 0)
                    {
//                  printf("writing energy to energy csv file peek_freq = %d \n", peek_freq); 

/*** Record data on progress ****/
                  fprintf(energy_csv_fp,"%d , seed: , %d  ,bld , %d , trials , %d  ,Accptd , %d  ,",
                                                                           entry_no, seed_number, 
                                                                           build_number, iloop+1, num_mc_accepted);
/*** Record system total energies ***/
                  fprintf(energy_csv_fp,"Einter , %9.4f , kcal/mol , (av , %9.4f , stddev: , %9.4f ,), total: , %9.4f , (av , %9.4f , stddev  ,%9.4f,), ",
                                                                           total_inter, av_einter, stddev_einter, total_energy, av_etot, stddev_etot ); 

/*** Record per molecule energies ***/
                  fprintf(energy_csv_fp," PER MOLECULE (for %d mols) : Einter , %9.4f , kcal/mol , (av , %9.4f , stddev: , %9.4f ,), total: , %9.4f , (av , %9.4f , stddev  ,%9.4f,),",
                                                                            *p_num_seed_mols, total_inter/(double) *p_num_seed_mols, av_einter/ (double) *p_num_seed_mols , 
                                                                           stddev_einter/ (double) *p_num_seed_mols , total_energy/ (double) *p_num_seed_mols , 
                                                                           av_etot/(double) *p_num_seed_mols , stddev_etot/ (double) *p_num_seed_mols ); 

/*** Finish with success/fail info **/
                  fprintf(energy_csv_fp," success/trials:, %6.4f\n", (double) num_mc_accepted/ (double) (iloop+1)); 

                      peek_now_energy = peek_freq;
                    }
                  else peek_now_energy--;
             /*   fprintf(energy_csv_fp, "%d, %10.6f\n", num_mc_accepted, total_energy); */

    }


/****************************************************************************/
/*** For template building reject template will be TRUE if bumps detected ***/
/****************************************************************************/
/******* Test to see if we have reached the good template criterion *********/
/****************************************************************************/

            if (need_test && (have_built || have_made_ring) && !reject_template)
              {
//               printf("No intra bump or box problems here\n");

                 if (which_test == STERIC_COST)
                    {

                       *p_average_cc = cc_tot/p_which_demarc->num;

                       printf("Average close contact distance for current template = %10.6f\n",
                                                                               *p_average_cc);
		                animate_this = TRUE;

                       if (*p_average_cc < stop_ctf)
                          {
                             printf(" Thats a good enough template for now!!!!\n");
                             finished= TRUE;
                          }
                    }
                else if (which_test == NON_BOND_COST || which_test == ENERGY_COST )
                    {
                      calculate_energy(&pore[0],num_pore_atoms, 
                                       guest_ptrs, *p_num_seed_mols,
                                       p_guest_demarc,
                                       p_molmol_list, num_molmol_list,
                                       p_kvecs, p_kvec2, p_gvec2,
                                       num_kvecs, p_cos_sum, p_sin_sum, 
                                       p_need_grad, p_grad,
                                       have_comb_rules, is_empty_pore);

/*** do not include intra vdw in total energy after a build this is in the cost function for MC moves ***/
                       if (have_built || have_made_ring )
                          {
                             total_energy=0.0;
                             for (imol=0; imol < *p_num_seed_mols; imol++)
                               {
                                    total_energy+=  interaction_energy[imol].non_bonded
                                                   +interaction_energy[imol].charges
                                                   +interaction_energy[imol].hbond 
                                                   +interaction_energy[imol].restraint;
                               }
                             printf("DEBUG>> Updating best_total from best_energy = %10.6f\n", *p_best_energy);

                             best_total = *p_best_energy; 
                          }
                       else
                          {
                             total_energy=0.0;
                             for (imol=0; imol < *p_num_seed_mols; imol++)
                               {
                                 total_energy+= interaction_energy[imol].non_bonded
                                               +interaction_energy[imol].charges
                                               +interaction_energy[imol].hbond 
                                               +interaction_energy[imol].restraint
                                               +intra_energy[imol].total;
                               }
                             best_total = *p_best_energy + best_intra;
                          }

/***** Stop re-setting of reject_template if this is just MC run ( prob_test >= 1.0 ) ****/
                   printf("DEBUG>> prob_test = %10.6f\n", prob_test);
                   if (prob_test < 1.0)
                     {
                        if (DEBUG) printf("Now new total is %10.6f compared with %10.6f from last time so ", total_energy, best_total);

                       if(total_energy > best_total) 
                         {
                           printf(" this has gone up may reject\n");
                           reject_template=TRUE;
                         }
                       else 
                         {
                           printf(" this has gone down, should accept\n");
                           reject_template=FALSE;
                         }
                     }

/****** Do proper Monte Carlo at user defined temperature ******/

                   if ( reject_template )
                     {

                       printf("But give boltzmann a go first....\n");
                       boltz_fact = exp( - (total_energy - best_total)
                                                / (BOLTZ * temperature));

                       if (DEBUG) printf("Energy now = %10.6f best = %10.6f applying boltzman\n", 
                                                                    total_energy, best_total); 

                       if(real_random(1) < boltz_fact) reject_template=FALSE;
                     }
                   else
                     {
                       if (DEBUG) printf("Energy now = %10.6f best = %10.6f accepting....\n", 
                                                                      total_energy, best_total);
                     }
                }
            }

/****************************************************************************/
/*** Now deal with acceptance or rejection of the latest structural change **/
/****************************************************************************/
/*** Code re-structured July 09 Dave Willock. *******************************/
/*** earlier version made two tests for !reject_template with the ***********/
/*** additional condition for builds etc below.                   ***********/
/*** This method is neater and allows print_energy to pick up     ***********/
/*** correct version of build_number.                             ***********/
/****************************************************************************/
                
/**** Bug fix March 07: replace else with new test so that up hill accepts from above are dealt with ***/

            if (need_test && (have_built || have_made_ring) && !reject_template)
                     {

                 if (have_built || have_made_ring || *p_just_position)
                   {
                      DEBUG=TRUE;
                      if (DEBUG) printf("PASSED ");
                      if (DEBUG) fprintf(output_fp,"PASSED\n");
                      if (DEBUG && have_built) printf("A build has worked. ");
                      if (DEBUG && have_made_ring) printf("A ring maker move has worked. ");
                      if (DEBUG && *p_just_position) printf("Just positioning. ");

//          if (DEBUG)  printf("\nDEBUG: %d Hydrogens in new guest structure, indices and weights.........\n",
//                                                                                *(p_num_guest_hyds+which_built_mol)+1);                       
//          if (DEBUG)  printf("DEBUG: sum_hyd_weights %d\n", *(p_sum_hyd_weights+which_built_mol));                             
//          for (ihyd=0; ihyd <= *(p_num_guest_hyds+which_built_mol); ihyd++)
//             {
//                printf("DEBUG>> hyd index : %d hyd_weight : %d\n", guest_hyd_list_ptrs[which_built_mol][ihyd],
//                                                                   guest_hyd_weights_ptrs[which_built_mol][ihyd]);
//             }

                      DEBUG=FALSE;

/*** Rest warnings for atoms clashing ***************************************/
                      told_of_clash=FALSE;

/****************************************************************************/
/*********** Initialise animation files for this template *******************/
/****************************************************************************/
                      if (logfile_needed)
                         {
                            fclose(anim_read_fp);
                            fclose(anim_show_fp); 
                         }
                      if (animate_flag !=FALSE)
                         {
                            fclose(anim_file_fp);
                            fclose(anim_pdb_fp);
                            fclose(anim_xyz_fp);
                         }
                      build_number++;

                      if (animate_flag != 0)
                          {
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
                             if (use_number > 0)
                               {
                                 printf("ERROR: This call wacko, use_number = %d\n", use_number);
                                 exit(0);
                               }
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
/***DEBUG DEBUG DEBUG ***/
                       printf("Starting new animation file\n");

                       initialise_animation(&this_animation_file[0], p_animation_file,
                                            seed_type, seed_number, use_number,
                                            build_number);
/*************************************/
/** Additions for pdb file writer ****/
/*************************************/

                       if (seed_type != MOLE)
                         {
                           sprintf(pdb_filename,"%s%d_%d_%d.pdb",p_animation_file,
                                               use_number, seed_number, build_number);
                           sprintf(xyz_filename,"%s%d_%d_%d.xyz",p_animation_file,
                                               use_number, seed_number, build_number);
                         }
                       else
                         {
                           sprintf(pdb_filename,"%s%d_%d.pdb",p_animation_file,
                                                seed_number, build_number);
                           sprintf(xyz_filename,"%s%d_%d.xyz",p_animation_file,
                                                seed_number, build_number);
                         }

                          if (open_file(&anim_pdb_fp, pdb_filename, "w")== EXIT_FAILURE)
                             {
                                printf("ERROR>> Cannot open pdb file %s for writing\n", pdb_filename);
                                exit(0);
                             }

                          if (open_file(&anim_xyz_fp, xyz_filename, "w")== EXIT_FAILURE)
                             {
                                printf("ERROR>> Cannot open xyz file %s for writing\n", xyz_filename);
                                exit(0);
                             }
  
                          super[0]=1;
                          super[1]=1;
                          super[2]=1;

/*********** Always animate first frame, added Nov 13th 2006 Dave Willock ***/

                          animate_this = TRUE;
                       }


/************************ Update Statistics ********************************/
                         statistics[action].accepted++;
                         statistics[action].accepted_after_build += have_built;

/****************************************************************************/
/******* Update template types **********************************************/
/******* Added the p_type variable 28th Aug. 1997 since the whole array *****/
/******* was not being tested !!  DJW ***************************************/
/****************************************************************************/

                         if (have_made_ring)
                           {
// Needs updating for multi-molecule version to identify which molecule has had a ring 
// added.
//
//                           p_type= p_template_types;
//                           for (itype=0; itype <= *p_num_template_types; itype++)
//                             {
//                               if (p_type->atom_type[0] == 'H' &&
//                                   p_type->atom_type[0] == '\0' )
//                                 {
//                                    p_type->num -= 2;
//                                 }
//                               p_type++;
//                             }
                           }
                         else if (have_built)
                           {
//                            p_atom= p_template+num_old_template_atoms;
//                            relabel(p_atom, num_to_rot, &num_elems_used[0], FALSE);
//
//                            start_types= (p_frag_types_list+last_frag_added)->start;
//
//                            add_types( p_template_types, p_num_template_types,
//                                       p_frag_types+start_types, 
//                                       (p_frag_types_list+last_frag_added)->num);

/**** Keep track of the fragments added ****/
                              composition[last_frag_added]++;
 
                              fprintf(output_fp, "fragments added so far:\n");
                              for (ifrag= 1; ifrag <= number_of_fragments; ifrag++)
                                {
                                  fprintf(output_fp, "fragment %d : %d\n", 
                                                    ifrag, composition[ifrag]);
                                }

/**** Allow further building *****/

                              have_built = FALSE;
                           }

                       }

                       *p_best_energy= 0.0; best_restraint=0.0; best_intra=0.0;
                       for (imol=0; imol < *p_num_seed_mols; imol++)
                         {
                           *p_best_energy+= interaction_energy[imol].non_bonded
                                           +interaction_energy[imol].charges
                                           +interaction_energy[imol].hbond;

                           best_restraint+= interaction_energy[imol].restraint;

                           best_intra += intra_energy[imol].total;
                         }

/*                       printf("\nEnergy for ensemble is : %10.6f\n", total_energy); */

                /*************************************************************/
                /*** KJ May 07 - now if this was the best ever template,   ***/
                /***             then store this and the energy            ***/
                /*************************************************************/
                if (optimise_post_mc)
                {
                /******was this the best ever template?******/
                if (best_total < best_ever)
                   {
                    if(DEBUG)
                        {
                printf ("this template %d being made, attempt no %d\n", no_annealed_temps, iloop);
                printf ("this template is lower in energy than the best ever template");
                printf(" - shall have to keep\n");
                printf ("this template = %10.6f, best_ever = %10.6f\n", best_total, best_ever);
                        }
                   best_ever = best_total;

                    /***** and copy the template *******/

                   } /* end if (best_total < best_ever) */

                } /* end if(optimise_post_mc)*/

/*******************************************************/

                       if (which_test == STERIC_COST)
                          {
                             printf("Average close contact distance for current template = %10.6f\n",
                                                                                   *p_average_cc);
                          }
                       else if (which_test == NON_BOND_COST || which_test == ENERGY_COST)
                         {
                            if (verbose)
                              {
                                if (!peek_now) fprintf(output_fp, "Energy summary after accepting at step %d\n", 
                                                                               iloop);
                                if (is_poly)
                                  {
                                     fprintf(output_fp, "inside reject_template, poly\n");

                                     p_demarc=p_guest_demarc;
                                     for (imol=0; imol < *p_num_seed_mols; imol++)
                                       {
                                          if (symm_set) index=imol*(num_symm_ops+2);
                                                                        else index=imol;

                                          fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

                                          print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                                                       guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                                                       imol, build_number);

                                          p_demarc++;
                                       }

                                    if (have_built)
                                      {
                                        fprintf(output_fp, "Density of system = %10.6f a.m.u. / cubic A\n", 
                                                       (total_mass_temp+pore_total_mass)/cell_volume);

                                        fprintf(output_fp, "                 or %10.6f g / cubic cm\n", 
                                                       1.66*(total_mass_temp+pore_total_mass)/cell_volume);
                                      }
                                  }
                             else
                               {
/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/
                                 fprintf(output_fp,"DEBUG>> Normal print of energy\n");
                                 p_demarc=p_guest_demarc;
                                 for (imol=0; imol < *p_num_seed_mols; imol++)
                                   {
                                      if (symm_set) index=imol*(num_symm_ops+2);
                                                                    else index=imol;

                                      fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

                                      print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                                                   guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                                                   imol, build_number);

                                      p_demarc++;
                                   }

                               }
                           }
                 
/*                         printf("Last accepted energy before this point = %10.6f ", best_total); */
/*                         printf(" Energy now= %10.6f\n", total_energy); */
                      }


/****************Only animate if peek says you should ***********************/

                       if (!peek_now) 
                         {
                            if (DEBUG) printf("peek_now zeroed out\n");
		                     animate_this = TRUE;

/*******************open file************************************************/

                            if (!(peek_fp = fopen(PEEK_FILENAME,"w")))
                              {
                                 fprintf(output_fp, 
                                         "ZEBEDDE ERROR: Unable to open Peek file for writing\n");

                                 fflush(output_fp);
                                 exit(-1);
                              }
                            else
                              {
                                 print_peek_file(guest_ptrs, *p_num_seed_mols, p_guest_demarc, 
                                                 iloop, peek_fp);
                              }
                          }
                       
                       if (total_energy < stop_ctf)
                         {
                             printf(" Thats a good enough template for now!!!!\n");
                             finished= TRUE;
                         }
                     }

/****************************************************************************/
/*** For rejected templates restore original ********************************/
/****************************************************************************/

        else if (need_test && reject_template && (have_built || have_made_ring))
          {
            DEBUG=TRUE;
            if (DEBUG)
             {
                printf("FAILED current energy %10.6f best %10.6f\n", total_energy, best_total);
                printf("FAILED have to reinstate guest %d\n", which_built_mol);
                fprintf(output_fp,"FAILED\n");
             }
            DEBUG=FALSE;

/*********************************************************************************/
/******** Reset to last template now should do this with realloc and unshuffle ***/
/*********************************************************************************/
/*********************************************************************************/
/*** THIS IS WHAT HAS TO BE UNDONE !!!! ******************************************/
/*********************************************************************************/

/*********************************************************************************/
/** Realloc arrays back to original sizes ****************************************/
/*********************************************************************************/

     p_demarc=p_guest_demarc+which_built_mol;
     *p_demarc= *p_old_demarc;

//     printf("reallocing space for guest %d to size %d\n", which_built_mol,  p_demarc->num );

     guest_ptrs[which_built_mol]= (atom*) realloc( guest_ptrs[which_built_mol], p_demarc->num * sizeof(atom));

/***************************************************************************/
/***** Restore information from old template to carry out rejection ********/
/***************************************************************************/

//    printf("DEBUG: Copying atoms back from old version...........\n");		 		   

/***************************************************************************/
/**** Reinstate the guest itself *******************************************/
/***************************************************************************/

    p_old_atom= p_old_guest; 

    for (iatom=0; iatom <=p_demarc->end; iatom++)
      {
         guest_ptrs[which_built_mol][iatom] = *p_old_atom;

//         if (iatom < 10) printf("DEBUG>> Reinstated  old atom %s %10.6f  %10.6f  %10.6f \n", p_old_atom->label,
//                                                                                             p_old_atom->x,
//                                                                                             p_old_atom->y,
//                                                                                             p_old_atom->z);
         p_old_atom++;
     }

//    printf("DEBUG>> freeing the old_guest\n");
    free(p_old_guest); 

/***************************************************************************/
/*** Link list is now more than inter-fragment so need num_links saved too */
/*** added Aug. 09 Dave Willock                                            */
/***************************************************************************/

//     if (DEBUG)  printf("DEBUG: Reallocing links list.........................\n");		 		   

     num_links_list[which_built_mol]=num_old_links; 
    
/**** Since the no-links value is -1 add 2 here to make sure some space is allocated ****/

     p_links_list_ptrs[which_built_mol]= (links*) realloc( p_links_list_ptrs[which_built_mol], 
                                                               (num_links_list[which_built_mol]+2) * sizeof(links));
 
//     if (DEBUG)  printf("DEBUG: Copying %d links back from old version...........\n", 
//                                                                    num_links_list[which_built_mol]);		 		   
     p_this_link=p_old_links; 
     for (ilink=0; ilink <= num_links_list[which_built_mol]; ilink++)            
       {                                                   
          p_links_list_ptrs[which_built_mol][ilink]= *p_this_link;	          
//          printf("Start: %d End: %d\n", p_this_link->start, p_this_link->end);
          p_this_link++;
       }                                                     

//    printf("DEBUG>> freeing the old_links\n");
    free(p_old_links); 
/******************************************************************************/
/******** Reinstate old hydrogen                        ***********************/
/******** and hydrogen weights to carry out rejection   ***********************/
/******************************************************************************/

//     if (DEBUG)  printf("DEBUG: Reallocing hydrogen list and H weights list...\n");		 		   

          *(p_num_guest_hyds+which_built_mol) = num_old_hyds;    
          *(p_sum_hyd_weights+which_built_mol) = sum_old_hyd_weights; 

          guest_hyd_list_ptrs[which_built_mol] = (int*) realloc(  guest_hyd_list_ptrs[which_built_mol],
                                                              ( *(p_num_guest_hyds+which_built_mol) + 1)* sizeof(int));

          guest_hyd_weights_ptrs[which_built_mol] = (int*) realloc(  guest_hyd_weights_ptrs[which_built_mol],
                                                              ( *(p_num_guest_hyds+which_built_mol) + 1)* sizeof(int));

          if (DEBUG)  printf("DEBUG: Copying hydrogen %d indices and weights back.........\n", 
                                                                                *(p_num_guest_hyds+which_built_mol)+1);		 		   
          if (DEBUG)  printf("DEBUG: Copied back sum_hyd_weights as %d\n", *(p_sum_hyd_weights+which_built_mol));		 		   

          p_int= p_old_hyd_list;
          p_weight= p_old_hyd_weights;
          for (ihyd=0; ihyd <= *(p_num_guest_hyds+which_built_mol); ihyd++) 
             {                                        
                guest_hyd_list_ptrs[which_built_mol][ihyd] = *p_int;
                guest_hyd_weights_ptrs[which_built_mol][ihyd] = *p_weight; 

//                printf("DEBUG>> hyd index : %d hyd_weight : %d  %s elem %s\n", guest_hyd_list_ptrs[which_built_mol][ihyd],
//                                                                               guest_hyd_weights_ptrs[which_built_mol][ihyd],
//                                                                               guest_ptrs[which_built_mol][guest_hyd_list_ptrs[which_built_mol][ihyd]].label, 
//                                                                               guest_ptrs[which_built_mol][guest_hyd_list_ptrs[which_built_mol][ihyd]].elem);
                                                    
                p_int++;                             
                p_weight++;                           
             }                                         

//    printf("DEBUG>> freeing the old_hyds and old_hyd_weight arrays.\n");
    free(p_old_hyd_list); 
    free(p_old_hyd_weights); 
//    printf("DEBUG>> reinstate 3\n");

/***************************************************************************/
/*** Reinstate old types list **********************************************/
/***************************************************************************/

//    printf("DEBUG: Reallocing types list...\n");		 		   

     *p_num_this_types = num_old_types;

     guest_types_ptrs[which_built_mol] = (atom_number*) realloc( guest_types_ptrs[which_built_mol],
                                                                       (*p_num_this_types) * sizeof(atom_number));

     p_type= p_old_types;
     for (itype=0; itype < *p_num_this_types; itype++)
        {
           guest_types_ptrs[which_built_mol][itype] = *p_type;
           p_type++; 
        }

//     printf("Reinstated %d types:\n",  *p_num_this_types);

//     for (itype=0; itype<num_old_types; itype++)
//      {
//         printf("%d %s %d\n", itype, guest_types_ptrs[which_built_mol][itype].atom_type, 
//                                     guest_types_ptrs[which_built_mol][itype].num);
//      }


/*** Update intra-molecular interactions and links list for just this molecule *****/

         assign_stretch(guest_ptrs[which_built_mol], p_demarc->end);

         just_count=TRUE;
         assemble_intra_lists(&guest_ptrs[which_built_mol], 1, p_demarc, 
                              num_symm_ops, &p_angles_list_ptrs[which_built_mol],
                              &p_torsions_list_ptrs[which_built_mol], &p_vdw_list_ptrs[which_built_mol],
                              &p_links_list_ptrs[which_built_mol],
                              &num_angles_counted[which_built_mol], &num_torsions_counted[which_built_mol], 
                              &num_vdws_counted[which_built_mol], &num_links_counted[which_built_mol], just_count);

/**** Need to malloc angles_, torsions_ and vdw_ lists and the link_atoms ***/
/**** Note we only need think about the symmetry unique molecules         ***/

//         printf("Memory re-allocation info for intra-molecular potentials:\n");
//         printf("......just_counting for the rejection of an unwanted addition to molecule %d\n", which_built_mol);
//         printf("molecule    angles    torsions    links    vdws\n");
//               printf(" %d     %d     %d       %d      %d\n", which_built_mol, num_angles_counted[which_built_mol]+1,
//                                                                    num_torsions_counted[which_built_mol]+1,
//                                                                    num_links_counted[which_built_mol]+1,
//                                                                    num_vdws_counted[which_built_mol]+1);

               if (num_angles_counted[which_built_mol] > 0 )
                  p_angles_list_ptrs[which_built_mol]= 
                       (angle_interact_list*) realloc( p_angles_list_ptrs[which_built_mol],
                                                      (10+num_angles_counted[which_built_mol])*sizeof(angle_interact_list));
               else
                   p_angles_list_ptrs[which_built_mol]= (angle_interact_list*) realloc( p_angles_list_ptrs[which_built_mol],
                                                                                        sizeof(angle_interact_list));

               if (num_torsions_counted[which_built_mol] > 0 )
                  p_torsions_list_ptrs[which_built_mol]= 
                      (torsion_interact_list*) realloc( p_torsions_list_ptrs[which_built_mol],
                                                (10+num_torsions_counted[which_built_mol])*sizeof(torsion_interact_list));
               else
                  p_torsions_list_ptrs[which_built_mol]= 
                      (torsion_interact_list*) realloc( p_torsions_list_ptrs[which_built_mol], sizeof(torsion_interact_list));

               if (num_vdws_counted[which_built_mol] > 0 )
                 p_vdw_list_ptrs[which_built_mol]= 
                       (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_built_mol],
                                                    (10+num_vdws_counted[which_built_mol])*sizeof(vdw_interact_list));
               else
                 p_vdw_list_ptrs[which_built_mol]= 
                       (vdw_interact_list*) realloc(p_vdw_list_ptrs[which_built_mol], sizeof(vdw_interact_list));

               if (num_links_counted[which_built_mol] > 0 )
                             p_links_list_ptrs[which_built_mol]= 
                                  (links*) realloc(p_links_list_ptrs[which_built_mol], 
                                                         (10+num_links_counted[which_built_mol])*sizeof(links));
               else
                             p_links_list_ptrs[which_built_mol]= 
                                  (links*) realloc(p_links_list_ptrs[which_built_mol], sizeof(links));

               if (     p_angles_list_ptrs[which_built_mol] == NULL
                   || p_torsions_list_ptrs[which_built_mol] == NULL
                   ||      p_vdw_list_ptrs[which_built_mol] == NULL
                   ||    p_links_list_ptrs[which_built_mol] == NULL )
                 {
                   printf("ERROR: Problem assigning memory for intra molecular potentials or links during rejection step.\n");
                   exit(0);
                 }
                                                          
               num_angles_list[which_built_mol]= num_angles_counted[which_built_mol];
               num_torsions_list[which_built_mol]= num_torsions_counted[which_built_mol];
               num_links_list[which_built_mol]= num_links_counted[which_built_mol];
               num_vdws_list[which_built_mol]= num_vdws_counted[which_built_mol];

               just_count=FALSE;
//               printf("DEBUG: Now assembling intra lists for the just built molecule %d\n",which_built_mol);

               assemble_intra_lists(&guest_ptrs[which_built_mol], 1, p_demarc, 
                                    num_symm_ops, &p_angles_list_ptrs[which_built_mol],
                                    &p_torsions_list_ptrs[which_built_mol], &p_vdw_list_ptrs[which_built_mol],
                                    &p_links_list_ptrs[which_built_mol],
                                    &num_angles_list[which_built_mol], &num_torsions_list[which_built_mol], 
                                    &num_vdws_list[which_built_mol], &num_links_list[which_built_mol], just_count);

//               printf("DEBUG>> back after filling lists during rejection\n");

//               printf("Array filling info for intra-molecular potentials:\n");
//               printf("molecule    angles    torsions    links    vdws\n");
//               printf(" %d     %d     %d       %d      %d\n", which_built_mol, num_angles_list[which_built_mol]+1,
//                                                                    num_torsions_list[which_built_mol]+1,
//                                                                    num_links_list[which_built_mol]+1,
//                                                                    num_vdws_list[which_built_mol]+1);

               mem_prob=FALSE;
               if (num_angles_list[which_built_mol] > num_angles_counted[which_built_mol]) mem_prob=TRUE;
               if (num_torsions_list[which_built_mol] > num_torsions_counted[which_built_mol]) mem_prob=TRUE;
               if (num_links_list[which_built_mol] > num_links_counted[which_built_mol]) mem_prob=TRUE;
               if (num_vdws_list[which_built_mol] > num_vdws_counted[which_built_mol]) mem_prob=TRUE;

              if ( mem_prob )
                {
                   printf("ERROR: Poor matching of counted and assigned array sizes will cause memory problems\n");
                   printf("ERROR: Could be that a potential is not available for one of the atoms sets found.\n");
                   printf("ERROR: Update frc file.\n");
                   exit(0);
                }

//               printf("Current links list for molecule %d:\n", which_built_mol);
//               for (iii=0; iii <= num_links_list[which_built_mol]; iii++)
//                  {
//                     ist=p_links_list_ptrs[which_built_mol][iii].start;
//                     ied=p_links_list_ptrs[which_built_mol][iii].end;
//
//                     printf("%d %d %s %s\n", ist, ied,
//                                            guest_ptrs[which_built_mol][ist].label, 
//                                            guest_ptrs[which_built_mol][ied].label);
//                  }

/**** Assemble molmol interaction indices list ******/
               num_molmol_list = (num_seeds_with_symm+1)*(num_seeds_with_symm+2)/2;

//               printf("DEBUG>> mallocing %d space for molmol list\n", num_molmol_list);

               p_molmol_list= (interaction_indices*) malloc(num_molmol_list*sizeof(interaction_indices));
               
//               printf("Have %d molecules so expect %d in the pair interaction list\n", num_seeds_with_symm, num_molmol_list);

               num_molmol_list= assemble_molmol_list(p_molmol_list, *p_num_seed_mols,
                                                     num_molmol_list);

               if (num_molmol_list<0)
                 {
                   printf("ERROR: Found too few interactions in molmol list\n");
                   exit(0);
                 }
//               else printf("Found the expected number of entries in the molmol list\n");

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

/*********************************************************************************/
/*** Old version below, Old version below DELETE !!!!   **************************/
/*********************************************************************************/

//          *p_num_template_atoms= num_old_template_atoms;
//          *p_num_temp_atoms_with_symm= num_old_template_atoms_with_symm;

/*** Redo centre of mass added June 07, Dave Willock ***/

            centre_of_mass(&(p_which_cofm->v[0]), p_which_totmass, 
                            guest_ptrs[which_built_mol], p_demarc->end, -1 );

/*****************************************************************************/
/******** Restore old atom set ***********************************************/
/******** including symmetry images, added Aug 09 Dave Willock     ***********/
/*****************************************************************************/

//          p_atom= p_template;
//          for (iatom=0; iatom<= *p_num_temp_atoms_with_symm; iatom++)
//            {
//              *p_atom= old_template[iatom];
//              p_atom++;
//            }

/*******************************************************************/
/****Bug fix August 1997 DJW ***************************************/
/****Only do these things if have_built not if have_made_ring! *****/
/*******************************************************************/

	     if (have_built) 
              {
                have_built=FALSE;
                fragments_in_template--;
              }

/****************************************************************************/
/******* Restore old links Really redundant as assemble torsions makes list */
/****************************************************************************/
//           num_links= num_old_links;
//           for (ilink=0; ilink <= num_links; ilink++)
//             {
//                link_atoms[ilink]= old_links[ilink];
//             }

/****************************************************************************/
/******* Restore old hydrogen list ******************************************/
/******* and weights               ******************************************/
/****************************************************************************/

//          *p_num_template_hyds= num_old_template_hyds;
//          *p_sum_temp_hyd_weights= sum_old_temp_hyd_weights;
            
//          p_int = p_template_hyd_list;
//          p_weight= p_temp_hyd_weights;
//          for (ihyd=0; ihyd <= *p_num_template_hyds; ihyd++)
//           {
//              *p_int = old_template_hyd_list[ihyd];
//              *p_weight= old_temp_hyd_weights[ihyd];
//
//              p_int++;
//              p_weight++;
//           }

/****************************************************************************/
/**** Re-do intra molecular interaction lists *******************************/
/**** Note that stretch potential indexing is held in atom structure ********/
/**** so no need to call assign_stretch here.                        ********/
/****************************************************************************/

//          num_links=-1;

//          assemble_intra_lists(p_template, *p_num_template_atoms, &num_angles_list,
//                                  &num_torsions_list, &num_vdws_list);

        have_built= FALSE;
        have_made_ring= FALSE;
        have_minimized = FALSE;
     }

/*****************************************************************************/
/******** only animate if something happened in this loop ********************/
/*****************************************************************************/

     if (DEBUG)
       {
         printf("Would like to animate but peek_now = %d\n",peek_now);
         if (animate_this) printf("Animate flag set\n");
         if (!animate_this) printf("Animate flag not set\n");
       }

     if (animate_flag  != 0 && animate_this) 
	 {
             if (DEBUG) printf("Animating\n"); 

             if (which_test == STERIC_COST)
                {
                  sprintf(&info[0],"ANIMATE: pass %d action = %d Av. cl. cont= %7.3f",
                                iloop, action, *p_average_cc );
                }
              else if (which_test == NON_BOND_COST || which_test == ENERGY_COST)
                {
                  sprintf(&info[0],"ANIMATE:pass %d act= %d Int.E = %7.3f kcal, res.=%7.3f",
                                iloop, action, *p_best_energy, best_restraint );
                }


//        printf("DEBUG1>> Calling animate with %d guests\n", *p_num_seed_mols);
          fflush(stdout);

          animate( guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                   anim_file_fp, anim_read_fp, anim_show_fp, 
                   &info[0], &this_animation_file[0], &num_anime_frames, 
                   pore);

        write_pdb(anim_pdb_fp, guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                  &abc[0], &super[0], &latt_vec[0]);

        write_xyz(anim_xyz_fp, guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                  &abc[0], &super[0], &latt_vec[0]);

          peek_now = peek_freq;
          animate_this = FALSE;
          animated= TRUE;
        }

      iloop++;
    }


/************************************************************/
/**************** End of Main loop **************************/
/************************************************************/

/************************************************************/
/*** If testing and no animations yet do one to get energy **/
/************************************************************/
 if (*p_just_position && !animated)
   {
      if (which_test == STERIC_COST)
        {
           printf("Average close contact distance for current template = %10.6f\n",
                                                           *p_average_cc);

           sprintf(&info[0],"ANIMATE:pass %d act.=%d Av.cl.cont= %7.3f",
                              iloop, action, *p_average_cc );

        }
      else
        {
           sprintf(&info[0],"ANIMATE:pass %d act.=%d Int.E= %7.3f kcal(res = %7.3f)",
                                iloop, action, *p_best_energy, best_restraint );
        }

          printf("DEBUG2>> Calling animate with %d guests\n", *p_num_seed_mols);
          fflush(stdout);

          animate( guest_ptrs, *p_num_seed_mols, p_guest_demarc,
                   anim_file_fp, anim_read_fp, anim_show_fp, 
                   &info[0], &this_animation_file[0], &num_anime_frames, 
                   pore);
   }

printf("leaving main loop: num_modify_attempts = %d finished = %d iloop = %d\n",
                                   num_modify_attempts, finished, iloop);

/* printf("min_inpore test\n");  */
/* strcpy(discover_root, inpore_strategy_file);  */
/* if ((fullstop = strchr(discover_root,'.')) !=NULL) *fullstop = '\0';  */

/**************** close up the necesary files ***************/

if (logfile_needed)
   {
      fclose(anim_read_fp);
      fclose(anim_show_fp); 
   }
if (animate_flag !=FALSE) /*No need to close these as they were never opened! AJWL 08*/
{
fclose(anim_file_fp);
fclose(anim_pdb_fp);
fclose(anim_xyz_fp);
}

/***************************************************************************************************/
/************* if annealing run - give the user the final annealed structure (ie at OK)*************/
/********************************** of this template as gulp input *********************************/
/***************************************************************************************************/

if(optimise_post_mc && anneal_last)
{

fprintf(output_fp,"KJ SA: Annealed structure: %d, and the energy is:\n", no_annealed_temps);
fprintf(output_fp,"and the best was %10.6f\n", best_ever); 



/**** calculate the energetics for best_ever_template ****/
                      calculate_energy(&pore[0],num_pore_atoms, 
                                       guest_ptrs, *p_num_seed_mols,
                                       p_guest_demarc,
                                       p_molmol_list, num_molmol_list,
                                       p_kvecs, p_kvec2, p_gvec2,
                                       num_kvecs, p_cos_sum, p_sin_sum, 
                                       p_need_grad, p_grad,
                                       have_comb_rules, is_empty_pore);

/**** then print these energetics to the output file ****/
                                     fprintf(output_fp, "inside optimise_post_mc, poly\n");

               p_demarc=p_guest_demarc;
               for (imol=0; imol < *p_num_seed_mols; imol++)
                 {
                    if (symm_set) index=imol*(num_symm_ops+2);
                                                  else index=imol;

                    fprintf(output_fp, "Energy information for molecule %d\n", imol);

/**** madalung, the second to last logical, should be a user keyword, removed to reduce output **/

                    print_energy(output_fp,  &interaction_energy[imol], &intra_energy[imol],
                                 guest_ptrs[index], p_demarc->end, FALSE, FALSE, FALSE,
                                 imol, build_number);

                    p_demarc++;
                 }

//printf("KJDB: neighbours now annealed....\n");
//print_neighbours(&best_ever_template[0], *p_num_template_atoms, stdout);


fprintf(output_fp,"SA: Annealed structure: %d, energy above here\n", no_annealed_temps);

print_dashes(80,output_fp);

/*** now save the best template position ******/
if (DEBUG) printf("Creating a gulp input file of the annealed structure, template no: %d\n", no_annealed_temps);

final_annealed_gulp=TRUE;
}

/**** Close files that have been used in this run ***/
if (energy_file_is_open)
  {
     fflush(output_fp);  /** debugging problem with output **/
     fflush(energy_csv_fp);
     printf("Closing energy csv file\n");
     fclose(energy_csv_fp);
     printf("file_closed..\n");
     energy_file_is_open=FALSE;
  }

  if (need_monitors)
   {

     for (host_mon=0; host_mon<=num_host_mon; host_mon++)
       {
          host_ind = p_monitored->host_list[host_mon];

          for (guest_mon=0; guest_mon<=num_guest_mon; guest_mon++)
           {
               guest_ind = p_monitored->guest_list[guest_mon];

               fprintf(csv_fp,"%d, ",p_monitored->entry[host_mon][guest_mon]);
           }
          fprintf(csv_fp,"\n");
       }
     fclose(csv_fp);
   }
/*********** be nice and give some info and statistics on the run **********/

printf("Going to print stats\n");
print_statistics(&statistics[0], NUMBER_OF_ACTIONS, &action_weights[0], output_fp);
printf("returned print stats\n");

print_dashes(80,output_fp);
fprintf(output_fp,"\n\t\t\t\tFinished this pass in %.5f sec\n",timer(*p_start_time));

printf(" returning to main\n");
return finished;
}
		       


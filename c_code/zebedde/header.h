/************************************/
/* THIS IS MY STUFF                 */
/************************************/

#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

/********************** read variables********************/

EXTERNAL int read_new_line;
EXTERNAL int line_no;


/********************** files pointers********************/

EXTERNAL FILE *input_fp; /*input file */
EXTERNAL FILE *output_fp;
EXTERNAL FILE *pore_fp; /* pore coordinate file */
EXTERNAL FILE *seed_fp; /* seed coordinate file */
EXTERNAL FILE *fragment_fp; /* fragment coordinate file */
EXTERNAL FILE *template_fp; /* good template file pointer */
EXTERNAL FILE *strategy_fp; /* strategy file for discover */
EXTERNAL FILE *forcefield_fp; /* forcefield_file for discover */
EXTERNAL FILE *discover_fp; /* discover files (car/mdf) for minimizer */
EXTERNAL FILE *gooduns_fp;  /* accepted templates arc file */

EXTERNAL FILE *gulp_pots_fp; /*potentials file to be used for gulp */
EXTERNAL FILE *extra_gulp_fp; /*file used for extra lines in gulp */
EXTERNAL FILE *pre_opt_fp;    /*file to store template pre_opt in anneal*/


/******* Animation file pointers **********************************/

EXTERNAL FILE *anim_file_fp;
EXTERNAL FILE *anim_read_fp;
EXTERNAL FILE *anim_show_fp;

EXTERNAL FILE *gooduns_read_fp;
EXTERNAL FILE *gooduns_show_fp;

/************ files and temporary file read variables ***************/

EXTERNAL char inputfile[FILELEN_MAX], outputfile[FILELEN_MAX];
EXTERNAL char pore_file[FILELEN_MAX], seed_file[FILELEN_MAX];

EXTERNAL char gulp_pots_file[FILELEN_MAX];

EXTERNAL char fragment_file[FILELEN_MAX];
EXTERNAL char template_file[FILELEN_MAX];
EXTERNAL char template_strategy_file[FILELEN_MAX];
EXTERNAL char  inpore_strategy_file[FILELEN_MAX];
EXTERNAL char anneal_strategy_file[FILELEN_MAX];
EXTERNAL char gooduns[FILELEN_MAX];
EXTERNAL char these_gooduns[FILELEN_MAX];
EXTERNAL char defaults_file[FILELEN_MAX];

EXTERNAL char discover_path[BUFFER];            /* path for discover */
EXTERNAL char template_min_car[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char template_min_mdf[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char inpore_min_car[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char inpore_min_mdf[FILELEN_MAX]; /* for intermediate discover use */

EXTERNAL char discover_forcefield_name[FILELEN_MAX];
EXTERNAL char forcefield_library[FILELEN_MAX]; /*user potentials file */

/****************************************************************************/
/********** Variables used to keep track of potentials **********************/
/****************************************************************************/

EXTERNAL int num_potential_types;
EXTERNAL int num_distinct_pot_types;    /* For double pots with no combining rules ! DJW Feb 99 */
EXTERNAL int num_hbonds;
EXTERNAL int num_bond_incs;             /* For holding bond increments  ! DJW Mar 99 */
EXTERNAL int num_stretches;
EXTERNAL int num_angles;
EXTERNAL int num_torsions;
EXTERNAL int num_allowed_torsions;
EXTERNAL int num_equivalences;

EXTERNAL angle_members angle_warned[MAX_ANGLE_WARN_LIST]; /* hold list of potential sets that have been mentioned */
EXTERNAL int num_angle_warnings;                          /* to the user as unassignable                          */

EXTERNAL int final_gulp;
EXTERNAL int final_annealed_gulp;
EXTERNAL int gulp_input_count;

EXTERNAL char minimizer_name[20];      	/* name of minimization tech to use */
EXTERNAL char title[BUFFER];
EXTERNAL char buffer[BUFFER];
EXTERNAL char dummy_head[BUFFER];   /* temp storage for biosym headers */
EXTERNAL char pore_title[BUFFER];
EXTERNAL char seed_title[BUFFER];
EXTERNAL char fragment_title[BUFFER];


/* command line for mopac from strategy file*/

EXTERNAL char mopac_cmdline_molecule[BUFFER];  
EXTERNAL char mopac_cmdline_inpore[BUFFER]; 
EXTERNAL char gulp_cmdline_molecule[BUFFER];
EXTERNAL char gulp_cmdline_inpore[BUFFER];
EXTERNAL char mopac_root[BUFFER];
EXTERNAL char gulp_root[BUFFER];
EXTERNAL char mopac_path[BUFFER];  /* path for mopac*/
EXTERNAL char gulp_path[BUFFER];

#ifdef MAIN
EXTERNAL int max_mopac_atoms =99999;  /*silly until we work out what the max is*/
#else
EXTERNAL int max_mopac_atoms;  
#endif

/********************* pore/template/fragment structures*********************/

/*EXTERNAL atom pore[MAX_PORE_ATOM];*/
EXTERNAL atom_number atom_limit[NUM_ELEMENTS];/* struct with max elem. concentration*/
EXTERNAL int have_forbidden_bonds;		/* logical for forbidden bonds */
EXTERNAL int num_forbidden_bonds;			/* number of forbidden bonds */
EXTERNAL int have_conc_limits; 			/* logical for concentration limit */
EXTERNAL int num_conc_limits; 			/* number of   concentration limit */
EXTERNAL dihedral allowed_torsions[MAX_ALLOWED_TORS]; /* user defined flexible torsions */
EXTERNAL bond forbidden_bond[MAXTEMPLATE];   /* struct for forbidden bonds */
/*EXTERNAL atom frag_lib[MAXTEMPLATE];*/
EXTERNAL int frag_lib_atno;
EXTERNAL int action_weights[10];  /* weights for action selection */
/*EXTERNAL int frag_weights[MAXTEMPLATE];*/  /* weights for fragment selection */
/*EXTERNAL int seed_weights[MAXTEMPLATE];*/  /* weights for fragment seed selection */
EXTERNAL int num_frag_weights;
EXTERNAL int num_seed_weights;
EXTERNAL int num_action_weights;
EXTERNAL int sum_action_weights;
EXTERNAL int sum_frag_weights;
EXTERNAL int sum_seed_weights;
EXTERNAL int frag_weights_given;
EXTERNAL int seed_weights_given;
EXTERNAL int action_weights_given;
EXTERNAL int num_pore_atoms;
EXTERNAL int total_frag_atoms;      /* total no of frgament atoms */
EXTERNAL int max_templates;
EXTERNAL int number_of_fragments;
EXTERNAL int number_of_members[MAXFRAGMENTS]; /* array of no of atoms per frag */
EXTERNAL int member_start[MAXFRAGMENTS]; /* reference list to start of each frag*/
EXTERNAL double closest_contact[MAXTEMPLATE]; /* contact to pore for each atom*/
EXTERNAL int centralise_template; /*flag: centre seed at c of mass of pore */
EXTERNAL int which_fragment_as_seed; /* which frgament from library to use as seed */

EXTERNAL int optimise_post_mc;      /* sets to optimise after an monte carlo run */
EXTERNAL int want_final_gulp;       /* does the user want a final gulp run? */

EXTERNAL int anneal_last;	    /* whether you are on the final 0K part of annealing */
EXTERNAL int no_annealed_temps;	    

/*************************************************************/
/****** AJWL additions 2008 **********************************/
/*************************************************************/

EXTERNAL double dock_energy;       /*docking rejection energy*/
EXTERNAL int car_output;            
EXTERNAL FILE *docked_template_fp; /*car file with docked template*/
EXTERNAL FILE *docked_zeb_input_fp; /*file for generating new input file with docked template*/
EXTERNAL FILE *condor_submit_fp; /*file for condor input*/
EXTERNAL char temp_car_output[FILELEN_MAX]; /*temp place to store the car file for printing*/
EXTERNAL char condor_submit_filename[FILELEN_MAX]; /*temp place to store condor submit filename*/
EXTERNAL int want_condor;           /* does the user want to run condor? */
EXTERNAL int car_file_count;
EXTERNAL char zebedde_input_filename[FILELEN_MAX];
EXTERNAL char os_linux;
EXTERNAL char os_windows;
EXTERNAL int cmd_line_read; /* Option to allow ZEBEDDE to read pore and seed from command line. Added Exxon 12/9/08*/
EXTERNAL int car_with_pore; /* Print out a car file with the pore. Added Exxon 12/9/08*/
EXTERNAL int exxon; /*Are we running on ExxonMobil cluster*/
EXTERNAL FILE *exxon_submit_fp;
EXTERNAL char exxon_submit_filename[FILELEN_MAX];
EXTERNAL int want_final_car;
EXTERNAL int anneal_run; /* We are running an annealing after dock*/
EXTERNAL int num_dock;
EXTERNAL int total_docked;
EXTERNAL int finished_dock;
EXTERNAL int ucl;
EXTERNAL int told_of_clash;
/*********************periodic structure info ***********************/

/********************************************************************/
/***** pbc    : flags that periodic boundary conditions are *********/
/*****          to be used                                  *********/
/***** abc    : holds a b c alpha beta gamma Angstroms      *********/
/*****                                     and degrees      *********/
/***** latt_vec : holds cartessian vectors for a b c        *********/
/***** recip_latt_vec : holds recip. space a* b* c*         *********/
/*****                                                      *********/
/***** use_xyz, use_zyx : logical to set axis convention    *********/
/********************************************************************/

EXTERNAL int pbc;
EXTERNAL int use_xyz;
EXTERNAL int use_zyx;
EXTERNAL double abc[6];
EXTERNAL double latt_vec[9];
EXTERNAL double real_latt_sizes[3];
EXTERNAL double recip_latt_vec[9];
EXTERNAL double recip_latt_sizes[3];


/*********************stuff to do with symmetry *********************/

EXTERNAL symm_ops symm[MAX_SYMMOPS]; /* the symmetry operations themselves */
EXTERNAL int symm_set;      /* flag to say symmetry set */
EXTERNAL int num_symm_ops;  /* number of symmetry operations defined (0=1!) */

/*********************Box limits for non periodic systems ***********/
EXTERNAL double box_limits[6];
EXTERNAL int user_box;		/* flag for user supplied or not */
EXTERNAL double box_fraction; /* fraction of extents to use */

/********************* Flags for energy and output *********************/
EXTERNAL int initial_minimize_template;
            /*flag: do a minimize template before we start or not */
EXTERNAL int initial_minimize_inpore;
            /*flag: do a minimize template inpore before we start or not */


EXTERNAL int verbose; /* verbose output flag */
EXTERNAL int steric;
EXTERNAL int non_bonded;
EXTERNAL int charges;
EXTERNAL double vdw_scale;
EXTERNAL double stop_ctf; /*stopping bit */
EXTERNAL double ch_ctf;
EXTERNAL double nb_ctf, nb_ctf_2; /* non-bond and charge cutoff distances */
EXTERNAL double hb_ctf, hb_ctf_2; /* hbond cutoff distances for AMBER DJW March 99 */
EXTERNAL double ring_ctf, ring_ctf_2; /* ring formation cutoff distances */

EXTERNAL double max_shake_step;
EXTERNAL double  max_rock_step;
EXTERNAL double  temperature; /* System temperature for MC process */

EXTERNAL double  mc_temp_step;  /*if optimising after MC, decrement temp by this */
EXTERNAL int  mc_modi_step;	/*     and carry out this many actions at that temp */

EXTERNAL double current_temp;      /* temp being used in next stage of annealing opt */

EXTERNAL int  num_modify_attempts;
EXTERNAL int  num_shake_attempts;
EXTERNAL int  num_rock_attempts;

EXTERNAL int animate_flag; /* flag for writing animation */
EXTERNAL int num_anime_frames; /* count of frames in animation */
EXTERNAL int num_goodun_frames; /* count of frames in animation */

EXTERNAL energy interaction_energy[MAX_MOLS];
EXTERNAL internal_energy intra_energy[MAX_MOLS];

/******************** Restraint data ******* DJW July 98 ******/

EXTERNAL int line_restraints;
EXTERNAL int line_hold;
EXTERNAL int num_line_restraints;
EXTERNAL lines line_atom[MAX_LINES];
EXTERNAL tethers tether;
EXTERNAL int have_tethers;

/******************** DEBUG flag ******************************/

EXTERNAL int DEBUG;

/******************** Control data ****************************/

EXTERNAL double prob_test;
EXTERNAL int logfile_needed;

/******************** Allow forced dihedrals, Dave Willock July 98 ****/

EXTERNAL int force_dihedrals;
EXTERNAL int num_for_dihedrals;
EXTERNAL dihedral set_diherals[MAX_SET_DI]; 

/******************** Constrain atoms in c2disco optimisations ********/

EXTERNAL int mini_fix_atoms;
EXTERNAL int num_mini_fix;
EXTERNAL neigh_set fix_atoms[MAXTEMPLATE];

#undef EXTERNAL

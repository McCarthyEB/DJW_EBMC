#define BOOLEAN int
#define TRUE 1
#define FALSE 0

#define NOT_SET "Not_Set"

#define MIN(A,B) (A < B ? A : B)
#define MAX(A,B) (A > B ? A : B)

/******** DEFAULTS WHICH MAY BE INSTALLATION SPECIFIC *****************/
/******** DEFAULTS WHICH MAY BE INSTALLATION SPECIFIC *****************/
/******** DEFAULTS WHICH MAY BE INSTALLATION SPECIFIC *****************/

#define PROB_TEST_DEFAULT 0.01   /* (real) default test probability */
#define MODIFY_TRY_DEFAULT 500   /* (int) default max no of modify attempts*/
#define STOP_CUTOFF_DEFAULT 30.0  /* (real) default stop cutoff distance */
#define RING_CUTOFF_DEFAULT 3.50  /* (real) default ring cutoff distance */
#define CH_CUTOFF_DEFAULT 30.0  /* (real) default charge cutoff distance */
#define NB_CUTOFF_DEFAULT 12.0  /* (real) default non-bond cutoff distance */
#define ATTEMPTS_DEFAULT 1      /* (int) default for shake/rock tries */
#define SHAKE_STEP_DEFAULT 0.1  /* (real) default for shake step */
#define ROCK_STEP_DEFAULT 1.0   /* (real) default for rock step */
#define VDW_SCALE_DEFAULT 1.0   /* (real) default vdw scaling parameter */

#define BOX_FRACTION_DEFAULT 0.8   /* (real) default box fraction parameter */


/*****default filenames for the intermediate discover files******/
#define DISCOVER_FORCEFIELD_DEFAULT "$BIOSYM_LIBRARY/cff91_czeo.bin"
#define TEMPLATE_STRATEGY_DEFAULT "temp_template_min.inp"
#define TEMPLATE_CAR_DEFAULT "temp_template_min.car"
#define TEMPLATE_MDF_DEFAULT "temp_template_min.mdf"
#define INPORE_STRATEGY_DEFAULT "temp_inpore_min.inp"
#define INPORE_CAR_DEFAULT "temp_inpore_min.car"
#define INPORE_MDF_DEFAULT "temp_inpore_min.mdf"
#define GOODUNS_DEFAULT    "good_templates"
#define BOX_OUTPUT_DEFAULT "box"

/*****filename for the peek file******/
#define PEEK_FILENAME "zebedde_peek.cor"

/*****default filenames for analysis run i.e. none *******/
#define NO_ANALYSE "NO_ANALYSE"

/*****default command line for MOPAC*****/
#define DEFAULT_MOPAC_COMMANDLINE "PM3 NOINTER XYZ GEO-OK MMOK "
#define DEFAULT_MOPAC_OUTPUT "mopac_minimize"

/*****default command line for GULP *****/
#define DEFAULT_GULP_COMMANDLINE "opti molq qok"
#define DEFAULT_GULP_OUTPUT "gulp_minimize"

/******* defaults for animations ********************************/

#define DEFAULT_ANIMATION_FILE "animation"
#define DEFAULT_ANIMATION_TYPE BIOSYM_ANIMATION

/****************************************************************/
/***** key words for potential reader ***************************/
/****************************************************************/

#define INFO_LINE            "@"
#define TITLE_LINE           "!"
#define ILLUSTRATION_LINE    ">"
#define JUST_RETURN          "\n"

/****************************************************************/
/***** seed types for seed selection  ***************************/
#define MOLE 1
#define ARCH  2
#define FRAG 3
/***note abbreviated to avoid clash with reader***/

/****************************************************************/
/***** Non-bond potential keywords ******************************/
/****************************************************************/

#define POT_COMBINATION "combination"
#define POT_TYPE        "type"

#define R_EPS           "r-eps"
#define A_B             "A-B"
#define BUCK            "buck"
#define SIXTH_POWER     "sixth-power"
#define ARITHMETIC      "arithmetic"
#define GEOMETRIC       "geometric"
#define NONE            "none"

/**********************************************************/
/****** Potential file headings for potential types *******/
/****** Can now do pcff or cvff Dave Willock Mar.97 *******/
/****** and hopefully AMBER !   Dave Willock Jan.99 *******/
/**********************************************************/

#define VERSION "#version" 
#define PCFF_STRING "pcff"
#define CVFF_STRING "cvff"
#define CFF91_STRING "cff91"
#define AMBER_STRING "amber"
#define OIE_STRING "oie"

#define PCFF 0
#define CVFF 1
#define CFF91 2
#define AMBER 3
#define OIE   4

#define EQUIVALENCE_PCFF "#equivalence"

#define EQUIVALENCE_CVFF "#equivalence"

#define EQUIVALENCE_CFF91 "#equivalence"

#define H_BOND_AMBER "#hydrogen_bond(10-12)"
#define EQUIVALENCE_AMBER "#equivalence"
#define BOND_INCREMENTS "#bond_increments"

/*******************************************/
/*** OIE uses both buck and 12_6 so need ***/
/*** to recongnise potentials by function***/
/*** rather than by author               ***/
/*** Dave Willock Jan 06                 ***/
/*******************************************/

#define NON_BOND_VDW_BUCK   "#nonbond(buck)"
#define NON_BOND_VDW_12_6   "#nonbond(12-6)"
#define NON_BOND_VDW_9_6    "#nonbond(9-6)"

#define H_BOND_OIE   "#hydrogen_bond(10-12)"
#define EQUIVALENCE_OIE   "#equivalence"
#define BOND_INC_OIE   "#bond_increments"

#define FAILED_VDW -1 
#define NB_12_6     0 
#define NB_9_6      1 
#define NB_BUCK     2 

#define FAILED_STRETCH    -1 
#define QUARTIC_STRETCH    0 
#define MORSE_STRETCH      1
#define QUADRATIC_STRETCH  2

#define QUARTIC_STRETCH_STRING    "#quartic_bond"
#define MORSE_STRETCH_STRING      "#morse_bond"
#define QUADRATIC_STRETCH_STRING  "#quadratic_bond"

#define FAILED_ANGLE      -1 
#define QUARTIC_ANGLE      0 
#define QUADRATIC_ANGLE    1

#define QUARTIC_ANGLE_STRING    "#quartic_angle"
#define QUADRATIC_ANGLE_STRING  "#quadratic_angle"

#define FAILED_TORSION    -1 
#define TORSION_1          0 
#define TORSION_3          1

#define TORSION_1_STRING    "#torsion_1"
#define TORSION_3_STRING  "#torsion_3"

/**********************************************************/
/****** own error codes ***********************************/
/**********************************************************/

#define END_OF_FILE -10 /* As defined in read_line */

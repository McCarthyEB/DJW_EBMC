#define MAX_PORE_ATOM 10000
#define MAX_LINE_LEN 256
#define BUFFER  512
#define PORE_GROUP pore
#define END_CAR "end"
#define NUM_POTENTIALS 28
#define MAX_POTS 2000     /* Maximum number of potential types */
#define NUM_ELEMENTS 105  /* update when more found we found deuterium! and lone pairs */

#define MAXTEMPLATE 10000 /* used as max size of template*/
#define MAXTEMPLATE3 30000 /* this needs to be 3xMAXTEMPLATE for gradient vector */
#define MAXFRAGMENTS 200 /* num of diff fragments in library */
#define MAX_SYMMOPS  50 /* maximum number of symmetry operators allowed */
#define MAX_MOLS 500 /* Maximum number of descrete molecules in the template */
#define MAX_MOLS3 1500 /* Maximum number of descrete molecules times 3 to define vectors like c_of_m */
#define FILELEN_MAX 512

/****************************************************************/
/***** Maximum Atom types allowed for fragment library and pore */
/****************************************************************/

#define MAX_TYPES 200

/****************************************************************/
/***** Maximum array dimension for users of ends in searches ****/
/****************************************************************/

#define MAX_ENDS 1000

/****************************************************************/
/***** Maximum list length for pair interaction lists ***********/
/***** PAIR is used for inter-molecular,              ***********/
/***** VDW for intra-molecular van der Waals          ***********/
/****************************************************************/

#define MAX_PAIR_LIST 50000
#define MAX_PAIR_LIST3 150000
#define MAX_ANGLE_LIST 100000
#define MAX_ANGLE_WARN_LIST 300
#define MAX_TORSION_LIST 150000
#define MAX_ALLOWED_TORS 20
#define MAX_VDW_LIST 500000

/****************************************************************/
/***** Maximum lenght of path variable                ***********/
/****************************************************************/

#define ZEB_PATH_MAX 300

/****************************************************************/
/***** Maxima for Ewald sum                           ***********/
/****************************************************************/

#define MAX_KVECS  5000       /* Maximum number of k-vectors in Ewald sum         */
#define MAX_KVECS_ELEM 15000 /* Maximum number of k-vector elements in Ewald sum */
                             /* i.e. 3xMAX_KVECS                                 */

/****************************************************************/
/**** Maximum allowed fixed dihedrals ***************************/
/****************************************************************/

#define MAX_SET_DI 100

/****************************************************************/
/**** Maximum allowed line restraints ***************************/
/****************************************************************/

#define MAX_LINES 10

/****************************************************************/
/**** Maximum close contact average   ***************************/
/****************************************************************/

#define MAX_AVE_CC 1000

/****************************************************************/
/*** Maximum length of hydrogen list in local copy used in ******/
/*** build.c                                               ******/
/****************************************************************/

#define MAX_HYDS 100

/****************************************************************/
/*** Maximum number of atoms in a monitored list           ******/
/*** used in make_a_template.c                             ******/
/****************************************************************/

#define MAX_MONIT 1000

/****************************************************************/
/** Maximum number of neighbours allowed for any atom ***********/
/****************************************************************/
#define MAX_ATOM_NEIGHS 5

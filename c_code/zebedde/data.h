/*************************************************************************/
/* data.h :                                                              */  
/* Changes: 12/3/96 Dewi shortened standard Si-O bond to 1.6 from 1.7    */
/* Changes:  3/7/96 Dewi added     standard Al-O bond at 1.8             */
/* Changes:  5/3/97 Dave added     standard  P-O bond at 1.6             */
/*************************************************************************/
#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

typedef struct
{
  char atom[6];
  char pot[6];
  char pot2[6];
  double a;
  double b;
  double c;
  double sqrt_a;
  double sqrt_b;
  double a3;
  double a6;
  int use_comb;
  int type;
} potential_list;

typedef struct
{
  int iatm1;
  int iatm2;
  int iatm3;
  int iangpot;
} angle_interact_list;

typedef struct
{
  int iatm1;
  int iatm2;
  int ivdwpot;
} vdw_interact_list;

typedef struct
{
  int iatm1;
  int iatm2;
  int iatm3;
  int iatm4;
  int itorpot;
} torsion_interact_list;

typedef struct
{
int  which;
char atom1[6];
char atom2[6];
double stretch_A;
double stretch_B;
double stretch_C;
double stretch_D;
}intra_pair_pot_list;

typedef struct
{
int  which;
char atom1[6];
char atom2[6];
char atom3[6];
double A;
double B;
double C;
double D;
}intra_angle_pot_list;

typedef struct
{
int  which;
char atom1[6];
char atom2[6];
char atom3[6];
char atom4[6];
double A;
double B;
double C;
double D;
double E;
double F;
}intra_torsion_pot_list;

typedef struct
{
  char name[30];
}name_list;

typedef struct
{
char type[6];
char nonbond[6];
char stretch[6];
char angle[6];
char torsion[6];
char oop[6];
}equivalences;

#define EQUIV_NONBOND 0
#define EQUIV_STRETCH 1
#define EQUIV_ANGLE   2
#define EQUIV_TORSION 3
#define EQUIV_OOP     4

typedef struct
{
  char elem[3];
  int  number;
  double mass;
  double covalent;
  double vdw;
} periodic_table_data;    /* for vdw and other information */


typedef struct
{
  char label1[8];
  char label2[8];
  double bond_length;
}bond_dist;

typedef struct
{
  char label[3];
  double value;
}fractions;

/* more complex coord for Template with bits for valencies etc */
/* under development!*/
typedef struct
{
  char name[6];
  double x_coord;
  double y_coord;
  double z_coord;
  int    free_valence;
  int    include;

} complex_coord;

typedef struct
{
  char type[25];
  char combination[25];
  char hbond[25];
  int version;
}potential_info;

/***** If in main fill data arrays **********************************/
#ifdef MAIN

EXTERNAL int num_bond_list = 27;
EXTERNAL bond_dist bond_table[] = {
{"C", "C", 1.51},
{"C", "N", 1.46},
{"C", "Cl", 1.80},
{"N", "N", 1.46},
{"N", "H", 1.15},
{"N", "HA", 1.15},
{"N", "HB", 1.15},
{"N", "L", 1.15},
{"C","H", 1.1},
{"C","HA", 1.1},
{"C","HB", 1.1},
{"O","H", 1.1},
{"O","HA", 1.1},
{"O","HB", 1.1},
{"O","L", 1.1},
{"N", "D", 1.15},
{"C","D", 1.15},
{"O","D", 1.0},
{"Si","O",1.65},
{"Al","O",1.80},
{"P","O",1.60},
{"S","O",1.60},
{"S","L",1.20},
{"S","H",1.40},
{"S","D",1.40},
{"S","C",2.00},
{"C", "O", 1.43},
{"N", "O", 1.21},
};

EXTERNAL int num_fractions = 5;
EXTERNAL fractions fraction_list[] = {
{"1/4", 0.25},
{"3/4", 0.75},
{"1/2", 0.5},
{"1/3", 1.0/3.0},
{"2/3", 2.0/3.0},
};

/*	Element, Z, mass, covalent radius, van der Waals radius*/
/*	Allow HA and HB as synonyms for H                      */

EXTERNAL periodic_table_data period_table[] = {
{ "Ac", 89,  227.028,  1.88,  2.81  },
{ "Ag", 47,  107.868,  1.59,  2.37  },
{ "Al", 13,  26.982,  1.35,  2.01  },
{ "Am", 95,  243.000,  1.51,  2.25  },
{ "Ar", 18,  39.948,  0.00,  0.00  },
{ "As", 33,  74.922,  1.21,  2.20  },
{ "At", 85,  210.000,  0.00,  0.00  },
{ "Au", 79,  196.967,  1.50,  1.85  },
{ "B", 5,  10.811,  0.83,  1.70  },
{ "Ba", 56,  137.327,  1.34,  2.00  },
{ "Be", 4,  9.012,  0.35,  0.52  },
{ "Bi", 83,  208.980,  1.54,  2.30  },
{ "Bk", 97,  247.000,  0.00,  0.00  },
{ "Br", 35,  79.904,  1.21,  1.95  },
{ "C", 6,  12.011,  0.68,  1.75  },
{ "Ca", 20,  40.078,  0.99,  1.48  },
{ "Cd", 48,  112.411,  1.69,  2.52  },
{ "Ce", 58,  140.115,  1.83,  2.73  },
{ "Cf", 98,  251.000,  0.00,  0.00  },
{ "Cl", 17,  35.453,  0.99,  1.77  },
{ "Cm", 96,  247.000,  0.00,  0.00  },
{ "Co", 27,  58.933,  1.33,  1.99  },
{ "Cr", 24,  51.996,  1.35,  2.01  },
{ "Cs", 55,  132.905,  1.67,  2.49  },
{ "Cu", 29,  63.546,  1.52,  1.54  },
{ "Dy", 66,  162.500,  1.75,  2.61  },
{ "Er", 68,  167.260,  1.73,  2.58  },
{ "Es", 99,  252.000,  0.00,  0.00  },
{ "Eu", 63,  151.965,  1.99,  2.97  },
{ "F", 9,  18.998,  0.64,  1.30  },
{ "Fe", 26,  55.847,  1.34,  2.00  },
{ "Fm", 100,  257.000,  0.00,  0.00  },
{ "Fr", 87,  223.000,  0.00,  0.00  },
{ "Ga", 31,  69.723,  1.22,  1.82  },
{ "Gd", 64,  157.250,  1.79,  2.67  },
{ "Ge", 32,  72.610,  1.17,  1.75  },
{ "H", 1,  1.008,  0.23,  1.17  },
{ "HA", 1,  1.008,  0.23,  1.17  },
{ "HB", 1,  1.008,  0.23,  1.17  },
{ "D", 1,  2.000,  0.23,  1.17  },
{ "He", 2,  4.003,  0.00,  0.00  },
{ "Hf", 72,  178.490,  1.57,  2.34  },
{ "Hg", 80,  200.590,  1.70,  1.90  },
{ "Ho", 67,  164.930,  1.74,  2.60  },
{ "I", 53,  126.905,  1.40,  2.10  },
{ "In", 49,  114.820,  1.63,  2.43  },
{ "Ir", 77,  192.220,  1.32,  1.97  },
{ "K", 19,  39.098,  1.33,  1.99  },
{ "Kr", 36,  83.800,  1.89,  2.82  },
{ "La", 57,  138.906,  1.87,  2.79  },
{ "Li", 3,  6.941,  0.68,  1.01  },
{ "Lr", 103,  260.000,  0.00,  0.00  },
{ "Lu", 71,  174.967,  1.72,  2.57  },
{ "Md", 101,  258.000,  0.00,  0.00  },
{ "Mg", 12,  24.305,  1.10,  1.64  },
{ "Mn", 25,  54.938,  1.35,  2.01  },
{ "Mo", 42,  95.940,  1.47,  2.19  },
{ "N", 7,  14.007,  0.68,  1.55  },
{ "Na", 11,  22.990,  0.97,  1.45  },
{ "Nb", 41,  92.906,  1.48,  2.21  },
{ "Nd", 60,  114.240,  1.81,  2.70  },
{ "Ne", 10,  20.180,  0.00,  0.00  },
{ "Ni", 28,  58.690,  1.50,  1.81  },
{ "No", 102,  259.000,  0.00,  0.00  },
{ "Np", 93,  237.048,  1.55,  2.31  },
{ "O", 8,  15.999,  0.68,  1.40  },
{ "Os", 76,  190.200,  1.37,  2.04  },
{ "P", 15,  30.974,  1.05,  1.90  },
{ "Pa", 91,  231.036,  1.61,  2.40  },
{ "Pb", 82,  207.200,  1.54,  2.30  },
{ "Pd", 46,  106.420,  1.50,  2.24  },
{ "Pm", 61,  145.000,  1.80,  2.69  },
{ "Po", 84,  209.000,  1.68,  2.51  },
{ "Pr", 59,  104.908,  1.82,  2.72  },
{ "Pt", 78,  195.080,  1.50,  1.97  },
{ "Pu", 94,  244.000,  1.53,  2.28  },
{ "Ra", 88,  226.025,  1.90,  2.84  },
{ "Rb", 37,  85.468,  1.47,  2.19  },
{ "Re", 75,  186.207,  1.35,  2.01  },
{ "Rh", 45,  102.906,  1.45,  2.16  },
{ "Rn", 86,  222.000,  0.00,  0.00  },
{ "Ru", 44,  101.070,  1.40,  2.09  },
{ "S", 16,  32.066,  1.02,  1.80  },
{ "Sb", 51,  121.750,  1.46,  2.17  },
{ "Sc", 21,  44.956,  1.44,  2.15  },
{ "Se", 34,  78.960,  1.22,  2.00  },
{ "Si", 14,  28.086,  1.20,  2.00  },
{ "Sm", 62,  150.360,  1.80,  2.69  },
{ "Sn", 50,  118.710,  1.46,  2.18  },
{ "Sr", 38,  87.620,  1.12,  1.67  },
{ "Ta", 73,  180.948,  1.43,  2.13  },
{ "Tb", 65,  158.925,  1.76,  2.63  },
{ "Tc", 43,  98.000,  1.35,  2.01  },
{ "Te", 52,  127.600,  1.47,  2.20  },
{ "Th", 90,  232.038,  1.79,  2.67  },
{ "Ti", 22,  47.880,  1.47,  2.19  },
{ "Tl", 81,  204.383,  1.55,  2.31  },
{ "Tm", 69,  168.934,  1.72,  2.57  },
{ "U", 92,  238.029,  1.58,  2.36  },
{ "V", 23,  50.942,  1.33,  1.99  },
{ "W", 74,  183.850,  1.37,  2.04  },
{ "Xe", 54,  131.290,  0.00,  0.00  },
{ "Y", 39,  88.906,  1.78,  2.66  },
{ "Yb", 70,  173.040,  1.94,  2.90  },
{ "Zn", 30,  65.390,  1.45,  2.16  },
{ "Zr", 40,  91.224,  1.56,  2.33  },
{ "L", 0,   0.001,  1.00,  1.00  },
};


/***** If not in main fill only declare data arrays *******/
#else
EXTERNAL int num_bond_list;
EXTERNAL bond_dist bond_table[] ;
EXTERNAL int num_fractions;
EXTERNAL fractions fraction_list[];
EXTERNAL periodic_table_data period_table[];

#endif

/**********************************************************/
/***** Arrays for potential related data ******************/
/**********************************************************/

EXTERNAL potential_list potent[MAX_POTS];
EXTERNAL potential_list h_potent[MAX_POTS];
EXTERNAL potential_list bond_increments[MAX_POTS];
EXTERNAL types h_pot_types[100];
EXTERNAL int num_h_pot_types;
EXTERNAL types bond_inc_types[100];
EXTERNAL int num_bond_inc_types;
EXTERNAL potential_info pot_info ;
EXTERNAL intra_pair_pot_list intra_pair_potent[MAX_POTS];
EXTERNAL intra_angle_pot_list intra_angle_potent[MAX_POTS];
EXTERNAL intra_torsion_pot_list intra_torsion_potent[MAX_POTS];
EXTERNAL equivalences equivalence_list[MAX_POTS];

/**********************************************************/
/****** Some useful #defines to clarify code **************/
/****** Dave Willock started 24th May 1995   **************/
/**********************************************************/

/**********************************************************/
/****** Cost_functions that can be used *******************/
/**********************************************************/

enum { STERIC_COST=1, NON_BOND_COST, ENERGY_COST };

/**********************************************************/
/****** Main Code: actions ********************************/
/**********************************************************/

#define NUMBER_OF_ACTIONS 10

#define BUILD_ACTION              0
#define ROTATE_LAST_FRAG_ACTION   1
#define SHAKE_ACTION              2
#define ROCK_ACTION               3
#define TWIST_ACTION              4
#define RING_MAKER_ACTION         5
#define MINIMIZE_ACTION           6
#define MINIMIZE_INPORE_ACTION    7
#define CREATE_GUEST_ACTION       8
#define REMOVE_GUEST_ACTION       9

#ifdef MAIN

EXTERNAL name_list action_list[NUMBER_OF_ACTIONS] = {
{"Build"},
{"Rotate Last Fragment"},
{"Shake"},
{"Rock"},
{"Twist"},
{"Ring Maker"},
{"Minimise template"},
{"Minimise in pore"},
};

#else
EXTERNAL name_list action_list[NUMBER_OF_ACTIONS];
#endif

/**********************************************************/
/****** Animation file types ******************************/
/**********************************************************/

#define BIOSYM_ANIMATION 1
#define XMOL_ANIMATION_NOPORE   2
#define XMOL_ANIMATION_WITHPORE   3


/**********************************************************/
/****** Minimizer names ***********************************/
/**********************************************************/

#define DISCOVER_MINIMIZER "discover"
#define C2DISCOVER_MINIMIZER "c2discover"
#define MOPAC_MINIMIZER "mopac"
#define INTERNAL_MINIMIZER "internal"
#define MINIMIZER_DEFAULT "discover"
#define GULP_MINIMIZER "gulp"

#define DEFAULT_MOPAC_PATH "/usr/biosym/300/irix5r4/exe/insight/mopac"
#define DEFAULT_DISCOVER_PATH "/usr/biosym/300/biosym_bin/discover"
#define DEFAULT_GULP_PATH "~/mfi/gulp1.4"

#define ZEBEDDE_DEFAULTS "zebedde.def" /*name for defaults file */
#define FINAL_PATH "/progs/zebedde" /*last place to look for zebedde files */

#undef EXTERNAL

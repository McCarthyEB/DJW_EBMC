/* structres.h */
/* Mine first then those from orig (tim) after for reference */

typedef struct
{
 int steric;
 double non_bonded;			        /* non_bonding energy */
 double vdw_rep;                                /* van der Waals repulsive */
 double vdw_disp;                               /* van der Waals dispersive */
 double hbond;                                  /* hbond in AMBER etc */
 double charges;				/* coulombic energy */
 double restraint;                              /* restraint energy DJW July 98 */
 double guest_guest;                            /* part of total energy due to guest-guest interactions */
 double guest_host;                             /* part of total energy due to guest-host  interactions */
 int acceptance;
 double minimizer_init_total;			/* minimizer non_bonding energy */
 double minimizer_end_total;			/* minimizer total energy */
 double minimizer_init_nonbond;			/* minimizer non_bonding energy */
 double minimizer_end_nonbond;			/* minimizer total energy */
} energy;

/**** Structure for molecule-molecule interaction list ***/
typedef struct
{
  int imol;
  int jmol;
  int ind;
  int jnd;
} interaction_indices;

/* rudimentary statistics structure */
 typedef struct
 {
  int tries;
  int tries_after_build;
  int accepted;
  int accepted_after_build;
 }stats;

typedef struct
 {
   char host[7];                    /* host atom type to monitor                    */
   char guest[7];                   /* guest atom type to monitor                   */
   int host_list[MAX_MONIT];        /* list of host monitored sites                 */
   int guest_list[MAX_MONIT];       /* list of guest monitored sites                */
   int entry[MAX_MONIT][MAX_MONIT]; /* The actual data for monitored pairs of atoms */
 } monit;

typedef struct
{
double stretch;
double angle;
double torsion;
double vdw;
double charge;
double total;
} internal_energy;

 
/* structure for the atoms */
typedef struct
{
  char  label[7];
  double x;
  double y;
  double z;
  char pot[6];
  char group[5];
  char group_no[9];
  char elem[3];
  int mol;                      /* index saying which molecule this belongs to */
  double  part_chge; 
  int nb_list;                  /* index of non-bonding parameters */
  int hb_list;                  /* index of h-bonding parameters */
  int bi_list;                  /* index of bond increment list */
  double vdw;			/* van der waals radius of the atom */
  int num_neigh; 		/* number of neighbours */
  int neighb[MAX_ATOM_NEIGHS];  /*initialise -1=not bonded */
  int neighb_stretch_list[MAX_ATOM_NEIGHS];   /* index for intra stretch potential */
  double vdw_energy;            /* This atoms vdw energy contribution to the total */
  double electrostatic_pot;     /* The electrostatic potential at the atom position */
}atom;

typedef struct
{
  char atom1[5];
  char atom2[5];
} bond;

typedef struct
{
  char atom_type[5];
  int  num;
} atom_number;

typedef struct
{
  int start;
  int end;
} links;

typedef struct
{
double matrix[9];
double translation[3];
} symm_ops;

typedef struct
{
int start;
int end;
int num;
} list_partition;

typedef struct
{
char name[6];
} types;

typedef struct
{
char pot[6];
} pot_types;

typedef struct
{
   int mole1;
   char label1[7];
   int index1;
   int mole2;
   char label2[7];
   int index2;
} pair_list;

typedef struct
{
char A[6];
char B[6];
char C[6];
char D[6];
double phi;
}dihedral;

typedef struct
{
char atom1[6];
char atom2[6];
char atom3[6];
int found;
}angle_members;

typedef struct
{
char A[3];
char B[3];
double ref_x;
double ref_y;
double ref_z;
double k;
}lines;

typedef struct
{
char A[5];
char B[5];
int index_a;
int index_b;
double r0;
double k;
}tethers;

typedef struct
{
char central_elem[5];
int num_neigh;
/*** neigh indexed [elements of string][string index] ***/
char neigh[5][9];
}neigh_set;

typedef struct
{
  double av_inter;
  double av_total;
  double std_inter;
  double std_total;
  double acceptance_ratio;
} mc_record;

typedef struct
{
  double v[3];
} vec;

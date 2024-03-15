/************************************************************/
/* align_bonds                                              */
/* Do the geometry re-organisation for a build              */
/*                                                          */
/* Started DJW 14/7/95                                      */
/************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "constants.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

double size_vector(double *p_vector);

double standard_bond( char *atom1, char *atom2 );

void unit_vector(double *p_vector, double *p_size);

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms, int which_mol);

void set_dihedral(atom *p_molecule, atom *p_A, atom *p_B, atom *p_C, atom *p_D,
                  atom *p_new_bit, int num_new_atoms, int new_C, double phi);

void align_bonds(atom *p_template, atom *p_new_bit, atom *p_new_bonding_atom,
                 atom *p_new_hyd_picked, atom *p_temp_bonding_atom,
                 atom *p_temp_hyd_picked, int num_template_atoms,
                 int num_new_atoms)
{
#include "header.h"
double inter_atom_vec[3], new_bond_to_hyd[3], temp_bond_to_hyd[3];
double norm[3], origin[3];
double mag_new_bth, mag_temp_bth, arg, theta, standard;
double dot, size;

double offline_vec[3];
atom offline_atom, *p_neigh;
int bother_to_rotate;
int ineigh, neigh_index;
double aligned_dot;

int i_di, have_BC, have_ABC, have_ABCD;
int new_C=FALSE, ABCD_check= TRUE;

atom *p_A, *p_B, *p_C, *p_D;

/*******************************************************************************/
/*** Set all pointers for dihedral to the template start to test later we have */
/*** a valid call.                                                             */
/*******************************************************************************/

p_A= p_template;
p_B= p_template;
p_C= p_template;
p_D= p_template;

/********************************************************************************/
/*************Position the new bonding atom at the template's********************/
/*************hydrogen that is about to be lost**********************************/
/********************************************************************************/

join_atoms( p_new_bonding_atom, p_temp_hyd_picked, &inter_atom_vec[0]);
move_molecule( p_new_bit, num_new_atoms, &inter_atom_vec[0], -1);

/********************************************************************************/
/**********************find the vectors from the bonding atoms to their *********/
/**********************respective hydrogens**************************************/
/********************************************************************************/

join_atoms( p_new_bonding_atom, p_new_hyd_picked, &new_bond_to_hyd[0]);  
join_atoms( p_temp_hyd_picked, p_temp_bonding_atom, &temp_bond_to_hyd[0]);

mag_new_bth = size_vector( &new_bond_to_hyd[0]);
mag_temp_bth = size_vector( &temp_bond_to_hyd[0]);

/************************* now reorientate **************************************/

vec_cross( &new_bond_to_hyd[0],  &temp_bond_to_hyd[0], &norm[0]); 

bother_to_rotate = TRUE;

if (size_vector(&norm[0]) < 10e-6)
	{
	/**** Template and New bit already aligned along new bond *****/
	/**** Special case because unit_vector divide by zero     *****/

	bother_to_rotate = FALSE;

	aligned_dot = vec_dot(&new_bond_to_hyd[0],  &temp_bond_to_hyd[0]);
        
	if (aligned_dot < 0) 
	   {
	   bother_to_rotate = TRUE;
           /*** anti - aligned: Looks like and SN2 attack!!! ****/
           /* offline_atom is a dummy atom off the line of the SN2 attack! */
           offline_atom.x   = p_new_bonding_atom->x;
           offline_atom.y   = p_new_bonding_atom->y;
           offline_atom.z   = p_new_bonding_atom->z + 1.0;

	   join_atoms( p_new_bonding_atom, &offline_atom, &offline_vec[0]);
           vec_cross(&offline_vec[0], &temp_bond_to_hyd[0], &norm[0]);
           if (size_vector(&norm[0]) < 10e-6)
		{
		/** still aligned!! bugger. so add 1 to y **/
                offline_atom.x   = p_new_bonding_atom->x;
                offline_atom.y   = p_new_bonding_atom->y +1.0;
                offline_atom.z   = p_new_bonding_atom->z;
                join_atoms( p_new_bonding_atom, &offline_atom, &offline_vec[0]);
                vec_cross(&offline_vec[0], &temp_bond_to_hyd[0], &norm[0]);

		}

	   }
	
	}

if (bother_to_rotate)
       {
	  unit_vector(&norm[0], &size);

	  dot= vec_dot( &new_bond_to_hyd[0],  &temp_bond_to_hyd[0]);
	  
	  arg =  ( dot/(mag_new_bth*mag_temp_bth));
	  
	  if (arg > 1.0 && arg < 1.01) arg = 1;
	  if (arg < -1.0 && arg > -1.01) arg = -1;
	  if (arg > 1.01 || arg < -1.01) { printf("Angle out of range\n"); exit(1);}
	  theta= acos ( arg );
	  
	  origin[0]=  p_new_bonding_atom->x;
	  origin[1]=  p_new_bonding_atom->y;
	  origin[2]=  p_new_bonding_atom->z;
	  
	  rotate(p_new_bit, &norm[0], &origin[0], theta, num_new_atoms, -1 );
	  
        }

/************************************************************************************/
/************************* adjust geom to correct bond length ***********************/
/************************************************************************************/

unit_vector(&temp_bond_to_hyd[0], &size);

/* find correct bond length for this pair */
 
standard= standard_bond( &(p_new_bonding_atom->elem[0]),
                         &(p_temp_bonding_atom->elem[0]));


/************************************************************************************/
/**************** Warn user if a bond not in the data.h list is found ***************/
/************************************************************************************/

if (standard < 0)
  {
     printf("Warning: Bonds between %s and %s do not have a standard value in data.h\n",
             p_new_bonding_atom->elem, p_temp_bonding_atom->elem);
     printf("A value of 1.5 Angstroms has been used.\n\n");
     standard= 1.5;
  }

standard -= mag_temp_bth;
temp_bond_to_hyd[0] =  -standard*temp_bond_to_hyd[0];
temp_bond_to_hyd[1] =  -standard*temp_bond_to_hyd[1];
temp_bond_to_hyd[2] =  -standard*temp_bond_to_hyd[2];

move_molecule( p_new_bit, num_new_atoms, &temp_bond_to_hyd[0], -1);

/*******************************************************************/
/**** Enforce dihedral bond for the new bond if required ***********/
/**** The dihedrals are always referenced in terms of    ***********/
/**** atoms A-B-C-D with A-B in the template and C-D     ***********/
/**** in the new bit.                                    ***********/
/**** New section started July 98 Dave Willock           ***********/
/*******************************************************************/

if (force_dihedrals)
  {
/*-----------------------------------------------------------------*/
/***** Loop over the list of dihedral settings given ***************/
/*-----------------------------------------------------------------*/

    i_di=-1;
    have_ABCD= FALSE;

    while (!have_ABCD && i_di < num_for_dihedrals)
      {

         i_di++;
/*-----------------------------------------------------------------*/
/***** Test if the labels are already the right way round **********/
/*-----------------------------------------------------------------*/
         have_BC = FALSE;
         new_C= FALSE;

         if ( strcmp(p_temp_bonding_atom->elem, set_diherals[i_di].B) == 0)
            {
               p_B = p_temp_bonding_atom;
               if ( strcmp(p_new_bonding_atom->elem, set_diherals[i_di].C) == 0)
                  {
                    p_C = p_new_bonding_atom;
                    have_BC = TRUE;
                    new_C=TRUE;
                  }
            }

/*-----------------------------------------------------------------*/
/***** Test for bond reversed w.r.t. standard **********************/
/***** Note that if they are we swap round here ********************/
/*-----------------------------------------------------------------*/

         if ( strcmp(p_temp_bonding_atom->elem, set_diherals[i_di].C) == 0 )
            {
               p_C = p_temp_bonding_atom;
               if ( strcmp(p_new_bonding_atom->elem, set_diherals[i_di].B) == 0 )
                  {
                    p_B = p_new_bonding_atom;
                    have_BC = TRUE;
                  }
            }

/*-----------------------------------------------------------------*/
/***** Get A and D from neighbour lists ****************************/
/*-----------------------------------------------------------------*/

         if (have_BC)
           {
              have_ABC= FALSE;
              for (ineigh = 0; ineigh <= p_B->num_neigh; ineigh++)
                {
                   neigh_index= p_B->neighb[ineigh];
                   p_neigh= p_template+neigh_index;

                   if (strcmp(p_neigh->elem, &(set_diherals[i_di].A[0]) ) == 0)
                     {
                        p_A= p_neigh;
                        have_ABC= TRUE;
                     }
                }
              if (have_ABC)
                {
                   for (ineigh = 0; ineigh <= p_C->num_neigh; ineigh++)
                     {
                        neigh_index= p_C->neighb[ineigh];
                        p_neigh= p_template+neigh_index;

                        if (strcmp(p_neigh->elem, &(set_diherals[i_di].D[0]) ) == 0)
                          {
                             p_D= p_neigh;
                             have_ABCD= TRUE;
                          }
                     }
                }
           }
      }

/*-----------------------------------------------------------------*/
/***** Set the dihedral if one was found ***************************/
/***** Check that the dihedral is at least reasonable with *********/
/***** all atoms different, added Dec 2006 Dave Willock    *********/
/*-----------------------------------------------------------------*/

    if (have_ABCD)
      {
         ABCD_check = p_A != p_B && p_A != p_C && p_A != p_D;
         ABCD_check = ABCD_check && p_B != p_C && p_B != p_D;
         ABCD_check = ABCD_check && p_C != p_D;
 
         if (ABCD_check)
           {
             set_dihedral( p_template, p_A, p_B, p_C, p_D, 
                       p_new_bit, num_new_atoms, new_C, set_diherals[i_di].phi);
           }
         else
           {
             printf("ERROR >> At least two atoms are identical when setting a dihedral\n");
             printf("         This is a program bug that should be notified to the    \n");
             printf("         authors.                                                \n");
             exit(0);
           }
      }
  }

return; 
}

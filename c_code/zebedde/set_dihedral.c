/************************************************************/
/* set_dihedral                                             */
/* Routine to force a user defined dihedral on a newly      */
/* created bond.                                            */
/* Expects to be passed the molecule and pointers to the    */
/* dihedral atoms in the form A-B-C-D if C is in the new    */
/* part of the molecule (i.e. the bit to be rotated )       */
/* the integer new_C will be set TRUE                       */
/*                                                          */
/* Started DJW July 98                                      */
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

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);

void unit_vector(double *p_vector, double *p_size);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms);

void set_dihedral(atom *p_molecule, atom *p_A, atom *p_B, atom *p_C, atom *p_D,
                  atom *p_new_bit, int num_new_atoms, int new_C, double phi)
{
double origin[3], BC[3], BA[3], CD[3];
double norm1[3], norm2[3];
double dot,size;
double phinow, phi_change;

/***** Generate vectors in the ABC and BCD planes to get plane normals ****/

join_atoms( p_B, p_C, &BC[0]);
join_atoms( p_B, p_A, &BA[0]);
join_atoms( p_C, p_D, &CD[0]);

unit_vector(&BC[0], &size);
unit_vector(&BA[0], &size);
unit_vector(&CD[0], &size);

vec_cross(&BA[0], &BC[0], &norm1[0]);
vec_cross(&CD[0], &BC[0], &norm2[0]);

unit_vector(&norm1[0], &size);
unit_vector(&norm2[0], &size);

/***** Work out current angle from dot product *****************************/

dot= vec_dot( &norm1[0],  &norm2[0]);

phinow= acos ( dot );

/***** Work out if this the BCD plane is above or below the ABC ************/

dot= vec_dot( &BA[0], &norm2[0]);

if (dot < 0) phinow= -phinow;

phi_change= phi - phinow;

if (new_C)
  {
    origin[0]=  p_C->x;
    origin[1]=  p_C->y;
    origin[2]=  p_C->z;

  }
else
  {
    origin[0]=  p_B->x;
    origin[1]=  p_B->y;
    origin[2]=  p_B->z;

    BC[0]= -BC[0];
    BC[1]= -BC[1];
    BC[2]= -BC[2];
  }

rotate(p_new_bit, &BC[0], &origin[0], phi_change, num_new_atoms );

return;
}


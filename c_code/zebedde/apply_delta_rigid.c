/**********************************************************/
/***** Routine to apply a rigid body translation **********/
/***** and rotation for the minimiser! DJW Nov 96**********/
/**********************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms, int which_mol);

void unit_vector(double *p_vector, double *p_size);

void apply_delta_rigid(atom *p_molecule, int num_atoms, double *p_c_of_m, double *p_delta)
{
double theta;

printf("\nDelta passed to apply_delta_rigid: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               *p_delta, *(p_delta+1), *(p_delta+2), *(p_delta+3), *(p_delta+4), *(p_delta+5));
printf("Centre of mass passed to apply_delta_rigid: %10.6f %10.6f %10.6f\n",
                                            *p_c_of_m, *(p_c_of_m+1), *(p_c_of_m+2));

unit_vector(p_delta+3, &theta); 

printf("\nAfter unit_vector delta=  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               *p_delta, *(p_delta+1), *(p_delta+2), *(p_delta+3), *(p_delta+4), *(p_delta+5));
printf("Centre of mass after unit_vector: %10.6f %10.6f %10.6f\n",
                                            *p_c_of_m, *(p_c_of_m+1), *(p_c_of_m+2));

rotate(p_molecule, p_delta+3, p_c_of_m, theta, num_atoms, -1); 

printf("\nAfter rotate delta=  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               *p_delta, *(p_delta+1), *(p_delta+2), *(p_delta+3), *(p_delta+4), *(p_delta+5));
printf("Centre of mass after rotate: %10.6f %10.6f %10.6f\n",
                                            *p_c_of_m, *(p_c_of_m+1), *(p_c_of_m+2));

move_molecule(p_molecule, num_atoms, p_delta, -1);

printf("\nAfter move delta=  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               *p_delta, *(p_delta+1), *(p_delta+2), *(p_delta+3), *(p_delta+4), *(p_delta+5));
printf("Centre of mass after move: %10.6f %10.6f %10.6f\n",
                                            *p_c_of_m, *(p_c_of_m+1), *(p_c_of_m+2));

*p_c_of_m     += *p_delta;
*(p_c_of_m+1) += *(p_delta+1);
*(p_c_of_m+2) += *(p_delta+2);

printf("Centre of mass after its move: %10.6f %10.6f %10.6f\n",
                                            *p_c_of_m, *(p_c_of_m+1), *(p_c_of_m+2));

return;
}

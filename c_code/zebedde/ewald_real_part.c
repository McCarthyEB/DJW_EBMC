/****************************************************************************************/
/***** Real part of Ewald Sum ***********************************************************/
/***** Started 12 June 1996 Dave Willock ************************************************/
/****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include "maxima.h"
#include "structures.h"
#include "own_maths.h"
#include "ewald.h"

void min_image( double *x, double *y, double *z);

double own_erf(double x, double accuracy);

double ewald_real_part(atom *p_mol_atom, atom *p_crystal_cell, int num_cell_atoms, int self)
{
#include "header.h"
int crystal_index;
int num_needed[3];
int index_a, index_b, index_c;

double dx_orig, dy_orig, dz_orig;
double dx,dy,dz,potential_contrib, erf_arg;
double dxa, dya, dza, dxb, dyb, dzb;
double this_crys_bit;
double mag_vec, self_factor;

atom *p_crystal_atom;

potential_contrib=0;

/**************************************************************************************/
/****** work out real space sum limits ************************************************/
/**************************************************************************************/

num_needed[0] = 1+real_sum_max/real_latt_sizes[0];
num_needed[1] = 1+real_sum_max/real_latt_sizes[1];
num_needed[2] = 1+real_sum_max/real_latt_sizes[2];

/**************************************************************************************/
/****** Loop over the atom group producing the potential ******************************/
/**************************************************************************************/

if (self)
  {
    self_factor= kappa*erf_normalise;
    p_crystal_atom= p_crystal_cell-1;
    for (crystal_index=0; crystal_index <= num_cell_atoms; crystal_index++)
      {
        this_crys_bit=0;

        p_crystal_atom++;
        dx_orig = p_mol_atom->x - p_crystal_atom->x;
        dy_orig = p_mol_atom->y - p_crystal_atom->y;
        dz_orig = p_mol_atom->z - p_crystal_atom->z;

        min_image( &dx_orig, &dy_orig, &dz_orig);

/***************************************************************************************/
/****** Loop over lattice translations *************************************************/
/***************************************************************************************/

        for (index_a= -num_needed[0]; index_a <= num_needed[0]; index_a++)
          {
            dxa = dx_orig + index_a * latt_vec[0];
            dya = dy_orig + index_a * latt_vec[1];
            dza = dz_orig + index_a * latt_vec[2];

            for ( index_b= -num_needed[1]; index_b <= num_needed[1]; index_b++)
              {
                dxb = dxa+ index_b * latt_vec[3];
                dyb = dya+ index_b * latt_vec[4];
                dzb = dza+ index_b * latt_vec[5];

                for ( index_c= -num_needed[2]; index_c <= num_needed[2]; index_c++)
                  {
                    dx = dxb+ index_c * latt_vec[6];
                    dy = dyb+ index_c * latt_vec[7];
                    dz = dzb+ index_c * latt_vec[8];

                    mag_vec =  dx*dx + dy*dy + dz*dz ;  
         
/*********************************************************************************/
/****** Self interaction terms now need to avoid double counting *****************/
/*********************************************************************************/

                    if (mag_vec <= 1E-6)
                       {
                          this_crys_bit -= self_factor;
                       }

                    else if (mag_vec <= real_sum_max2)
                       {
                          mag_vec= sqrt(mag_vec);
                          erf_arg= kappa*mag_vec;
                          this_crys_bit += (1.0 - own_erf(erf_arg, erf_accuracy))/mag_vec; 
                       }
                 }
             }
         }
       potential_contrib += (p_crystal_atom->part_chge)*this_crys_bit;
     }
   potential_contrib = 0.5*potential_contrib; 
  }

/****** Not self cases ******/
else
  {
    p_crystal_atom= p_crystal_cell-1;
    for (crystal_index=0; crystal_index <= num_cell_atoms; crystal_index++)
      {
        this_crys_bit=0;
        p_crystal_atom++;
        dx_orig = p_mol_atom->x - p_crystal_atom->x;
        dy_orig = p_mol_atom->y - p_crystal_atom->y;
        dz_orig = p_mol_atom->z - p_crystal_atom->z;

        min_image( &dx_orig, &dy_orig, &dz_orig);

/***************************************************************************************/
/****** Loop over lattice translations *************************************************/
/***************************************************************************************/

        for (index_a= -num_needed[0]; index_a <= num_needed[0]; index_a++)
          {
            dxa = dx_orig + index_a * latt_vec[0];
            dya = dy_orig + index_a * latt_vec[1];
            dza = dz_orig + index_a * latt_vec[2];

            for ( index_b= -num_needed[1]; index_b <= num_needed[1]; index_b++)
              {
                dxb = dxa+ index_b * latt_vec[3];
                dyb = dya+ index_b * latt_vec[4];
                dzb = dza+ index_b * latt_vec[5];

                for ( index_c= -num_needed[2]; index_c <= num_needed[2]; index_c++)
                  {
                    dx = dxb+ index_c * latt_vec[6];
                    dy = dyb+ index_c * latt_vec[7];
                    dz = dzb+ index_c * latt_vec[8];

                    mag_vec =  dx*dx + dy*dy + dz*dz ;  
         
                    if (mag_vec <= 1E-6) 
                       {
                         printf ("Problem with zero vector in real space\n\n");
                         exit(EXIT_FAILURE);
                       }

/*********************************************************************************/
/****** Normal term for distinct atoms *******************************************/
/*********************************************************************************/

                    else if (mag_vec <= real_sum_max2)
                       {
                          mag_vec= sqrt(mag_vec);
                          erf_arg= kappa*mag_vec;
                          this_crys_bit += (1.0 - own_erf(erf_arg, erf_accuracy))/mag_vec; 
                       }
                 }
             }
         }
       potential_contrib += (p_crystal_atom->part_chge)*this_crys_bit;
     }
  }

/******** Doing full sum at the moment!! *****/

return potential_contrib;
}

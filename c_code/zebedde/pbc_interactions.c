#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

/*********************************************************************************/
/********  Routine to work out the separations for periodic   ********************/
/********  images of a given vector which are within cut off  ********************/
/********                                                     ********************/
/********  Dave Willock April 1996                            ********************/
/*********************************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *frac_a, double *frac_b, double *frac_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void min_image( double *x, double *y, double *z);

int pbc_interactions(double *p_dx, double *p_dy, double *p_dz, double cutoff,
                                   double cutoff_2, double *p_pos_separations,  
                                   double *p_pos_vectors, int count_flag,
                                   int is_self)

  {
#include "header.h"
double fract_a, fract_b, fract_c;
double base_fract_a, base_fract_b, base_fract_c;
double x,y,z;
double separation;
int num_interactions;
int max_a,max_b,max_c;
int a_vecs, b_vecs, c_vecs;

num_interactions= -1;

/******* Check that the min_image vector is within the cutoff ******/
/******* if not there can be no interaction within the cutoff ******/

min_image(p_dx, p_dy, p_dz);
separation= (*p_dx)*(*p_dx) + (*p_dy)*(*p_dy) + (*p_dz)*(*p_dz);

//printf("Min image separation: %10.6f\n", sqrt(separation));

if (separation <= cutoff_2) 
  {
    if(!is_self && count_flag==0)
    {
      *p_pos_separations= sqrt(separation);

/*******************************************************************/
/**** Remember inter-atomic vector too for gradient calculation ****/
/*******************************************************************/

      *p_pos_vectors= *p_dx;
      p_pos_vectors++;
      *p_pos_vectors= *p_dy;
      p_pos_vectors++;
      *p_pos_vectors= *p_dz;
      p_pos_vectors++;
  
      p_pos_separations++;
      num_interactions++;
    }
  }
else
  {
    return num_interactions;
  }

/*** Only build periodic images lists if they are needed ***/
if ( real_latt_sizes[0] < 2.0*cutoff || real_latt_sizes[1] < 2.0*cutoff || real_latt_sizes[2] < 2.0*cutoff)  
  {
/******* find the fractional co-ords of this the minimal vector ****/
 
cart_to_fract( *p_dx, *p_dy, *p_dz, &base_fract_a, &base_fract_b, &base_fract_c, 
                                                       &recip_latt_vec[0]);

/******* Use reciprocal space vectors to set the search limits ****/

max_a= recip_latt_sizes[0]*cutoff+1;
max_b= recip_latt_sizes[1]*cutoff+1;
max_c= recip_latt_sizes[2]*cutoff+1;
 
/****** Now assemble the list      ******************************/

  for (a_vecs= -max_a; a_vecs < max_a; a_vecs++)
    {
      fract_a= base_fract_a+ a_vecs;

      for (b_vecs= -max_b; b_vecs < max_b; b_vecs++)
        {
          fract_b= base_fract_b+ b_vecs;

          for (c_vecs= -max_c; c_vecs < max_c; c_vecs++)
            {
               if (!(a_vecs==0 && b_vecs==0 && c_vecs==0))
                 {
                    fract_c= base_fract_c+ c_vecs;

/****** convert back to cartessian  *****************************/
/****** and measure/check real size *****************************/

                    fract_to_cart( &x, &y, &z, fract_a, fract_b, fract_c, 
                                                                 &latt_vec[0] );
                    separation= x*x + y*y + z*z;
              
                    if (separation < cutoff_2)
                       {
                          if (num_interactions < MAX_PAIR_LIST)
                             {
                                if(count_flag==0)
                                {
                                  *p_pos_separations= sqrt(separation);

/*******************************************************************/
/**** Remember inter-atomic vector too for gradient calculation ****/
/*******************************************************************/

                                  *p_pos_vectors= *p_dx;
                                   p_pos_vectors++;
                                  *p_pos_vectors= *p_dy;
                                   p_pos_vectors++;
                                  *p_pos_vectors= *p_dz;
                                   p_pos_vectors++;

                                   p_pos_separations++;
                                }
                             }
                          num_interactions++;
                       }
                 }
            }
        }
    }
 }

/*******************************************************************/
/**** Safe guard added 4th July 1998 Dave Willock ******************/
/*******************************************************************/

if (num_interactions > MAX_PAIR_LIST)
   {
      printf("ERROR : The number of pair interactions required (%d) for the vdw interactions ", num_interactions);
      printf("        list exceeds the maximum allowed: %d\n", MAX_PAIR_LIST);
      printf("Either decrease the non-bond cutoff or re-parameterise the program\n");
      exit(0);
   }

return num_interactions;   
}

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

/*********************************************************************************/
/********  Routine to work out the minimum image vector for   ********************/
/********  periodic systems                                   ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *frac_a, double *frac_b, double *frac_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void min_image( double *x, double *y, double *z)
  {
#include "header.h"
double fract_a, fract_b, fract_c;
int done;

/******* find the fractional co-ords of this vector *************/
 
    cart_to_fract( *x, *y, *z, &fract_a, &fract_b, &fract_c, &recip_latt_vec[0]);

/******* adjust to min image ************************************/

    done= FALSE;
    while (!done)
      {
        done = TRUE;
        if (fract_a >= 0.5) { fract_a= fract_a- 1.0; done = FALSE; }
        if (fract_a < -0.5) { fract_a= fract_a+ 1.0; done = FALSE; }

        if (fract_b >= 0.5) { fract_b= fract_b- 1.0; done = FALSE; }
        if (fract_b < -0.5) { fract_b= fract_b+ 1.0; done = FALSE; }

        if (fract_c >= 0.5) { fract_c= fract_c- 1.0; done = FALSE; }
        if (fract_c < -0.5) { fract_c= fract_c+ 1.0; done = FALSE; }
      }

/****** convert back to cartessian ******************************/

    fract_to_cart( x, y, z, fract_a, fract_b, fract_c, &latt_vec[0] );

    return;   

  }

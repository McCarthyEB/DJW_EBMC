/*********************************************************************************/
/********  Routine to work out cartessian real and reciprocal ********************/
/********  lattice vectors from the a,b,c,alpha,beta,gamma    ********************/
/********  line of a biosym .car vector, assuming that the XYZ********************/
/********  convention is followed in the .mdf file.           ********************/
/********                                                     ********************/
/********  p_latt_vec and p_recip_latt_vec point to 9 comp.   ********************/
/********  arrays to be filled                                ********************/
/********  ax ay az bx by bz cx cy cz                         ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "own_maths.h"
#include "ewald.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void cart_latt_vecs( double *p_abc, double *p_latt_vec, double *p_real_latt_sizes, 
                     double *p_recip_latt_vec, double *p_recip_latt_sizes,
                     double *p_cell_volume, int num_atoms, int use_xyz, int use_zyx)
  {
    double bx_cx, az_bz, mag_a, mag_b, mag_c;
    double a_cross_b[3], b_cross_c[3], c_cross_a[3]; 
    double alpha, beta, gamma;
    double alpha_check, beta_check, gamma_check;
    double dot_set[9];
    double cell_length;

    int iloop;

/******** report convention to be used *********************************/

   if (use_xyz) printf("Assuming XYZ, i.e. BIOSYM convention\n");
   if (use_zyx) printf("Assuming ZYX, i.e. MSI convention\n");

/********* make alpha beta gamma in radians ****************************/

    mag_a = *p_abc;
    mag_b = *(p_abc+1);
    mag_c = *(p_abc+2);
    alpha = *(p_abc+3)/RAD_TO_DEG;
    beta  = *(p_abc+4)/RAD_TO_DEG;
    gamma = *(p_abc+5)/RAD_TO_DEG;

    if (use_xyz)
       {
/********* X indicates that a is to be aligned with the x axis *********/

          *p_latt_vec    = mag_a;
          *(p_latt_vec+1)= 0.0;
          *(p_latt_vec+2)= 0.0;

/********* XY indicates that the b axis is to lie in the plane formed ******/
/********* by the y-axis and the x-axis i.e. perp. to Z               ******/

          *(p_latt_vec+3)= mag_b * cos( gamma);
          *(p_latt_vec+4)= mag_b * sin( gamma);
          *(p_latt_vec+5)= 0.0;

/******** XYZ indicates that the c vector position is now defined **********/

          *(p_latt_vec+6)= mag_c * cos( beta); 

          bx_cx = ( *(p_latt_vec+3)) * ( *(p_latt_vec+6));

          *(p_latt_vec+7)= ( mag_b*mag_c*cos( alpha)  - bx_cx )/ ( *(p_latt_vec+4));

          *(p_latt_vec+8)= sqrt ( mag_c*mag_c - *(p_latt_vec+6) * *(p_latt_vec+6) 
                                              - *(p_latt_vec+7) * *(p_latt_vec+7)); 

      }

  if (use_zyx)
    {
/********* Z indicates that c is to be aligned with the z axis *********/
          *(p_latt_vec+6)= 0.0;
          *(p_latt_vec+7)= 0.0;
          *(p_latt_vec+8)= mag_c;

/********* ZY indicates that the b axis is to lie in the plane formed ******/
/********* by the y-axis and the z-axis i.e. perp. to X               ******/

          *(p_latt_vec+3)= 0.0;
          *(p_latt_vec+4)= mag_b * sin( alpha);
          *(p_latt_vec+5)= mag_b * cos( alpha);

/******** ZYX indicates that the a vector position is now defined **********/

          *(p_latt_vec+2)= mag_a * cos( beta); 

          az_bz = ( *(p_latt_vec+2)) * ( *(p_latt_vec+5));

          *(p_latt_vec+1)= ( mag_a*mag_b*cos( gamma)  - az_bz )/ ( *(p_latt_vec+4));

          *p_latt_vec    = sqrt ( mag_a*mag_a - *(p_latt_vec+1) * *(p_latt_vec+1) 
                                              - *(p_latt_vec+2) * *(p_latt_vec+2)); 

    }
/****** work out the magnitude of each lattice vector **********************/

    *p_real_latt_sizes = *p_latt_vec * *p_latt_vec +
                         *(p_latt_vec+1) * *(p_latt_vec+1) +
                         *(p_latt_vec+2) * *(p_latt_vec+2);

    *p_real_latt_sizes = sqrt(*p_real_latt_sizes);

    *(p_real_latt_sizes+1) = *(p_latt_vec+3) * *(p_latt_vec+3) +
                             *(p_latt_vec+4) * *(p_latt_vec+4) +
                             *(p_latt_vec+5) * *(p_latt_vec+5);

    *(p_real_latt_sizes+1) = sqrt(*(p_real_latt_sizes+1));
   
    *(p_real_latt_sizes+2) = *(p_latt_vec+6) * *(p_latt_vec+6) +
                             *(p_latt_vec+7) * *(p_latt_vec+7) +
                             *(p_latt_vec+8) * *(p_latt_vec+8);

    *(p_real_latt_sizes+2) = sqrt(*(p_real_latt_sizes+2));

/****** re-calculate alpha, beta, gamma as a check **************************/

    alpha_check = RAD_TO_DEG* acos( vec_dot( (p_latt_vec+3), (p_latt_vec+6))/
                                    (*(p_real_latt_sizes+1) * *(p_real_latt_sizes+2)));

     beta_check = RAD_TO_DEG* acos( vec_dot( p_latt_vec,     (p_latt_vec+6))/
                                    (*p_real_latt_sizes * *(p_real_latt_sizes+2)));

    gamma_check = RAD_TO_DEG* acos( vec_dot( (p_latt_vec+3), p_latt_vec)/
                                    (*(p_real_latt_sizes+1) * *p_real_latt_sizes));

    printf("\n        car file supplied a= %10.6f, b=%10.6f, c=%10.6f, alpha=%10.6f, beta=%10.6f, gamma=%10.6f\n",
                               mag_a, mag_b, mag_c, RAD_TO_DEG*alpha, RAD_TO_DEG*beta, RAD_TO_DEG*gamma);

    printf("cartessian forms generate a= %10.6f, b=%10.6f, c=%10.6f, alpha=%10.6f, beta=%10.6f, gamma=%10.6f\n",
                               *p_real_latt_sizes,
                               *(p_real_latt_sizes+1),
                               *(p_real_latt_sizes+2),
                               alpha_check, beta_check, gamma_check);

    printf("\nCartessian vectors for real space lattice:\n");
    printf("%10.6f %10.6f %10.6f\n", *p_latt_vec,     *(p_latt_vec+1), *(p_latt_vec+2));
    printf("%10.6f %10.6f %10.6f\n", *(p_latt_vec+3), *(p_latt_vec+4), *(p_latt_vec+5));
    printf("%10.6f %10.6f %10.6f\n\n", *(p_latt_vec+6), *(p_latt_vec+7), *(p_latt_vec+8));
   
    
/****** generate reciprocal lattice vectors ********************************/

    vec_cross( p_latt_vec,   p_latt_vec+3, &a_cross_b[0]); 
    vec_cross( p_latt_vec+3, p_latt_vec+6, &b_cross_c[0]); 
    vec_cross( p_latt_vec+6, p_latt_vec  , &c_cross_a[0]); 

    *p_cell_volume= vec_dot(p_latt_vec, &b_cross_c[0]); 

    printf(" a cross b = %10.6f %10.6f %10.6f \n", a_cross_b[0],a_cross_b[1],a_cross_b[2]);
    printf(" b cross c = %10.6f %10.6f %10.6f \n", b_cross_c[0],b_cross_c[1],b_cross_c[2]);
    printf(" c cross a = %10.6f %10.6f %10.6f \n", c_cross_a[0],c_cross_a[1],c_cross_a[2]);

    *p_recip_latt_vec    = b_cross_c[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+1)= b_cross_c[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+2)= b_cross_c[2]/ (*p_cell_volume);

    *(p_recip_latt_vec+3)= c_cross_a[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+4)= c_cross_a[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+5)= c_cross_a[2]/ (*p_cell_volume);

    *(p_recip_latt_vec+6)= a_cross_b[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+7)= a_cross_b[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+8)= a_cross_b[2]/ (*p_cell_volume);

/**************** Work out recip vector magnitudes ********************/

    *p_recip_latt_sizes =  *p_recip_latt_vec     * *p_recip_latt_vec
                         + *(p_recip_latt_vec+1) * *(p_recip_latt_vec+1)
                         + *(p_recip_latt_vec+2) * *(p_recip_latt_vec+2);

    *p_recip_latt_sizes = sqrt(*p_recip_latt_sizes);

    *(p_recip_latt_sizes+1) =  *(p_recip_latt_vec+3) * *(p_recip_latt_vec+3)
                             + *(p_recip_latt_vec+4) * *(p_recip_latt_vec+4)
                             + *(p_recip_latt_vec+5) * *(p_recip_latt_vec+5);
 
    *(p_recip_latt_sizes+1) = sqrt(*(p_recip_latt_sizes+1));

    *(p_recip_latt_sizes+2) =  *(p_recip_latt_vec+6) * *(p_recip_latt_vec+6)
                             + *(p_recip_latt_vec+7) * *(p_recip_latt_vec+7)
                             + *(p_recip_latt_vec+8) * *(p_recip_latt_vec+8);

    *(p_recip_latt_sizes+2) = sqrt(*(p_recip_latt_sizes+2));

/**************** Set up Ewald sum parameters *************************/
    
    printf("one_third = %10.6f, cell_volume= %10.6f\n\n",one_third, *p_cell_volume);
    printf("one_sixth = %10.6f\n\n",one_sixth);
    cell_length= pow(*p_cell_volume, one_third);

    printf("pi_tothehalf = %10.6f, num_atoms= %d, cell_length= %10.6f\n\n",
                                   pi_tothehalf,num_atoms,cell_length);

    kappa = pi_tothehalf * pow((num_atoms+1), one_sixth)/ cell_length;

/****************************************************************/
/***DEBUG This is the Ewald parameter if all is well altering ***/
/***DEBUG kappa should just divi up the work between the real ***/
/***DEBUG and reciprocal parts of the sum differently. Use to ***/
/***DEBUG test that the sum is working!!! Dave Willock June 97***/
/*      kappa = 0.5*kappa;                                      */
/***DEBUG END                                                 ***/
/****************************************************************/

    four_pi_over_vol= four_pi/(*p_cell_volume);
    printf("four_pi_over_vol= %10.6f\n",four_pi_over_vol);

    pi_sqrd_over_kappa_sqrd= pi*pi/(kappa*kappa);
    printf("pi_sqrd_over_kappa_sqrd= %10.6f\n",pi_sqrd_over_kappa_sqrd);

/***********************************************************************/
/***** factor of pi takes into account that my reciprocal space vecs ***/
/***** are simply a X b / cell_vol etc and Norgett uses 2pi times this */
/***** According to Catlow and Norgett best sum limits are           ***/
/*****           Gm = 2*kappa*f_param and Rm= f_param/kappa          ***/
/***********************************************************************/
   
    printf ("\nfparam = %10.6f, kappa = %10.6f\n\n",f_param, kappa);

    recip_sum_max = f_param * kappa / pi;
    real_sum_max  = f_param / kappa;

    recip_sum_max2 = recip_sum_max*recip_sum_max;
    real_sum_max2  = real_sum_max*real_sum_max;
    
/**************** Check by forming all dot products this should *******/
/**************** give the identity matrix ****************************/

    dot_set[0] = vec_dot(p_recip_latt_vec, p_latt_vec);    
    dot_set[1] = vec_dot(p_recip_latt_vec+3, p_latt_vec);    
    dot_set[2] = vec_dot(p_recip_latt_vec+6, p_latt_vec);    

    dot_set[3] = vec_dot(p_recip_latt_vec, p_latt_vec+3);    
    dot_set[4] = vec_dot(p_recip_latt_vec+3, p_latt_vec+3);    
    dot_set[5] = vec_dot(p_recip_latt_vec+6, p_latt_vec+3);    

    dot_set[6] = vec_dot(p_recip_latt_vec, p_latt_vec+6);    
    dot_set[7] = vec_dot(p_recip_latt_vec+3, p_latt_vec+6);    
    dot_set[8] = vec_dot(p_recip_latt_vec+6, p_latt_vec+6);    


       printf ("Check on recip/real lat vecs\n\n");

       for (iloop = 0; iloop < 9; iloop++)
         {
            printf("%10.6f  ", dot_set[iloop]);
            if (iloop == 2 || iloop == 5 || iloop == 8) printf("\n");
         }

    return;
  }

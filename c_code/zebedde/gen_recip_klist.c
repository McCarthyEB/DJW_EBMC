/****************************************************************************************/
/***** Generator for the k vector list for recip Ewald sum ******************************/
/***** Started 18 March 1998 Dave Willock  **********************************************/
/****************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "own_maths.h"
#include "ewald.h"

void  gen_recip_klist( double *p_kvecs, double *p_kvec2, double *p_gvec2, int *p_num_kvecs,
                       atom *p_host, int num_host_atoms, double *p_cos_sum, double *p_sin_sum)
{
#include "header.h"

int iatom, index_a, index_b, index_c;
int num_needed[3];

double k_sqrd, kx, ky, kz;
double kxa, kya, kza, kxb, kyb, kzb;
double g_dot_r_cell;

atom *p_this_host;

num_needed[0]= 1+recip_sum_max/recip_latt_sizes[0]; 
num_needed[1]= 1+recip_sum_max/recip_latt_sizes[1]; 
num_needed[2]= 1+recip_sum_max/recip_latt_sizes[2]; 

/*********************************************************************/
/********* Loop over reciprocal space vectors ************************/
/*********************************************************************/

for ( index_a= -num_needed[0]; index_a <= num_needed[0]; index_a++)
  {
    kxa= index_a*recip_latt_vec[0];
    kya= index_a*recip_latt_vec[1];
    kza= index_a*recip_latt_vec[2];

    for ( index_b= -num_needed[1]; index_b <= num_needed[1]; index_b++)
      {
        kxb = kxa+ index_b*recip_latt_vec[3];
        kyb = kya+ index_b*recip_latt_vec[4];
        kzb = kza+ index_b*recip_latt_vec[5];

        for ( index_c= -num_needed[2]; index_c <= num_needed[2]; index_c++)
          {
             if (*p_num_kvecs < MAX_KVECS)
               {
                 kx = kxb+ index_c*recip_latt_vec[6];
                 ky = kyb+ index_c*recip_latt_vec[7];
                 kz = kzb+ index_c*recip_latt_vec[8];
            
                 k_sqrd= kx*kx + ky*ky + kz*kz;

/********* Avoid k=0 term ********************************************/
            
                 if (k_sqrd > 1E-6 && k_sqrd <= recip_sum_max2 )
                   {
                      *p_kvecs = kx;
                      p_kvecs++;
                      *p_kvecs = ky;
                      p_kvecs++;
                      *p_kvecs = kz;
                      p_kvecs++;
           
                      *p_kvec2 = kx*kx + ky*ky + kz*kz;
                      p_kvec2++;
  
                      *p_gvec2 = four_pi_sqrd*k_sqrd;
                      p_gvec2++;

/*******************Since static host can get sin_sum and cos_sum once off!***/
                      *p_sin_sum=0;
                      *p_cos_sum=0;

                      p_this_host= p_host-1;
                      for (iatom=0; iatom <= num_host_atoms; iatom++)
                         {
                            p_this_host++;

                            g_dot_r_cell= two_pi*( kx* (p_this_host->x)
                                                 + ky* (p_this_host->y)
                                                 + kz* (p_this_host->z));

                            (*p_sin_sum) += (p_this_host->part_chge)* sin(g_dot_r_cell);
                            (*p_cos_sum) += (p_this_host->part_chge)* cos(g_dot_r_cell);
                         }

                      p_sin_sum++;
                      p_cos_sum++;
                      (*p_num_kvecs)++; 
                    }
               }
          }
      }
  }

if ( *p_num_kvecs > MAX_KVECS ) 
  {
    printf("ERROR : The number of k-vectors required for this problem exceeds the");
    printf(" maximum allowed\n");
    printf("%d have been requested when the maximum is %d\n", *p_num_kvecs, MAX_KVECS);
    printf("Lower the Ewald accuracy, or re-parameterise the program\n");
    exit(0);
  }

return;
}

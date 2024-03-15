/**************************************************************************/
/******* H-bond energy calculation assuming AMBER force field *************/
/******* Dave Willock started 4th March 1998                  *************/
/**************************************************************************/


#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

void min_image( double *x, double *y, double *z);

int pbc_interactions(double *p_dx, double *p_dy, double *p_dz, double cutoff,
                     double cutoff_2, double *p_pos_separations,  double *p_pos_vectors, 
                     int count_flag, int is_self);

double hbond_amber(atom *p_pore, int num_p_atoms,
                   atom *p_template, int num_t_atoms, int *p_need_grad,
                   double *p_grad, int is_self)
{
#include "header.h"
atom *p_tatom, *p_patom;

int i,j,pot_index_t, pot_index_p,this_pair;
int index_temp_atom,num_interactions; 
int index;

double dx,dy,dz,r,r2,r3,r5,r10,r12;
double hb_energy, hb_energy_rep, hb_energy_disp;
double pos_separations[MAX_PAIR_LIST], pos_vectors[3*MAX_PAIR_LIST];
double rep,disp,grad_factor;
double aaa,bbb;

double *p_grad_x, *p_grad_y, *p_grad_z;

/**************************************************************************/
/********** loop over template atoms :                 ********************/
/********** p_tatom = pointer to current template atom ********************/
/**************************************************************************/

hb_energy_rep = 0.0;
hb_energy_disp = 0.0;

p_tatom = p_template-1;     
for (index_temp_atom=0; index_temp_atom<=num_t_atoms ;index_temp_atom++)
   {
     p_tatom++;
     pot_index_t = p_tatom->hb_list;

     if (pot_index_t >= 0)
       {

          p_grad_x = p_grad+ 3*index_temp_atom;
          p_grad_y = p_grad_x+1;
          p_grad_z = p_grad_y+1;

/**************************************************************************/
/********** loop over pore atoms :                     ********************/
/********** p_patom = pointer to current pore atom     ********************/
/**************************************************************************/

          p_patom = p_pore-1;     

          for ( j=0; j<=num_p_atoms; j++)
             {
               p_patom++;
               pot_index_p = p_patom->hb_list;
             
               if (pot_index_p >= 0)
                  {

                     dx = p_tatom->x - p_patom->x;
                     dy = p_tatom->y - p_patom->y;
                     dz = p_tatom->z - p_patom->z;

/**************************************************************************/
/******* if this is a periodic pore assemble a list of interactions *******/
/******* of this type otherwise check if the pair is in range       *******/
/**************************************************************************/

                     if (pbc)
                       {
                          num_interactions = pbc_interactions( &dx, &dy, &dz, hb_ctf, hb_ctf_2,  
                                                               &pos_separations[0], &pos_vectors[0], FALSE, 
                                                               is_self);
                       }
                     else
                       {
                          r2 = dx*dx + dy*dy + dz*dz;
            
                          if (r2 <= hb_ctf_2)
                            {
                               num_interactions=0;
                               pos_separations[0]= sqrt(r2);
                               pos_vectors[0]= dx;
                               pos_vectors[1]= dy;
                               pos_vectors[2]= dz;
                            }
                          else
                            {
                               num_interactions=-1;
                            }
                       }

/*************************************************************************/
/******* Only bother with list if there are some members in it! **********/
/*************************************************************************/

                     if (num_interactions >= 0)
                       {

/*************************************************************************/
/******* Assume hbond parameters are in a sorted list          ***********/
/******* according to order_double_pot.c                       ***********/
/*************************************************************************/
/******* look up index according to atom potentials **********************/


                          if ( pot_index_t < pot_index_p )
                            {
                              index= pot_index_t*(num_h_pot_types+1);
                              for (i=1; i < pot_index_t; i++) index-=i;
                              index += pot_index_p-pot_index_t;
                            }
                          else
                            {
                              index= pot_index_p*(num_h_pot_types+1);
                              for (i=1; i < pot_index_p; i++) index-=i;
                              index += pot_index_t-pot_index_p;
                            }
          
                          aaa= h_potent[index].a;
                          bbb= h_potent[index].b; 

/*************************************************************************/

                          rep = 0;
                          disp= 0;
                          for (this_pair= 0; this_pair <= num_interactions; this_pair++) 
                            {
                               r = pos_separations[this_pair];

/******* Avoid very close contacts and report ****************************/

                               if (r <= 0.1) 
                                 {
                                   printf("Warning two atoms are at a separation less than 0.1A,\n");
                                   printf("Have discounted their contribution to the energy.\n");
                                   hb_energy_rep       += 10000;
                       
                                 }

/**** Ignore interactions with no parameters set i.e. none hydrogen bonds ****/

                               else if (aaa > 0.00001)
                                 {
                                   r2 = r*r;
                                   r3 = r2*r;
                                   r5 = r2 * r3;
                                   r10= r5 * r5;
                                   r12= r10* r2;

                
/******* Repulsion part of potential *************************************/

                                   rep = aaa/r12;

/******* Dispersion part of potential ************************************/

                                   disp= bbb/r10;

/****** Will calculate the derivatives ***********************************/

                                   if (*p_need_grad)
                                    {
                                       grad_factor= (10.0*disp-12.0*rep)/r2; 
                                       *p_grad_x += grad_factor* pos_vectors[3*this_pair]  ;
                                       *p_grad_y += grad_factor* pos_vectors[3*this_pair+1];
                                       *p_grad_z += grad_factor* pos_vectors[3*this_pair+2];
                                     }

                                   p_tatom->vdw_energy += rep-disp;
                                   p_patom->vdw_energy += rep-disp;
                                   hb_energy_rep       += rep; 
                                   hb_energy_disp      -= disp;
                                }
                            }
                        }
                    }
                 }
              }
           }

hb_energy= hb_energy_rep + hb_energy_disp;
return hb_energy;
}

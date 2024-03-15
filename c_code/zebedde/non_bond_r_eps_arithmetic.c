/**************************************************************************/
/******* Non-bond energy calculation assuming r-eps potential *************/
/******* format and arithmetic combination rules              *************/
/******* Dave Willock Added for AMBER forcefield Jan. 99.     *************/
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

double non_bond_r_eps_arithmetic(atom *p_pore, int num_p_atoms,
                                 atom *p_template, int num_t_atoms, int *p_need_grad,
                                 double *p_grad, double *p_pos_separations, 
                                 double *p_pos_vectors, int is_self)
{
#include "header.h"
atom *p_tatom, *p_patom;

int j,pot_index_t, pot_index_p,this_pair;
int index_temp_atom,num_interactions; 

double dx,dy,dz,r0,r02,r04,r06,r012,r,r2,r4,r6,r12, eps;
double nb_energy, nb_energy_rep, nb_energy_disp;
/* Replace parameter for array size with a malloc in calling routine (energy.c) JL Oct. 2007*/
/*double pos_separations[MAX_PAIR_LIST], pos_vectors[3*MAX_PAIR_LIST];*/
double rep,disp,grad_factor;

double *p_grad_x, *p_grad_y, *p_grad_z;

/**************************************************************************/
/********** loop over template atoms :                 ********************/
/********** p_tatom = pointer to current template atom ********************/
/**************************************************************************/

nb_energy_rep = 0.0;
nb_energy_disp = 0.0;

p_tatom = p_template-1;     
for (index_temp_atom=0; index_temp_atom<=num_t_atoms ;index_temp_atom++)
   {
     p_tatom++;
     pot_index_t = p_tatom->nb_list;

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
          pot_index_p = p_patom->nb_list;

          dx = p_tatom->x - p_patom->x;
          dy = p_tatom->y - p_patom->y;
          dz = p_tatom->z - p_patom->z;

/**************************************************************************/
/******* if this is a periodic pore assemble a list of interactions *******/
/******* of this type otherwise check if the pair is in range       *******/
/**************************************************************************/

          if (pbc)
            {
               num_interactions = pbc_interactions( &dx, &dy, &dz, nb_ctf, nb_ctf_2,  
                                                    &p_pos_separations[0], &p_pos_vectors[0], 
                                                    FALSE, is_self);
            }
          else
            {
               r2 = dx*dx + dy*dy + dz*dz;
            
               if (r2 <= nb_ctf_2)
                 {
                    num_interactions=0;
                    p_pos_separations[0]= sqrt(r2);
                    p_pos_vectors[0]= dx;
                    p_pos_vectors[1]= dy;
                    p_pos_vectors[2]= dz;
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
/******* Apply combining rules to the atom potent parameters   ***********/
/******* using the nb_list parameter each atom should have for ***********/
/******* referencing the potent list                           ***********/
/******* When reading r-eps potentials a= r; b= eps            ***********/
/*************************************************************************/

                eps = potent[pot_index_t].sqrt_b * 
                      potent[pot_index_p].sqrt_b; 

                r0 = 0.5*(potent[pot_index_t].a + potent[pot_index_p].a);
 
                r02 = r0 * r0;
                r04 = r02 * r02;
                r06 = r04 * r02;
                r012 = r06 * r06;

/*************************************************************************/

               rep = 0;
               disp= 0;
               for (this_pair= 0; this_pair <= num_interactions; this_pair++) 
                 {
                    r = p_pos_separations[this_pair];

/******* Avoid very close contacts and report ****************************/

                    if (r <= 0.1) 
                      {
                        if (!told_of_clash)
                          {
                            printf("WARNING: Two atoms are at a separation less than 0.1A,\n");
                            printf("Have capped their contribution to the energy at 10000.\n");
                            printf("Multiple warnings suppressed until next build success or accepted move if MC\n");
                            told_of_clash=TRUE;
                          }

                        p_tatom->vdw_energy += 10000;
                        nb_energy_rep       += 10000;
                       
                      }
                    else
                      {
                        r2 = r * r;
                        r4 = r2 * r2;
                        r6 = r4 * r2;
                        r12 = r6 * r6;

/******* Repulsion part of potential *************************************/

                        rep = eps * r012/r12;

/******* Dispersion part of potential ************************************/

                        disp= 2.0 * eps * r06/r6;

/****** Will calculate the derivatives ***********************************/

                        if (*p_need_grad)
                         {
                            grad_factor= 12.0*rep/r2;
                            *p_grad_x -= grad_factor* p_pos_vectors[3*this_pair]  ;
                            *p_grad_y -= grad_factor* p_pos_vectors[3*this_pair+1];
                            *p_grad_z -= grad_factor* p_pos_vectors[3*this_pair+2];

                            grad_factor= 6.0 * disp / r2 ;
                            *p_grad_x += grad_factor* p_pos_vectors[3*this_pair]  ;
                            *p_grad_y += grad_factor* p_pos_vectors[3*this_pair+1];
                            *p_grad_z += grad_factor* p_pos_vectors[3*this_pair+2];
                          }

                        p_tatom->vdw_energy += rep-disp;
                        p_patom->vdw_energy += rep-disp;
                        nb_energy_rep       += rep; 
                        nb_energy_disp      -= disp;
                     }
                 }
             }
         }
   }

nb_energy= nb_energy_rep + nb_energy_disp;
return nb_energy;
}

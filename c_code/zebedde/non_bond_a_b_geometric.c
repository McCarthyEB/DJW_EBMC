/**************************************************************************/
/******* Non-bond energy calculation assuming A-B(12-6) potential *********/
/******* format and geometric combination rules                   *********/
/******* Dave Willock March 1997                                  *********/
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

double non_bond_a_b_geometric(atom *p_pore, int num_p_atoms,
                              atom *p_template, int num_t_atoms, int *p_need_grad,
                              double *p_grad, double *p_pos_separations, double *p_pos_vectors, 
                              int is_self)
{
#include "header.h"
atom *p_tatom, *p_patom;

int j,pot_index_t, pot_index_p,this_pair;
int index_temp_atom,num_interactions; 

double dx,dy,dz,r,r2,r3,r6,r8,r12;
double nb_energy, nb_energy_rep, nb_energy_disp;
/* Replace parameter for array size with a malloc in calling routine (energy.c) JL Oct. 2007*/
/*double pos_separations[MAX_PAIR_LIST], pos_vectors[3*MAX_PAIR_LIST];*/
double rep,disp,grad_factor;
double aij, bij;

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
 
/******* Catch very close interactions *************************************/
               if (r2 < 0.01) return 10000;
            
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
/******* Geometric combining rules for cvff forcefields:       ***********/
/*******                                                       ***********/
/******* > E = Aij/r^12 - Bij/r^6                              ***********/
/******* > where  Aij = sqrt( Ai * Aj )                        ***********/
/******* >        Bij = sqrt( Bi * Bj )                        ***********/
/*************************************************************************/

               aij= potent[pot_index_t].sqrt_a*potent[pot_index_p].sqrt_a;
               bij= potent[pot_index_t].sqrt_b*potent[pot_index_p].sqrt_b;

/*************************************************************************/

               rep = 0;
               disp= 0;
               for (this_pair= 0; this_pair <= num_interactions; this_pair++) 
                 {
                    r = p_pos_separations[this_pair];

                    r2 = r*r;
                    r3 = r2*r;
                    r6 = r3 * r3;
                    r12= r6 * r6;
                
/******* Repulsion part of potential *************************************/

                    rep = aij/r12;

/******* Dispersion part of potential ************************************/

                    disp= bij/r6;

/****** Will calculate the derivatives ***********************************/

                    if (*p_need_grad)
                     {
                        grad_factor= 12.0*aij/r2;
                        *p_grad_x -= grad_factor* p_pos_vectors[3*this_pair]  ;
                        *p_grad_y -= grad_factor* p_pos_vectors[3*this_pair+1];
                        *p_grad_z -= grad_factor* p_pos_vectors[3*this_pair+2];

                        r8 = r6 * r2;
                        grad_factor= 6.0*bij/r2;
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

nb_energy= nb_energy_rep + nb_energy_disp;
return nb_energy;
}

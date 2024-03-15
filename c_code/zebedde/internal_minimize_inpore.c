/*********************************************************************/
/* internal_minimize_inpore.c : minimizes a molecule with respect    */	
/*                       to (fixed) pore using internal minimiser    */
/* started 19/6/ DWL updated  27/11/96 DJW                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

#define INTERNAL_SUCCESS 1

void print_product(double *p_matrix, double *p_vector, int dimension);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, 
                    atom *p_molecule, int num_atoms, int which_mol );

void apply_delta_rigid(atom *p_molecule, int num_atoms, 
                       double *p_c_of_m, double *p_delta);

void rigid_body_gradients(atom *p_molecule, int num_atoms,
                          double *p_c_of_m,
                          double *p_atom_grads, double *p_central_grad);

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *p_templ, int num_t_atoms,
                      double *p_kvecs, double *p_kvec2,
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad,
                      int have_comb_rules, int is_empty_pore);

void bfgs_update(double *p_G, double *p_gamma, double *p_delta, int dimension);

double bra_x_ket( double *p_a, double *p_b, int dimension );

void matrix_x_vector( double *p_matrix, double *p_ai, double *p_bi,
                      int dimension );

void internal_minimize_inpore(atom *p_pore, int num_p_atoms,
                              atom *p_molecule, int num_t_atoms, 
                              double *p_kvecs, double *p_kvec2,
                              double *p_gvec2, int num_kvecs,
                              double *p_cos_sum, double *p_sin_sum,
                              int *p_need_grad, double *p_grad,
                              int have_comb_rules, int is_empty_pore)

{
#include "header.h"
int index, index_grad, index_centre;
int iloop, inc, index_diag, need_grad;

atom moved_molecule[MAXTEMPLATE];

double  moved_c_of_m[3], gamma[6], delta[6];
double moved_grad[MAXTEMPLATE3], moved_central_grad[6];
double *p_this_grad, c_of_m[3], total_mass;
double central_grad[6], central_hessian[21];
double step_scale, grad_dot_delta;

printf("Entering internal minimiser!!!!\n");
printf("have %d pore atoms, %d template atoms\n", num_p_atoms, num_t_atoms);

printf("Need centre of mass for template to do torques\n"); 

centre_of_mass(&c_of_m[0], &total_mass, p_molecule, num_t_atoms, -1 );

printf("Centre of mass is at: %10.6f %10.6f %10.6f \n",
                         c_of_m[0], c_of_m[1], c_of_m[2] );
printf("current grad vector:\n");

p_this_grad= p_grad;
for (index_grad=0; index_grad <= 3*num_t_atoms+2; index_grad++)
  {
    printf ("%d => %10.6f\n",index_grad, *p_this_grad);
    p_this_grad++;
  }

/**** For rigid body minimisation work out central force and central torque ****/

rigid_body_gradients(p_molecule, num_t_atoms, &c_of_m[0], p_grad, &central_grad[0]); 

/****** Set unit starting hessian *******************/

for (index_centre=0; index_centre <21; index_centre++) central_hessian[index_centre]=0;

index_diag = 0;
inc=6;
for(index=0; index < 6; index++)
  {
    central_hessian[index_diag] = 1;
    index_diag += inc;
    inc--;
  }

printf("Centralised Gradients:\n");

print_product(&central_hessian[0], &central_grad[0], 6);

/***************************************************************/
/**** try out bfgs update routine                          *****/
/**** Start off with steepest descent to get second point! *****/
/***************************************************************/

step_scale= 0.1;
for (index=0; index < 6; index++)
   {
     delta[index]= -step_scale* central_grad[index];
   }

for (index=0; index <= num_t_atoms; index++) moved_molecule[index]= *(p_molecule+index);
for (index=0; index < 3;  index++) moved_c_of_m[index]= c_of_m[index];

printf("Initially moved_c_of_m = %10.6f %10.6f %10.6f \n", 
                              moved_c_of_m[0], moved_c_of_m[1], moved_c_of_m[2]);

printf("Taking a step into the dark:\n");
for (index=0; index < 6; index++)
   {
     printf("%d => %10.6f\n", index,  delta[index]);
   }

grad_dot_delta= bra_x_ket(&central_grad[0], &delta[0], 6);
printf("initial G.D product= %10.6f\n", grad_dot_delta);

apply_delta_rigid(&moved_molecule[0], num_t_atoms, &moved_c_of_m[0], &delta[0]);

printf("New Centre of mass: %10.6f %10.6f %10.6f \n", 
                              moved_c_of_m[0],moved_c_of_m[1],moved_c_of_m[2]);



/****************************************************************/
/**** Get new Gradients!!!  *************************************/
/****************************************************************/

calculate_energy(p_pore, num_p_atoms, &moved_molecule[0], num_t_atoms,
                 p_kvecs, p_kvec2, p_gvec2, num_kvecs, p_cos_sum, p_sin_sum,
                 &need_grad, &moved_grad[0], have_comb_rules, is_empty_pore);

printf("Grad at moved position:\n");
for (index_grad=0; index_grad <= 3*num_t_atoms+2; index_grad++)
  {
    printf ("%d => %10.6f\n",index_grad, moved_grad[index_grad]);
  }

rigid_body_gradients(&moved_molecule[0], num_t_atoms, &moved_c_of_m[0], 
                     &moved_grad[0], &moved_central_grad[0]); 

grad_dot_delta= bra_x_ket(&moved_central_grad[0], &delta[0], 6);
printf("after first steepest descent G.D= %10.6f\n", grad_dot_delta);

printf("New gradient:\n");
for (index=0; index < 6; index++)
   {
     printf("%d => %10.6f\n", index, moved_central_grad[index]);     
   }

for (index=0; index < 6; index++) 
                     gamma[index]= moved_central_grad[index] - central_grad[index];

printf("gamma:\n");
for (index=0; index < 6; index++)
   {
     printf("%d => %10.6f\n", index, gamma[index]);     
   }

/***************************************************************/
/*** do a few loops round using bfgs to minimise ***************/
/***************************************************************/

for (iloop=0; iloop < 10; iloop++)
  {

    bfgs_update(&central_hessian[0], &gamma[0], &delta[0], 6);

printf("After bfgs_update:\n");

print_product(&central_hessian[0], &central_grad[0], 6);

matrix_x_vector( &central_hessian[0], &central_grad[0], &delta[0], 6);

/**************************************************************/
/********* Apply the move *************************************/
/**************************************************************/

for (index=0; index < 6; index++) delta[index]= step_scale*delta[index];
apply_delta_rigid(&moved_molecule[0], num_t_atoms, &moved_c_of_m[0], &delta[0]);

printf("New Centre of mass: %10.6f %10.6f %10.6f \n",
                              moved_c_of_m[0],moved_c_of_m[1],moved_c_of_m[2]);

/****************************************************************/
/**** Get new Gradients!!!  *************************************/
/****************************************************************/

need_grad=TRUE;
calculate_energy(p_pore, num_p_atoms, &moved_molecule[0], num_t_atoms,
                 p_kvecs, p_kvec2, p_gvec2, num_kvecs, p_cos_sum, p_sin_sum,
                 &need_grad, &moved_grad[0], have_comb_rules, is_empty_pore);

rigid_body_gradients(&moved_molecule[0], num_t_atoms, &moved_c_of_m[0],
                     &moved_grad[0], &moved_central_grad[0]);

for (index=0; index < 6; index++)
                     gamma[index]= moved_central_grad[index] - central_grad[index];

grad_dot_delta= bra_x_ket(&moved_central_grad[0], &delta[0], 6);
printf("after %dth step descent G.D= %10.6f\n", iloop+2, grad_dot_delta);
  }

return;
}

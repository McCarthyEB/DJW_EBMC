/**************************************************************************/
/******* Routine to calculate the maximum number of vdw interactions ******/
/******* for this particular instance to allow MALLOC of      *************/
/******* pbc_interactions                                     *************/
/*******           James Landon and  Dave Willock Oct.2007    *************/
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

int count_vdw_pbc(atom *p_pore, int num_p_atoms,
                  atom *p_template, int num_t_atoms, int is_self)
{
#include "header.h"
atom *p_tatom, *p_patom;

int j;
int index_temp_atom,num_interactions, max_interactions; 
int pot_index_1, pot_index_2;
int ilist, told;

double dx,dy,dz, r2, r2min;
double dxmin,dymin,dzmin;
double pos_separation, pos_vector[3];

double *p_grad_x, *p_grad_y, *p_grad_z;

/**************************************************************************/
/********** loop over template atoms :                 ********************/
/********** p_tatom = pointer to current template atom ********************/
/**************************************************************************/

//printf("Counting interactions for include for num_p_atoms=%d and num_t_atoms=%d\n",
//                              num_p_atoms, num_t_atoms);

if (!pbc) return 1;

r2min = nb_ctf_2;
dxmin = nb_ctf;
dymin = nb_ctf;
dzmin = nb_ctf;

p_tatom = p_template-1;     
for (index_temp_atom=0; index_temp_atom<=num_t_atoms ;index_temp_atom++)
   {
     p_tatom++;
     /*pot_index_t = p_tatom->nb_list;*/

/**************************************************************************/
/********** loop over pore atoms :                     ********************/
/********** p_patom = pointer to current pore atom     ********************/
/**************************************************************************/

     if (!is_self)
       {
         for ( j=0; j<=num_p_atoms; j++)
            {
              p_patom=p_pore+j;

              dx = p_tatom->x - p_patom->x;
              dy = p_tatom->y - p_patom->y;
              dz = p_tatom->z - p_patom->z;

//              printf("Before min_image: %10.6f  %10.6f  %10.6f \n", dx, dy, dz);

              min_image(&dx, &dy, &dz);
 
              r2 = dx*dx + dy*dy + dz*dz;
//              printf("After min_image : %10.6f  %10.6f  %10.6f dist: %10.6f\n\n", dx, dy, dz, sqrt(r2));

              if (r2 < r2min)
                {
                  dxmin=dx;
                  dymin=dy;
                  dzmin=dz;
                  r2min = r2;
                }
             }
       }
    else
       {
//         printf("Molecule with own image....\n");
         for ( j=index_temp_atom+1; j<=num_p_atoms; j++)
            {
              p_patom=p_pore+j;

              dx = p_tatom->x - p_patom->x;
              dy = p_tatom->y - p_patom->y;
              dz = p_tatom->z - p_patom->z;

//              printf("Before min_image: %10.6f  %10.6f  %10.6f \n", dx, dy, dz);

              min_image(&dx, &dy, &dz);
 
              r2 = dx*dx + dy*dy + dz*dz;
//              printf("After min_image : %10.6f  %10.6f  %10.6f dist: %10.6f\n\n", dx, dy, dz, sqrt(r2));

              if (r2 < r2min)
                {
                  dxmin=dx;
                  dymin=dy;
                  dzmin=dz;
                  r2min = r2;
                }
             }
       }
   }

//printf("Minimum separation found = %10.6f A\n\n\n", sqrt(r2min));
/*************************************************************************/
/******* Find the number of interactions for the shortest vector *********/
/*************************************************************************/

 num_interactions = pbc_interactions( &dxmin, &dymin, &dzmin, nb_ctf, nb_ctf_2,  
                                                       &pos_separation, pos_vector, TRUE, is_self);

return num_interactions;
}

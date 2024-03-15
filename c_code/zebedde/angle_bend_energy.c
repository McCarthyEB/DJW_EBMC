/**************************************************************************/
/******* Intra-molecular angle bend energy calculation based on   *********/
/******* pre-assembled list.                                      *********/
/******* Dave Willock Aug 2006                                    *********/
/*******                                                          *********/
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

void min_image( double *x, double *y, double *z);

double angle_bend_energy(atom *p_molecule, int num_atoms,
                         angle_interact_list *p_angles_list, int num_angles_listed)
{
#include "header.h"
atom *p_atom1, *p_atom2, *p_atom3;

double vec1[3], vec2[3], theta, size1, size2;
double dot, dtheta, dtheta2, dtheta3, dtheta4;

int iangle, pot_index;

double energy;

energy = 0.0;

for (iangle=0; iangle<=num_angles_listed ;iangle++)
   {
      pot_index = p_angles_list->iangpot;

/******************************************************/
/** set atom pointers so angle 1-2-3 is what we want **/
/******************************************************/

      p_atom1 = p_molecule +  p_angles_list->iatm1;
      p_atom2 = p_molecule +  p_angles_list->iatm2;
      p_atom3 = p_molecule +  p_angles_list->iatm3;

      p_angles_list++;

/******************************************************/
/*** Work out angle ***********************************/
/******************************************************/

      vec1[0] = p_atom1->x - p_atom2->x;
      vec1[1] = p_atom1->y - p_atom2->y;
      vec1[2] = p_atom1->z - p_atom2->z;

      vec2[0] = p_atom3->x - p_atom2->x;
      vec2[1] = p_atom3->y - p_atom2->y;
      vec2[2] = p_atom3->z - p_atom2->z;

/******************************************************/
/** Need min-image vectors incase they break over *****/
/** pbc. Added Nov. 2006 Dave Willock             *****/
/******************************************************/
      min_image( &vec1[0], &vec1[1], &vec1[2]);
      min_image( &vec2[0], &vec2[1], &vec2[2]);

      dot = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];

      size1 = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2]);
      size2 = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2]);

      theta = acos(dot/(size1*size2));
//      if (DEBUG) printf("angle %s %s %s is %10.6f degrees\n",
//                                p_atom1->label, p_atom2->label,
//                                p_atom3->label, theta*RAD_TO_DEG);
      
/*******************************************************/
/*** Use appropriate functional form *******************/
/*******************************************************/

      if (intra_angle_potent[pot_index].which==QUADRATIC_ANGLE)
        {
           dtheta = (theta - intra_angle_potent[pot_index].A/RAD_TO_DEG);
           dtheta2= dtheta * dtheta;

           energy += intra_angle_potent[pot_index].B * dtheta2;

//           if (DEBUG) printf("Using quadratic energy term giving contribution %10.6f\n",
//                       intra_angle_potent[pot_index].B * dtheta2);
        }
      else if (intra_angle_potent[pot_index].which==QUARTIC_ANGLE)
        {
           dtheta = (theta - intra_angle_potent[pot_index].A/RAD_TO_DEG);
           dtheta2= dtheta * dtheta;
           dtheta3= dtheta2* dtheta;
           dtheta4= dtheta2* dtheta2;

//           if (DEBUG)
//             {
//                printf("dtheta: %10.6f rads sq: %10.6f cube: %10.6f pow4: %10.6f\n",
//                                                           dtheta, dtheta2, dtheta3, dtheta4);    
//
//                printf("coeffs: %10.6f %10.6f %10.6f \n",
//                             intra_angle_potent[pot_index].B,
//                             intra_angle_potent[pot_index].C,
//                             intra_angle_potent[pot_index].D);
//             }

           energy += intra_angle_potent[pot_index].B * dtheta2
                    +intra_angle_potent[pot_index].C * dtheta3
                    +intra_angle_potent[pot_index].D * dtheta4;

//           if (DEBUG) printf("Using quartic energy term giving contribution %10.6f\n",
//                                         intra_angle_potent[pot_index].B * dtheta2
//                                        +intra_angle_potent[pot_index].C * dtheta3
//                                        +intra_angle_potent[pot_index].D * dtheta4);
        }

//printf("DB>>ENERGY %s %s %s %lf\n",p_atom1->label,p_atom2->label,p_atom3->label, energy);

   }
  return energy;
}

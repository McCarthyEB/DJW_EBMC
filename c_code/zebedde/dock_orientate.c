/***************************************************************/
/******* Rigid body rotation to randomly orientate docked ******/
/******* molecule                                         ******/
/***************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "data.h"
#include "global_values.h"
#include "own_maths.h"
#include "constants.h"

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms, int which_mol);

double real_random(int done);

void unit_vector(double *p_vector, double *p_size);

void dock_orientate(atom *p_template, int num_template_atoms,
                    double *p_cofm) 
{
#include "header.h"
   int isign; 

   double axis[3], theta;
   double size;

/***************************************************************************/
/***** Start of executable lines *******************************************/
/***************************************************************************/

/***************************************************************************/
/******** set random angle and axis ****************************************/
/***************************************************************************/

         theta = real_random(1)*pi;
         isign =  real_random(1) < 0.5;
         if (isign) theta= -1.0*theta;

         printf("Dock orientate rotating by %10.6f\n", theta);

/********************************************************/
/*** This is limited to rots around X,Y or Z ************/
/*** could do a random vector in space by assigning *****/
/*** random components.                             *****/
/*** Must go to rotate with axis as a unit vector   *****/
/********************************************************/
         axis[0]= real_random(1);
         axis[1]= real_random(1);
         axis[2]= real_random(1);

         unit_vector(&axis[0], &size);

         printf("Selected axis %5.2f %5.2f %5.2f\n",
                                     axis[0],axis[1],axis[2]);

         rotate(p_template, &axis[0],
                          p_cofm, theta, num_template_atoms, -1);

return ;
}


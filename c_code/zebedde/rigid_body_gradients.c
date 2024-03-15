/****************************************************/
/*** Generate centre of mass gradients from atom ****/
/*** gradients. Dave Willock 27th Nov. 96        ****/
/****************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"

void rigid_body_gradients(atom *p_molecule, int num_atoms, 
                          double *p_c_of_m,            
                          double *p_atom_grads, double *p_central_grad)
{
int index_atoms, index_centre;

double *p_this_central_grad;
double *p_this_atom_grad;
double a[3];

atom *p_this_atom;

p_this_central_grad= p_central_grad;
for (index_centre=0; index_centre < 6; index_centre++) 
  {
    *p_this_central_grad=0;
    p_this_central_grad++;
  }

p_this_atom_grad= p_atom_grads;
p_this_atom= p_molecule;
for (index_atoms=0; index_atoms <= num_atoms; index_atoms++)
  {
    p_this_central_grad= p_central_grad;

/****** Central translational derivatives **************/
    *p_this_central_grad     += *p_this_atom_grad;
    *(p_this_central_grad+1) += *(p_this_atom_grad+1);
    *(p_this_central_grad+2) += *(p_this_atom_grad+2);


/****************************************************/
/****** Central rotational derivatives **************/
/****** a x grad !!!                   **************/
/****************************************************/
/****** get the vector a = centre of mass to atom ***/
/****************************************************/

    a[0] = p_this_atom->x - *p_c_of_m;
    a[1] = p_this_atom->y - *(p_c_of_m+1);
    a[2] = p_this_atom->z - *(p_c_of_m+2);

    printf("Rotat. div: %10.6f %10.6f %10.6f \n", a[0], a[1], a[2]);

/****************************************************/
/***** Do the cross product !! **********************/
/****************************************************/

    *(p_this_central_grad+3) += a[1]* (*(p_this_atom_grad+2)) - a[2]* (*(p_this_atom_grad+1));
    *(p_this_central_grad+4) += a[2]* (*p_this_atom_grad) - a[0]* (*(p_this_atom_grad+2));
    *(p_this_central_grad+5) += a[0]* (*(p_this_atom_grad+1)) - a[1]* (*p_this_atom_grad);

    p_this_atom_grad +=3;
    p_this_atom++;

}

return;
}

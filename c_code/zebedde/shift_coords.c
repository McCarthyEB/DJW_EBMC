/**************************************************************************/
/***** Subroutine to shift all the coordinates into the same box **********/
/***** Started AJWL 11/08                                        **********/
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void fract_to_cart(double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void cart_to_fract(double cart_x,  double cart_y,  double cart_z,
                    double *frac_a, double *frac_b, double *frac_c,
                    double *p_recip_latt_vec );

void shift_coords(atom *p_molecule, int num_atoms)
{
#include "header.h"
double fract_a,fract_b,fract_c;
int iatom;
int done;
atom *p_atom;
p_atom = p_molecule;


//printf("DB numatoms %i\n",num_atoms);
num_atoms++;
for (iatom=1;iatom<=num_atoms;iatom++)
{
//printf("DB iatom %i\n",iatom);
//printf("DB label %s \n",p_atom->label);
//printf("DB xyz>> %lf %lf %lf\n", p_atom->x,p_atom->y,p_atom->z);

cart_to_fract(p_atom->x,p_atom->y,p_atom->z, &fract_a, &fract_b, &fract_c, &recip_latt_vec[0]);

//printf("DB frac>> %lf %lf %lf\n", fract_a, fract_b, fract_c);

    done= FALSE;
    while (!done)
      {
        done = TRUE;
        if (fract_a >= 0.5) { fract_a= fract_a- 1.0; done = FALSE; }
        if (fract_a < -0.5) { fract_a= fract_a+ 1.0; done = FALSE; }

        if (fract_b >= 0.5) { fract_b= fract_b- 1.0; done = FALSE; }
        if (fract_b < -0.5) { fract_b= fract_b+ 1.0; done = FALSE; }

        if (fract_c >= 0.5) { fract_c= fract_c- 1.0; done = FALSE; }
        if (fract_c < -0.5) { fract_c= fract_c+ 1.0; done = FALSE; }
      }


fract_to_cart(&p_atom->x,&p_atom->y,&p_atom->z, fract_a, fract_b, fract_c, &latt_vec[0] );
//printf("DB new xyz>> %lf %lf %lf\n", p_atom->x,p_atom->y,p_atom->z);
p_atom=p_molecule+iatom;


}
return;
}

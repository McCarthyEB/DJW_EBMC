/*******************************************************/
/***** calculate_box.c : Calculates the size of a*******/
/*****                   box which fits inside a *******/
/*****                   molecules extents       *******/
/***** Dewi 17/5/95                            *********/
/*******************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

void calculate_box(atom *p_molecule, int num_atoms, double *p_box_limits)
{
#include "header.h"

int iloop,i;
atom *p_atom;
double xrange,yrange,zrange;

printf("calculate_box\n");
for (i=0;i<=5;i++) printf("%f \n",p_box_limits[i]);

for (iloop=0 ; iloop <= num_atoms; iloop++)
	{
	p_atom = p_molecule + iloop;

	p_box_limits[0] = MIN(p_box_limits[0], p_atom->x);
	p_box_limits[1] = MAX(p_box_limits[1], p_atom->x);
	p_box_limits[2] = MIN(p_box_limits[2], p_atom->y);
	p_box_limits[3] = MAX(p_box_limits[3], p_atom->y);
	p_box_limits[4] = MIN(p_box_limits[4], p_atom->z);
	p_box_limits[5] = MAX(p_box_limits[5], p_atom->z);
	}

for (i=0;i<=5;i++) printf("%f \n",p_box_limits[i]);

/****calculate how much to lop off each end******/
xrange= (p_box_limits[1] - p_box_limits[0]) * ((1-box_fraction)/2);
yrange= (p_box_limits[3] - p_box_limits[2]) * ((1-box_fraction)/2);
zrange= (p_box_limits[5] - p_box_limits[4]) * ((1-box_fraction)/2);
printf("ranges : %f %f %f %f %f\n",xrange,yrange,zrange,box_fraction,1-box_fraction);
p_box_limits[0] = p_box_limits[0] + xrange;
p_box_limits[1] = p_box_limits[1] - xrange;

p_box_limits[2] = p_box_limits[2] + yrange;
p_box_limits[3] = p_box_limits[3] - yrange;

p_box_limits[4] = p_box_limits[4] + zrange;
p_box_limits[5] = p_box_limits[5] - zrange;

for (i=0;i<=5;i++) printf("%f \n",p_box_limits[i]);

return;
}

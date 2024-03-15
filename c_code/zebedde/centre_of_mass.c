/**************************************************************************/
/* work out the centre of mass for a molecule or list of atoms            */
/* Adapted March 07 to allow atom mol index testing for doing molecules   */
/* as parts of the list.   Dave Willock                                   */
/* Adapted March 09 to ensure that all atoms are at min_image of first    */
/* in list before doing cofm. Dave Willock                                */
/* Corrected July 09 to allow for molecules large with respect to half    */
/* shortest pbc distance. Dave Willock                                    */
/**************************************************************************/
/**************************************************************************/
#include <stdio.h>
#include "maxima.h"
#include "structures.h"

/* protype list for this routine */

double atomic_mass_list( char *element );

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol);

void min_image( double *x, double *y, double *z);
/*---------------------------------------------------------------------------*/

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol )
{
  int icomp,this_atom;
  int ifirst;
  atom *p_atom, *p_first;

  double atomic_mass, x,y,z;
/*----------------------------------------------------------------------------*/

   *p_c_of_m= 0.0;
   *(p_c_of_m+1)= 0.0; 
   *(p_c_of_m+2)= 0.0;

/* assign atomic masses and locate mass weighted total positon */

   *p_total_mass= 0.0;

/*** deal with min_image                             ***/
/*** Gather molecule routine added July 09           ***/
/*** uses neighbour tree to shift atoms to min_image ***/
/*** of their neighbours in sequence so that large   ***/
/*** molecules are shifted correctly.                ***/
/*** Routine added by Dave Willock                   ***/

//   printf("DEBUG>> Calling gather_molecule from centre of mass\n");
   gather_molecule(p_molecule, num_atoms, -1);

/*** For negative values of which_mol do cofm for whole list ****/

   if (which_mol < 0)
     {
       for (this_atom=0; this_atom <= num_atoms; this_atom++)
         {
            p_atom= p_molecule+ this_atom;

            atomic_mass = atomic_mass_list( &(p_atom->elem[0]) );

            *p_c_of_m     += atomic_mass * (p_atom->x);
            *(p_c_of_m+1) += atomic_mass * (p_atom->y);
            *(p_c_of_m+2) += atomic_mass * (p_atom->z);

            *p_total_mass += atomic_mass;
        }
     }

/*** For postive values of which_mol do cofm only for atoms ***/
/*** belonging to that molecule                             ***/

   else
     {
       for (this_atom=0; this_atom <= num_atoms; this_atom++)
         {
            p_atom= p_molecule+ this_atom;
 
            if (p_atom->mol == which_mol)
              {
                atomic_mass = atomic_mass_list( &(p_atom->elem[0]) );

                *p_c_of_m     += atomic_mass * (p_atom->x);
                *(p_c_of_m+1) += atomic_mass * (p_atom->y);
                *(p_c_of_m+2) += atomic_mass * (p_atom->z);

                *p_total_mass += atomic_mass;
             }
        }
     }

/* normalise to total mass */

   for (icomp=0; icomp < 3; icomp++)
                        *(p_c_of_m+icomp)= *(p_c_of_m+icomp) / *p_total_mass;

}

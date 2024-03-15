/*******************************************************************/
/******** sort_out_molecule:                                       */
/********                generates neighbours and hydrogen indexes */
/********                calls centre of mass and potential        */
/********                assignment routines                       */
/********                 djw/dwl 28/2                             */
/******** Addition to look for multiple molecules. March 2007 DJW  */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void find_hydrogens(atom *p_molecule, int *p_hydrogen_list,
                    int *p_num_hydrogens, int num_atoms,
                    int *p_have_AB);

void generate_neighbours( atom *p_molecule, int num_atoms,
                          atom_number *p_types, int *p_num_types,
                          int use_pbc);

void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, 
                    atom *p_molecule, int num_atoms, int which_mol );

void assign_stretch(atom *p_molecule, int num_atoms);

void assign_nonbond(atom *p_molecule, int num_atoms);

void assign_hbond(atom *p_molecule, int num_atoms);

void assign_bond_inc(atom *p_molecule, int num_atoms);

void steric_setup(atom *p_molecule, int num_atoms);

void min_image( double *x, double *y, double *z);

int find_mol(atom *p_molecule, int num_atoms);

void demarcate_molecule(atom *p_molecule, int num_atoms, list_partition *p_demarc);

void sort_out_molecule(atom *p_molecule, int num_atoms, atom_number *p_types,
		       list_partition *p_types_demarc, int *p_hyd_list, int *p_num_hyds, 
                       int use_pbc, vec *p_c_of_m, double *p_total_mass, 
                       int is_pore, int look_for_tethers, int *p_have_AB,
                       int *p_num_mols, list_partition *p_demarc)
{
#include "header.h"

int iatom, imol, found_tether;
double dx, dy, dz;

list_partition *p_this_demarc;
list_partition *p_this_types_demarc;

atom *p_atom, *p_this_mol;

atom_number *p_this_types;

/******************************************************************/
/******** generate neighbour list for the template ****************/
/******************************************************************/

if (num_atoms >0 )
    {
	  generate_neighbours( p_molecule, num_atoms,
                               p_types, &(p_types_demarc->num), use_pbc);

          *p_num_mols = find_mol(p_molecule, num_atoms);

          printf("This structure actually contains %d molecules\n", *p_num_mols);

          if (*p_num_mols >= MAX_MOLS)
            {
              printf("This number exceeds the current maximum of %d, increase MAX_MOLS in maxima.h\n", MAX_MOLS);
              exit(0);
            }

          if (*p_num_mols > 1)
            {
/*** Multiple molecules detected case ****/
               printf("Going to demarcate molecule\n");
               demarcate_molecule( p_molecule, num_atoms, p_demarc);

               p_this_demarc= p_demarc;
               p_this_types_demarc=p_types_demarc;
               p_this_types_demarc->start = 0;

               for (imol=0; imol < *p_num_mols; imol++)
                 {
                    p_this_mol=p_molecule+p_this_demarc->start;
                    p_this_types=p_types+p_this_types_demarc->start;

                    generate_neighbours( p_this_mol, p_this_demarc->num-1,
                                         p_this_types, &(p_this_types_demarc->num), use_pbc);

/****************************************************************************/
/********* Find centre of Mass **********************************************/
/****************************************************************************/
                    centre_of_mass(&(p_c_of_m->v[0]), p_total_mass, p_this_mol, (p_this_demarc->num-1), -1); 

                    p_c_of_m++;
                    p_total_mass++;

                    p_this_types_demarc->end = p_this_types_demarc->start+p_this_types_demarc->num-1;

                    p_this_types_demarc++;
                    p_this_types_demarc->start = (p_this_types_demarc-1)->end + 1;

                    p_this_demarc++;
                 }
            }
          else
            {
/*** Single molecule detected case ****/
               printf("Dealing with as a single molecule\n");

               p_demarc->start=0;
               p_demarc->num=num_atoms+1;
               p_demarc->end=num_atoms;

               p_types_demarc->start = 0;
               p_types_demarc->end   = p_types_demarc->num-1;
/****************************************************************************/
/********* Find centre of Mass **********************************************/
/****************************************************************************/

               centre_of_mass(&(p_c_of_m->v[0]), p_total_mass, p_molecule, num_atoms, -1);
               
            }

          printf("Going to steric_setup\n");

          steric_setup(p_molecule, num_atoms);


          printf("Back from steric_setup\n");

          if (is_pore) printf("this is the pore\n");
          else printf("this is NOT the pore\n");

          if (use_pbc) printf("using pbc\n");
          else printf("NOT using pbc\n");
/****************************************************************************/
/*** If this is not the pore force all atom co-ordinates to their min_image */
/*** with the first in the list                                             */
/*** Update added Feb 06 Dave Willock                                       */
/****************************************************************************/

          if (!is_pore && use_pbc)
            {
              printf("Moving all atoms to same unit cell\n");

              p_atom  = p_molecule;
              printf("Atom 0 %s now at %10.6f %10.6f %10.6f\n",
                                  p_atom->label, p_atom->x,
                                  p_atom->y, p_atom->z);

              for (iatom=1; iatom <= num_atoms; iatom++)
                {
                   p_atom++;
                   dx= p_atom->x - p_molecule->x;
                   dy= p_atom->y - p_molecule->y;
                   dz= p_atom->z - p_molecule->z;
                  
                   min_image( &dx, &dy, &dz);

                   p_atom->x = p_molecule->x + dx;
                   p_atom->y = p_molecule->y + dy;
                   p_atom->z = p_molecule->z + dz;

                   printf("Atom %d %s now at %10.6f %10.6f %10.6f\n",
                                  iatom, p_atom->label, p_atom->x,
                                  p_atom->y, p_atom->z);
                }
            }

/****************************************************************************/
/********* Assign inter molecular potentials ********************************/
/****************************************************************************/

          if (non_bonded) 
             {
                printf("assigning non-bonds\n");
                assign_nonbond( p_molecule, num_atoms);
                printf("assigning h-bonds\n");
                assign_hbond(p_molecule, num_atoms);
                printf("assigning bond incs\n");
                assign_bond_inc(p_molecule, num_atoms);
             }


          if (!is_pore)
            { 

/****************************************************************************/
/************** Find hydrogen atom indexes for template *********************/
/****************************************************************************/

               printf("Finding hydrogens\n");
               find_hydrogens( p_molecule, p_hyd_list, p_num_hyds, num_atoms,
                               p_have_AB);

               printf("done --- Finding hydrogens\n");
/****************************************************************************/
/************** Find atom indicies for tethering in the pore ****************/
/****************************************************************************/
            
               if (have_tethers && look_for_tethers)
                  {
                     found_tether= FALSE;
                     p_atom= p_molecule;

                     for (iatom=0; iatom <= num_atoms; iatom++)
                        {
                           if (strcmp(p_atom->label, tether.B) == 0)
                             {
                                found_tether= TRUE;
                                tether.index_b= iatom;
                                break;
                             }
                           p_atom++;
                        }

                     if (!found_tether) 
                        {
                           fprintf(output_fp, "ERROR>> failed to find molecule tethering atom\n");
                           exit(EXIT_FAILURE);
                        }
                  }

/****************************************************************************/
/********* Assign intra molecular potentials ********************************/
/****************************************************************************/

           if (!steric)
              {
               printf("Assign stretch\n");
               assign_stretch(p_molecule, num_atoms);
              }
            }

/****************************************************************************/
/************** Find atom indicies for tethering in the seed ****************/
/****************************************************************************/

         else
            {
               if (have_tethers && look_for_tethers)
                  {
                     found_tether= FALSE;
                     p_atom= p_molecule;
                     for (iatom=0; iatom <= num_atoms; iatom++)
                        {
                           if (strcmp(p_atom->label, tether.A) == 0)
                             {
                                found_tether= TRUE;
                                tether.index_a= iatom;
                                break;
                             }
                           p_atom++;
                        }

                     if (!found_tether)
                        {
                           fprintf(output_fp, "ERROR>> Failed to find molecule tethering atom\n");
                           exit(EXIT_FAILURE);
                        }

                  }
            }
    }
else
    {
	printf("ERROR: No atoms to sort out in sort_out_molecule\n");
    }

printf("Returning from sort_out_molecule\n");
return;
}


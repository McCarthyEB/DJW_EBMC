#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

/* routine to apply symmetry operations to sets of atoms */
/* Updated December 2010 for array of pointers structure */
/* this routine now recieves the reference molecule      */
/* and a pointer to the start of the array for the symm  */
/* image. Makes a single symm image. Dave & Filippo      */

/* required functions ------------------------------------*/

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void matrix_vec_mult( double *p_matrix, double *p_vector, double *p_answer);

void min_image( double *x, double *y, double *z);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol);
/*--------------------------------------------------------*/

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm )
{
#include "header.h"
 
  atom *p_atom, *p_new_atom;

  double base_coords[3], new_coords[3];
  double want_pos[3], cell_centre[3];

  int iatom;

/**** Make sure all the atoms of the molecule are in the same unit cell ****/

 gather_molecule(p_molecule, num_atoms, -1);

/****** loop over each atom in the starting molecule **************/
/****** Debugged Aug 09 Dave Willock, only use fractional *********/
/****** co-ords for the translation component of the operation ****/
/****** earlier version probably failed for non-orthogonal cells **/

/****** Loop over atoms **********/
     for (iatom= 0; iatom <= num_atoms; iatom++)
        {
          p_atom= p_molecule+iatom;
          p_new_atom= p_image +iatom;

/*** copy identity of atom to the new image ***/
          *p_new_atom= *p_atom;

/*** copy co-ordinates for application of symmetry *****/
          base_coords[0]= p_atom->x;
          base_coords[1]= p_atom->y;
          base_coords[2]= p_atom->z;

/******* apply the symmetry matrix ********************************/
/******* result held in new_coords ********************************/

          matrix_vec_mult( &(p_symm->matrix[0]), &base_coords[0], &new_coords[0]);

/******* apply the symmetry translation ********************************/
/******* Note that the symmetry translation refers to fractional       */
/******* co-ordinates so the new_coords after matrix operation are     */
/******* converted to fractional, held in base_coords, the translation */
/******* is applied and then the atom updated.                         */
/***********************************************************************/

      if (pbc)
        {
           cart_to_fract( new_coords[0], new_coords[1], new_coords[2],
                          &base_coords[0], &base_coords[1], &base_coords[2],
                          &recip_latt_vec[0]);

          new_coords[0] = base_coords[0] + p_symm->translation[0]; 
          new_coords[1] = base_coords[1] + p_symm->translation[1]; 
          new_coords[2] = base_coords[2] + p_symm->translation[2]; 

/*** update co-ordinates of new atom **********/
          fract_to_cart( &(p_new_atom->x), &(p_new_atom->y),
                         &(p_new_atom->z), new_coords[0], 
                         new_coords[1], new_coords[2],
                         &latt_vec[0]);

            }
          else
            {
               p_new_atom->x= new_coords[0];
               p_new_atom->y= new_coords[1];
               p_new_atom->z= new_coords[2];
            }
        }

/**********************************************************************/
/**** Move each molecule back to its minimum image position     *******/
/**** Note that the "which_mol = -1" used here works because    *******/
/**** the atom list is restricted to one symmetry image.        *******/
/**** Assuming there is only one molecule when symmetry is on ! *******/
/**********************************************************************/

/*** find centre of cell **********************/
          fract_to_cart( &cell_centre[0], &cell_centre[1],
                         &cell_centre[2], 0.5, 0.5, 0.5,
                         &latt_vec[0]);

     for (iatom= 0; iatom <= num_atoms; iatom++)
        {
          p_new_atom= p_image +iatom;

          want_pos[0]= p_new_atom->x - cell_centre[0];  
          want_pos[1]= p_new_atom->y - cell_centre[1];  
          want_pos[2]= p_new_atom->z - cell_centre[2];   

          min_image( &want_pos[0], &want_pos[1], &want_pos[2]);   

          p_new_atom->x= cell_centre[0] + want_pos[0];  
          p_new_atom->y= cell_centre[1] + want_pos[1];  
          p_new_atom->z= cell_centre[2] + want_pos[2];  
       }

       gather_molecule(p_image, num_atoms, -1);

 return;
}


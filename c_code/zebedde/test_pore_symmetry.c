/*********************************************************************/
/* test_pore_symmetry.c : check that the symmetry operations supplied*/
/*                        work for the pore and so are valid to use. */
/*                        moved out of make_a_template               */
/*                        Feb. 2013 Dave Willock.                    */
/*********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm );

void min_image( double *x, double *y, double *z);

int test_pore_symmetry(FILE *fp, atom *p_pore, int num_atoms, 
                       symm_ops *p_symm, int num_symm_ops)
{
int isymm, iii, jjj, kkk; 
int *p_flags, *p_this_flag;
double *vec_comp, min_dist, max_image_dist;
double dist2;
double dx, dy, dz;
atom *pore_image;
atom *p_this_atom, *p_this_image;

/**** First test that the symmetry operation are valid for the pore ****/
       fprintf(fp,"Testing proposed operations on pore\n\n");

/**** malloc the image array to make space for copy ***/

       if (num_atoms >= 0)
         {
            printf("Mallocing image array for pore with size: %d\n", num_atoms+1);
            pore_image=(atom*)malloc((num_atoms+1)*sizeof(atom));

            if (pore_image == NULL)
               {
                  printf("ERROR: Problem with memory assignment for pore symmetry image %d entries long.\n",
                                                              num_atoms+1);
                  exit(0);
               }

            printf("Mallocing flags array for pore with size: %d\n", num_atoms+1);
            p_flags=(int*)malloc((num_atoms+1)*sizeof(int));

            if (p_flags == NULL)
               {
                  printf("ERROR: Problem with memory assignment for flags array for pore symmetry test %d entries long.\n",
                                                              num_atoms+1);
                  exit(0);
               }
         }
       else
         {
            printf("WARNING: Cannot do symmetry on pore as it is empty!\n");
            return TRUE;
         }

/*** Now check if symmetry operations work ***/
       if (num_atoms < 0 )
          {
             fprintf(fp, "Skipping test as pore array is too small...\n");
          }
       else
          {
             for (isymm=0; isymm <= num_symm_ops; isymm++)
                {
                    fprintf(fp,"Image %d\n", isymm+1);
                    printf("Image %d\n", isymm+1);

/**** Pass symm_ops one at a time for the pore  ************************/

                    do_symmetry( p_pore, num_atoms,
                                 pore_image, p_symm );

                    p_symm++;

/**** Compare transformed and original lattice *************************/
/**** could speed up for big pores by flagging matches *****************/

                    p_this_flag=p_flags;
                    for (jjj=0; jjj <= num_atoms; jjj++) 
                      {
                        *p_this_flag=FALSE;
                        p_this_flag++;
                      }

/*** loop over atoms ***/
                    p_this_atom=p_pore;
                    for (jjj=0; jjj <= num_atoms; jjj++)
                       {
/*** loop over images ***/
                         p_this_image=pore_image;
                         for (kkk= 0; kkk<=num_atoms; kkk++)
                            {
                              if (!*p_this_flag)
                                {
                                  dx = p_this_atom->x - p_this_image->x;
                                  dy = p_this_atom->y - p_this_image->y;
                                  dz = p_this_atom->z - p_this_image->z;
                                  p_this_image++;

                                  min_image(&dx, &dy, &dz);

                                  dist2= dx*dx + dy*dy + dz*dz;
                                  if (kkk==0 || dist2 < min_dist) min_dist = dist2;

                                  if (dist2 < 0.00001) *p_this_flag=TRUE;
                                }
                            }

                         if (jjj==0 || min_dist > max_image_dist) max_image_dist = min_dist;
                         p_this_atom++;
                       }

                     if ( max_image_dist < 0.00001 )
                       {
                         fprintf(fp,"Pore obeys this symmetry operation,");
                         fprintf(fp," greatest difference between original and transformed lattice : %10.6f\n",
                                                   sqrt(max_image_dist));
                       }
                    else
                       {
                         fprintf(fp,"ERROR: Pore does not obey this symmetry operation,");
                         fprintf(fp," greatest difference between original and transformed lattice : %10.6f\n",
                                                   sqrt(max_image_dist));
                         exit(0);
                       }
                }
          }

/**** Free up memory ****/
printf("Freeing memory for pore_image and related flags\n");
free(pore_image); free(p_flags);

return TRUE;
}

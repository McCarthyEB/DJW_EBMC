/************************************************************/
/* calculate_intra_energy.c                                 */
/* Calculates the self interaction energy for a molecule    */
/*                                                          */
/* Started DWL 27/11/94                                     */
/* Updated for angle, torsion and vdw terms Sept. 06 DJW    */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

void min_image( double *x, double *y, double *z);

double non_bond_r_eps_sixth(atom *p_pore, int num_p_atoms,
                            atom *p_molecule, int num_t_atoms, int *p_need_grad,
                            double *p_grad, int is_intra, 
                            vdw_interact_list *p_vdw_list, int num_vdws_list,
                            double *p_pos_separations, double *p_pos_vectors);

double angle_bend_energy(atom *p_molecule, int num_atoms, 
                         angle_interact_list *p_angles_list, int num_angles_listed);

double torsion_energy(atom *p_molecule, int num_atoms, 
                      torsion_interact_list *p_torsions_list, int num_torsions_listed);

void calculate_intra_energy( atom *guest_ptrs[], int num_guests,
                             list_partition *p_guest_demarc,
                             angle_interact_list *p_angles_list_ptrs[],
                             int *p_num_angles_list, 
                             torsion_interact_list *p_torsions_list_ptrs[],
                             int *p_num_torsions_list, 
                             vdw_interact_list *p_vdw_list_ptrs[],
                             int *p_num_vdws_list)
{

#include "header.h"
int iatom,jatom, neigh_index, pot_index;
int imol, ind;
double dx,dy,dz,r, dummy[2];
double r2,r3,r4, exp_sqrd;
double *dum_seps, *dum_vecs;

atom *p_atom, *p_neigh;

/*DEBUG=TRUE;*/

// printf("Arrived in intra_energy with %d guests\n", num_guests);

/*******************************************************************************/
/******** Loop over all atoms in the molecule using neighbours info to *********/
/******** pick out atom pairs                                          *********/
/*******************************************************************************/

 for (imol=0; imol<num_guests; imol++)
   {
     if (symm_set) ind=imol*(num_symm_ops+2);
                                             else ind=imol;

     intra_energy[imol].stretch = 0.0;

     p_atom= guest_ptrs[ind];
     for (iatom=0; iatom <= p_guest_demarc->end; iatom++)
       {
//       if (DEBUG) printf("Atom %d looking at neighbours of %s (pot= %s)\n", p_atom->label,
//                                                                            p_atom->pot);


         for (neigh_index=0; neigh_index <= p_atom->num_neigh ; neigh_index++)
          {
            jatom = p_atom->neighb[neigh_index];

/*******************************************************************************/
/******** Avoid double counting by only allowing lower indicies ****************/
/*******************************************************************************/

            if (jatom > iatom)
               {
                  p_neigh = guest_ptrs[ind] + jatom;
                  pot_index = p_atom->neighb_stretch_list[neigh_index];

                  dx = p_atom->x - p_neigh->x;
                  dy = p_atom->y - p_neigh->y;
                  dz = p_atom->z - p_neigh->z;
                  
                  if (pbc) min_image( &dx, &dy, &dz);

                  r = dx*dx + dy*dy + dz*dz;
                  r = sqrt(r);

//    if (DEBUG) printf("%s and %s are %10.6f apart\n", p_atom->label, p_neigh->label, r);

                  if (intra_pair_potent[pot_index].which == QUARTIC_STRETCH)
                     {

/*************************************************************************/
/****** Calculate this ones contribution to the stretch assuming    ******/
/****** we are using cff91 quartic stretch potential:               ******/
/****** E = K2 * (R - R0)^2  +  K3 * (R - R0)^3  +  K4 * (R - R0)^4 ******/
/*************************************************************************/

                         r = r - intra_pair_potent[pot_index].stretch_A;

                         r2= r*r;
                         r3= r*r2;
                         r4= r2*r2;               
              
                         intra_energy[imol].stretch += intra_pair_potent[pot_index].stretch_B*r2
                                                     + intra_pair_potent[pot_index].stretch_C*r3
                                                     + intra_pair_potent[pot_index].stretch_D*r4;

 //                    if (DEBUG)
 //                      {
 //                          printf("Using QUARTIC STRETCH\n");
 //                          printf("reference bond length: %10.6f, contribution %10.6f\n",
 //                                     intra_pair_potent[pot_index].stretch_A,
 //                                     intra_pair_potent[pot_index].stretch_B*r2
 //                                   + intra_pair_potent[pot_index].stretch_C*r3
 //                                   + intra_pair_potent[pot_index].stretch_D*r4);
 //                      }

                     }

/*********************************************************************/
/*** Morse potential calculation with form: **************************/
/*** E = D * (1 - exp(-ALPHA*(R - R0)))^2   **************************/
/*********************************************************************/

                  else if (intra_pair_potent[pot_index].which == MORSE_STRETCH)
                     {
                         r = r - intra_pair_potent[pot_index].stretch_A;
                
                         exp_sqrd= 1-exp(-intra_pair_potent[pot_index].stretch_C*r);
                         exp_sqrd= exp_sqrd*exp_sqrd;

                         intra_energy[imol].stretch += intra_pair_potent[pot_index].stretch_B*exp_sqrd; 
                     }

/*********************************************************************/
/*** Quadratic bond calculation with form:  **************************/
/*** E = K2 * (R - R0)^2                    **************************/
/*********************************************************************/

                  else if (intra_pair_potent[pot_index].which == QUADRATIC_STRETCH)
                     {
                         r = r - intra_pair_potent[pot_index].stretch_A;

                         r2= r*r;
              
                         intra_energy[imol].stretch += intra_pair_potent[pot_index].stretch_B*r2;

//                     if (DEBUG)
//                       {
//                           printf("Using QUADRATIC STRETCH\n");
//                           printf("reference bond length: %10.6f, contribution %10.6f\n",
//                                      intra_pair_potent[pot_index].stretch_A,
//                                      intra_pair_potent[pot_index].stretch_B*r2);
///////                  }
                     }
                  
               }
          }
//    if (DEBUG) printf("Stretch now %10.6f\n", intra_energy[imol].stretch);
      p_atom++;
    }
/*********************************************************************/
/** Can now do angle potential contribution as a sum over the angles */
/** list. Put as a new subroutine!                                   */
/*********************************************************************/

//  printf("stretches done..\n");
//  printf("....now for %d angles for molecule %d index %d\n", *p_num_angles_list, imol, ind);
  if (*p_num_angles_list >= 0)
    {
      intra_energy[imol].angle= angle_bend_energy(guest_ptrs[ind], p_guest_demarc->end, 
                                                  p_angles_list_ptrs[imol], *p_num_angles_list);
    }
  else
      intra_energy[imol].angle= 0.0;

//  printf("angles done...now for %d torsions\n", *p_num_torsions_list);
  if (*p_num_torsions_list >= 0)
    {
      intra_energy[imol].torsion= torsion_energy(guest_ptrs[ind], p_guest_demarc->end, 
                                                 p_torsions_list_ptrs[imol], *p_num_torsions_list);
    }
  else
      intra_energy[imol].torsion= 0.0;

/**** TESTING with pcff forcefield, need to make decision on which routine to ****/
/**** call eventually!                                                        ****/

//  printf("torsions done...now for %d vdws\n", *p_num_vdws_list);
  if (*p_num_vdws_list >= 0)
    {
       intra_energy[imol].vdw = non_bond_r_eps_sixth(guest_ptrs[ind], p_guest_demarc->end, 
                                                     guest_ptrs[ind], p_guest_demarc->end, 
                                                     FALSE, &dummy[0], TRUE,
                                                     p_vdw_list_ptrs[imol], *p_num_vdws_list, 
                                                     dum_seps, dum_vecs);
    }
  else
      intra_energy[imol].vdw = 0.0;

//  printf("vdws done .... now summing\n");
  intra_energy[imol].total = intra_energy[imol].stretch + intra_energy[imol].angle
                            +intra_energy[imol].torsion + intra_energy[imol].vdw;

  p_num_angles_list++;
  p_num_torsions_list++;
  p_num_vdws_list++;

  p_guest_demarc++;
 }
return;
}

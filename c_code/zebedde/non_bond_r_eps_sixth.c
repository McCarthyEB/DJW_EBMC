/**************************************************************************/
/******* Non-bond energy calculation assuming r-eps potential *************/
/******* format and sixth power combination rules             *************/
/******* as used by cff91 and related forcefields             *************/
/******* Dave Willock started 2nd December 1995               *************/
/*******              added list processing May 1996          *************/
/*******              added intra version for calculating     *************/
/******* self interaction (not periodic images) contribution  *************/
/******* to intra moelcular energy. Using pre-assembled list  *************/ 
/******* of atom pairs which has any cases that should be     *************/
/******* ignored (e.g. 1-3 interactions in angles) excluded.  *************/
/******* This term appears in the atoms own vdw record and is *************/
/******* held in intra_energy as well.                        *************/
/*******                           August 2006.               *************/
/******* Added suppression of multiple warnings about close   *************/
/******* atoms via "told" flag.                               *************/
/*******                           Dave Willock March 2007    *************/
/******* Changed to a global told_of_clash logical to give    *************/
/******* and then suppress warnings to one per energy calc.   *************/
/******* and a local clash logical for catching clashes.      *************/
/*******                           Dave Willock July 2009     *************/
/******* For multi-molecule cases added "is_self" flag so     *************/
/******* the molecule with its own periodic images can be     *************/
/******* included.                 Dave Willock, May 2013     *************/
/**************************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

void min_image( double *x, double *y, double *z);

int pbc_interactions(double *p_dx, double *p_dy, double *p_dz, double cutoff,
                     double cutoff_2, double *p_pos_separations,  double *p_pos_vectors, 
                     int count_flag, int is_self);

double non_bond_r_eps_sixth(atom *p_pore, int num_p_atoms,
                            atom *p_molecule, int num_t_atoms, int *p_need_grad,
                            double *p_grad, int is_intra,
                            vdw_interact_list *p_vdw_list, int num_vdws_listed,
                            double *p_pos_separations, double *p_pos_vectors,  
                            int is_self)
{
#include "header.h"
atom *p_tatom, *p_patom, *p_atom1, *p_atom2;

int j,pot_index_t, pot_index_p,this_pair;
int index_temp_atom,num_interactions; 
int pot_index_1, pot_index_2;
int ilist, clash;

double dx,dy,dz,r0,r06,r09,r,r2,r3,r6,r8,r9, eps, sum_sixth;
double nb_energy, nb_energy_rep, nb_energy_disp;
/* Replace parameter for array size with a malloc in calling routine (energy.c) JL Oct. 2007*/
/*double pos_separations[MAX_PAIR_LIST], pos_vectors[3*MAX_PAIR_LIST];*/
double rep,disp,grad_factor;

double *p_grad_x, *p_grad_y, *p_grad_z;

nb_energy_rep = 0.0;
nb_energy_disp = 0.0;

//  printf("In non_bond_r_eps_sixth with %d pore atoms %d molecule atoms\n", num_p_atoms, num_t_atoms);
//  printf("%d interactions in vdws list\n", num_vdws_listed);

if (is_intra)
  {
//    printf("This is an intra molecular calculation\n");
/**************************************************************************/
/********** loop over molecule atoms :                 ********************/
/********** p_tatom = pointer to current molecule atom ********************/
/**************************************************************************/

    clash= FALSE;

    for (ilist=0; ilist<=num_vdws_listed; ilist++)
      {
        p_atom1= p_molecule + p_vdw_list->iatm1;
        p_atom2= p_molecule + p_vdw_list->iatm2;
        p_vdw_list++;

        pot_index_1= p_atom1->nb_list;
        pot_index_2= p_atom2->nb_list;
 
        dx = p_atom1->x - p_atom2->x;
        dy = p_atom1->y - p_atom2->y;
        dz = p_atom1->z - p_atom2->z;

        min_image( &dx, &dy, &dz);
        
        r2 = dx*dx + dy*dy + dz*dz;

/*************************************************************************/
/******* Apply combining rules to the atom potent parameters   ***********/
/******* using the nb_list parameter each atom should have for ***********/
/******* referencing the potent list                           ***********/
/******* When reading r-eps potentials a= r; b= eps            ***********/
/*************************************************************************/

        sum_sixth = potent[pot_index_1].a6 + potent[pot_index_2].a6;

        eps = 2.0 * potent[pot_index_1].sqrt_b * 
                    potent[pot_index_2].sqrt_b *
                    potent[pot_index_1].a3 *
                    potent[pot_index_2].a3 / sum_sixth;

        r06 = sum_sixth / 2.0;

        r0 = pow(r06,(1.0/6.0));
 
        r09 = r0 * r0 * r0 * r06;

/*************************************************************************/

       rep = 0;
       disp= 0;

       r = sqrt(r2);

/******* Avoid very close contacts and report ****************************/

       if (r <= 0.1) 
         {
            if (!told_of_clash)
              { 
                printf("Warning two intra-molecular atoms are at a separation less than 0.1A,\n");
                printf("Have capped their contribution to the energy at 10000.");
                printf("Multiple warnings suppressed until next build success or accepted move if MC\n");
                told_of_clash=TRUE;
              }
            nb_energy_rep       += 10000;
            clash = TRUE;

/*** No use carrying on if clash present ***/
/*** Send back high energy for rejection ***/

            if (clash) return nb_energy_rep;
          }
        else
          {
            r3 = r2*r;
            r6 = r3 * r3;
            r9 = r6 * r3;
                
/******* Repulsion part of potential *************************************/

            rep = 2.0 * eps * r09/r9;

/******* Dispersion part of potential ************************************/

            disp= 3.0 * eps * r06/r6;

            nb_energy_rep       += rep; 
            nb_energy_disp      -= disp;

          }
      }
    return nb_energy_rep+nb_energy_disp;
  }

clash= FALSE;

/**************************************************************************/
/********** loop over molecule atoms :                 ********************/
/********** p_tatom = pointer to current molecule atom ********************/
/**************************************************************************/

//printf("DEBUG>> host %d atoms with guest %d atoms...\n", num_p_atoms,num_t_atoms); 

nb_energy_rep = 0.0;
nb_energy_disp = 0.0;

p_tatom = p_molecule-1;     
for (index_temp_atom=0; index_temp_atom<=num_t_atoms ;index_temp_atom++)
   {
     p_tatom++;
     pot_index_t = p_tatom->nb_list;

     if (*p_need_grad)
       {
         p_grad_x = p_grad+ 3*index_temp_atom;
         p_grad_y = p_grad_x+1;
         p_grad_z = p_grad_y+1;
       }

/**************************************************************************/
/********** loop over pore atoms :                     ********************/
/********** p_patom = pointer to current pore atom     ********************/
/**************************************************************************/

     p_patom = p_pore-1;     

     for ( j=0; j<=num_p_atoms; j++)
        {
          p_patom++;
          pot_index_p = p_patom->nb_list;

//          printf("DB>> Interaction pair %s (%d) with %s (%d)\n",
//                  p_tatom->label, pot_index_t, p_patom->label, pot_index_p); 

          dx = p_tatom->x - p_patom->x;
          dy = p_tatom->y - p_patom->y;
          dz = p_tatom->z - p_patom->z;

//          printf("\nAtom %d at %10.6f  %10.6f  %10.6f \n", index_temp_atom, p_tatom->x,
//                                                                            p_tatom->y,
//                                                                            p_tatom->z);
//
//          printf("Atom %d at %10.6f  %10.6f  %10.6f \n\n", j, p_patom->x,
//                                                              p_patom->y,
//                                                              p_patom->z);

/**************************************************************************/
/******* if this is a periodic pore assemble a list of interactions *******/
/******* of this type otherwise check if the pair is in range       *******/
/**************************************************************************/

          if (pbc)
            {
               num_interactions = pbc_interactions( &dx, &dy, &dz, nb_ctf, nb_ctf_2,  
                                                              p_pos_separations, p_pos_vectors, 
                                                              FALSE, is_self);

               if (DEBUG && num_interactions >= 0 && *p_pos_separations < 2.0) 
                   {
                       printf("vector sent: %10.6f %10.6f %10.6f\n", dx, dy, dz);

                       printf("%s %10.6f %10.6f %10.6f %s %10.6f %10.6f %10.6f, r= %10.6f\n",
                                      p_tatom->label, p_tatom->x, p_tatom->y, p_tatom->z,
                                      p_patom->label, p_patom->x, p_patom->y, p_patom->z,
                                      *p_pos_separations);
                   }
            }
          else
            {
               r2 = dx*dx + dy*dy + dz*dz;

               if (!is_self && r2 <= nb_ctf_2)
                 {
                    num_interactions=0;
                    p_pos_separations[0]= sqrt(r2);
                    p_pos_vectors[0]= dx;
                    p_pos_vectors[1]= dy;
                    p_pos_vectors[2]= dz;
                 }
               else
                 {
                    num_interactions=-1;
                 }
            }


/*************************************************************************/
/******* Only bother with list if there are some members in it! **********/
/*************************************************************************/

          if (num_interactions >= 0)
            {
//                printf("DEBUG: num_interactions %d\n", num_interactions);

/*************************************************************************/
/******* Apply combining rules to the atom potent parameters   ***********/
/******* using the nb_list parameter each atom should have for ***********/
/******* referencing the potent list                           ***********/
/******* When reading r-eps potentials a= r; b= eps            ***********/
/*************************************************************************/

                sum_sixth = potent[pot_index_t].a6 + potent[pot_index_p].a6;

                eps = 2.0 * potent[pot_index_t].sqrt_b * 
                            potent[pot_index_p].sqrt_b *
                            potent[pot_index_t].a3 *
                            potent[pot_index_p].a3 / sum_sixth;

                r06 = sum_sixth / 2.0;

                r0 = pow(r06,(1.0/6.0));
 
                r09 = r0 * r0 * r0 * r06;

/*************************************************************************/

               rep = 0;
               disp= 0;
               for (this_pair= 0; this_pair <= num_interactions; this_pair++) 
                 {
                    r = *(p_pos_separations+this_pair);

//                    printf("DEBUG>> processing atoms at distance %10.6f\n", r);

/******* Avoid very close contacts and report ****************************/

                    if (r <= 0.1) 
                      {
                        if (!told_of_clash)
                          { 
                            printf("Warning two inter-molecular atoms are at a separation less than 0.1A,\n");
                            printf("Have capped their contribution to the energy at 10000.");
                            printf("Multiple warnings suppressed until next build success or accepted move if MC\n");
                            told_of_clash= TRUE;
                          }

                        clash = TRUE;
                        p_tatom->vdw_energy += 10000;
                        nb_energy_rep       += 10000;
/*** No use carrying on if clash present ***/
/*** Send back high energy for rejection ***/
                   
                        if (clash) return nb_energy_rep;
                      }
                    else
                      {
                        r2 = r*r;
                        r3 = r2*r;
                        r6 = r3 * r3;
                        r9 = r6 * r3;

                
/******* Repulsion part of potential *************************************/

                        rep = 2.0 * eps * r09/r9;

/******* Dispersion part of potential ************************************/

                        disp= 3.0 * eps * r06/r6;

/****** Will calculate the derivatives ***********************************/

                        if (*p_need_grad)
                         {
                            grad_factor= 9.0*rep/r2;
                            *p_grad_x -= grad_factor* p_pos_vectors[3*this_pair]  ;
                            *p_grad_y -= grad_factor* p_pos_vectors[3*this_pair+1];
                            *p_grad_z -= grad_factor* p_pos_vectors[3*this_pair+2];

                            r8 = r6 * r2;
                            grad_factor= 18.0*eps*r06/r8;
                            *p_grad_x += grad_factor* p_pos_vectors[3*this_pair]  ;
                            *p_grad_y += grad_factor* p_pos_vectors[3*this_pair+1];
                            *p_grad_z += grad_factor* p_pos_vectors[3*this_pair+2];
                          }

                        p_tatom->vdw_energy += rep-disp;
                        p_patom->vdw_energy += rep-disp;
                        nb_energy_rep       += rep; 
                        nb_energy_disp      -= disp;

//           printf("Atom %d with atom %d at distance %10.6f rep=%10.6f, disp=%10.6f\n",
//                                 index_temp_atom, j, r, rep, disp);
                     }
                 }
             }
         }
   }

nb_energy= nb_energy_rep + nb_energy_disp;
//printf("DEBUG>> returning energy %10.6f..\n", nb_energy);

return nb_energy;
}

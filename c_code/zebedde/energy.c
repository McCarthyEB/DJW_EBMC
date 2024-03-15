/************************************************************/
/* energy.c                                                 */
/* Calculates the interaction between the pore and template */
/* Mode of calculation determined by input file             */
/*                                                          */
/* Parameters:                                              */
/*  struct pore, struct template                            */
/* Returns:                                                 */
/*  energy of interaction, acceptance flag                  */
/*                                                          */
/* Started DWL 27/11/94                                     */
/* Updated for AMBER forcefield February 1999 DJW           */
/*                                                          */
/* Updated for multi-molecules March 2013 Dave Willock      */
/* Now includes guest-guest interactions in the             */
/* inter-molecular energy, via the molmol list              */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "maxima.h"
#include "ewald.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"
#include "own_maths.h"

double timer(time_t reference_t);

void min_image( double *x, double *y, double *z);

void ewald_sum(atom *p_molecule, int num_molecule_atoms,
               atom *p_crystal_cell, int num_cell_atoms,
               double *p_kvecs, double *p_kvec2,
               double *p_gvec2, int num_kvecs,
               double *p_cos_sum, double *p_sin_sum,
               int zero_first, int self, int using_host);

void isol_coul(atom *p_mol_A, int num_mol_A_atoms,
               atom *p_mol_B, int num_mol_B_atoms, int zero_first, int self);

double non_bond_r_eps_sixth(atom *p_pore, int num_p_atoms,
                            atom *p_molecule, int num_t_atoms, int *p_need_grad,
                            double *p_grad, int is_intra, 
                            vdw_interact_list *p_vdw_list, int num_vdws_listed,
                            double *p_pos_separations, double *p_pos_vectors, int is_self);

double non_bond_r_eps_arithmetic(atom *p_pore, int num_p_atoms,
                                 atom *p_template, int num_t_atoms, int *p_need_grad,
                                 double *p_grad, double *p_pos_separations, double *p_pos_vectors,
                                 int is_self);

double non_bond_a_b_geometric(atom *p_pore, int num_p_atoms,
                              atom *p_template, int num_t_atoms, int *p_need_grad,
                              double *p_grad, double *p_pos_separations, double *p_pos_vectors,
                              int is_self);

double non_bond_a_b_none(atom *p_pore, int num_p_atoms,
                         atom *p_template, int num_t_atoms, int *p_need_grad,
                         double *p_grad, double *p_pos_separations, double *p_pos_vectors,
                         int is_self);

int count_vdw_pbc(atom *p_pore, int num_p_atoms,
                  atom *p_template, int num_t_atoms,
                  int is_self);

void unit_vector(double *p_vector, double *p_size);

void find_line_atoms(int line, atom *p_mol, int num_atoms, int *p_iA, int *p_iB);

double hbond_amber(atom *p_pore, int num_p_atoms,
                   atom *p_template, int num_t_atoms, int *p_need_grad,
                   double *p_grad);

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *guest_ptrs[], int num_guests,
                      list_partition *p_guest_demarc,
                      interaction_indices *p_molmol_list,
                      int num_molmol_list,
                      double *p_kvecs, double *p_kvec2,
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad,
                      int have_comb_rules, int is_empty_pore)
{
#include "header.h"
int i,j, iatom, imol, jmol, index, jdex, isymm, jstart, jsymm;
int ilist, ind, jnd, zero=0;
int zero_first, self_term;
int using_host, list_size;
int pairs_to_malloc, debug_atom;
int not_intra=FALSE;

double dx,dy,dz,r,v_sum, contrib;
double dot, line_now[3], theta;
double gg_inter, size;
double *pos_separations, *pos_vectors;

list_partition *p_demarc, *p_i_demarc, *p_j_demarc;

interaction_indices *p_this_list;

atom *p_patom;
atom *p_tatom;
atom *p_atom;
atom *p_image;
atom *p_A, *p_B;

vdw_interact_list *dum_vdw_list;
int irestr, iA, iB;

//printf("Calculating Energy for %d guests\n", num_guests);

for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
  {
    interaction_energy[imol].steric = FALSE;
    interaction_energy[imol].charges = 0.0;
    interaction_energy[imol].non_bonded = 0.0;
    interaction_energy[imol].hbond = 0.0;
    interaction_energy[imol].restraint = 0.0;
  }

/********************************************************************************/
/****** From just steric information we can calculate the sum of close **********/
/****** contact distances, between guests and pore, only need the      **********/
/****** base molecule of any symmetry related set as the rest must be  **********/
/****** in identical environments.                                     **********/
/********************************************************************************/

if (steric == TRUE && !is_empty_pore) 
  {
    p_demarc=p_guest_demarc;    
    for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
      {
        if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

        for (iatom=0; iatom<= p_demarc->end; iatom++) 	/* loop over this guests atoms */
          {
            closest_contact[iatom]= 100.0;
            p_tatom = &guest_ptrs[index][iatom];          /* pointer to current guest atom */

            for (j=0; j<= num_p_atoms; j++)             /* loop over pore atoms */
              {
                p_patom = p_pore+j;                     /* pointer to current pore atom */

/********************************************************************************/
/********************** calc separation *****************************************/
/********************************************************************************/

               dx = p_tatom->x - p_patom->x;
               dy = p_tatom->y - p_patom->y;
               dz = p_tatom->z - p_patom->z;

/********************************************************************************/
/************ if this is a periodic pore use the nearest image convention *******/
/********************************************************************************/

               if (pbc) min_image( &dx, &dy, &dz);

               r = dx*dx + dy*dy + dz*dz;
               r = sqrt(r);

/********************************************************************************/
/********* update closest contact template -> pore distance *********************/
/********************************************************************************/

               if ( r < closest_contact[iatom] ) closest_contact[iatom]=r;

/********************************************************************************/
/********* calc sum vderWaals radii *********************************************/
/********************************************************************************/

               v_sum = p_patom->vdw + p_tatom->vdw; 

//                  v_sum = v_sum * vdw_scale;              /* do not mult by scaling factor */

               if (r < v_sum) interaction_energy[imol].steric = TRUE; 
            } 
         }
       p_demarc++;
      }
   }

/*=============================================================================*/
/******** end steric ***********************************************************/
/*=============================================================================*/
/*=============================================================================*/
/******** calculate non bond ***************************************************/
/*=============================================================================*/

if (non_bonded) 
  {
/*******************************************************************************/
/******** Decide which potential type we are using and what combination ********/
/******** rules                                                         ********/
/*******************************************************************************/

/*******************************************************************************/
/**** Interaction of each guest with the pore **********************************/
/*******************************************************************************/

    if (!is_empty_pore)
      {
//         printf("Calculating interaction with pore....\n");
         p_demarc=p_guest_demarc;    
         for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
           {
             if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;
    
        if (DEBUG) printf("DB>> Trying to get non-bond energy for molecule %d\n", imol);

/**** Count number of periodic images required and malloc accordingly **********/
/**** for guest with pore                                             **********/

             list_size= count_vdw_pbc(p_pore, num_p_atoms,
                                      guest_ptrs[index], p_demarc->end, FALSE);
                  
             if(list_size < 1)      pairs_to_malloc=100;
                                    else pairs_to_malloc=list_size;

//             printf("cut_off= %10.6f, cut_off2=%10.6f, list_size = %d pairs_to_malloc = %d\n", 
//                             nb_ctf, nb_ctf_2, list_size, pairs_to_malloc);

             pos_separations=(double*)malloc(pairs_to_malloc*sizeof(double));
             pos_vectors=(double*)malloc(3*pairs_to_malloc*sizeof(double));
/**** mallocing finished ****/

/**** zero the atom vdw_energies ***/
            if (imol == 0)
              {
                p_atom= p_pore;
                for (iatom=0; iatom <= num_p_atoms; iatom++)
                  {
                    p_atom->vdw_energy=0;
                    p_atom++;
                  }
              }
           
            p_atom= guest_ptrs[index];
            for (iatom=0; iatom <= p_demarc->end; iatom++)
             {
               p_atom->vdw_energy=0;
               p_atom++;
             }

//      if (DEBUG)
//       {
//        printf("DB>> Energy before doing anything = %10.6f\n",interaction_energy[imol].non_bonded );
//        if (*p_need_grad) printf("DB>> Need Gradient is TRUE\n"); else printf("DB>> Need Gradient is FALSE\n");
//       }

/**** Now work out interaction energies ****/
            if ( strcmp(pot_info.type, R_EPS) == 0 )
              {
                if (strcmp(pot_info.combination, SIXTH_POWER) == 0 )
                  {
                        interaction_energy[imol].non_bonded +=  
                                  non_bond_r_eps_sixth(p_pore, num_p_atoms, 
                                                       guest_ptrs[index], p_demarc->end,
                                                       p_need_grad, p_grad, not_intra, 
                                                       dum_vdw_list, zero,
                                                       pos_separations, pos_vectors, FALSE);

//              printf("DB>> Non-bond energy after interacting molecule %d with the pore %10.6f\n",
//                                                               imol, interaction_energy[imol].non_bonded );

                  }
                else if (strcmp(pot_info.combination, ARITHMETIC) == 0 )
                  {
                         interaction_energy[imol].non_bonded +=
                                 non_bond_r_eps_arithmetic(p_pore, num_p_atoms,
                                                           guest_ptrs[index], p_demarc->end,
                                                           p_need_grad, p_grad, 
                                                           pos_separations, pos_vectors, FALSE);
       
//     if (DEBUG) printf("DB>> Non-bond energy after interacting molecule %d with the pore %10.6f\n",
//                                                                 imol, interaction_energy[imol].non_bonded );
                  }
                else
                  {
                     printf ("\nERROR: Attempt to use unknown combination ");
                     printf ("rules for this potential type\n when calculating vdw energy.\n");
                     printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                   pot_info.type, pot_info.combination);
                     printf ("Program terminating, goodbye.\n");
                     exit(EXIT_FAILURE);
                  }
            }
         else if ( strcmp(pot_info.type, A_B) == 0 )
           {    
              if (strcmp(pot_info.combination, GEOMETRIC) == 0 )
                {
                   interaction_energy[imol].non_bonded +=  
                          non_bond_a_b_geometric(p_pore, num_p_atoms, 
                                                 guest_ptrs[index], p_demarc->end,
                                                 p_need_grad, p_grad, 
                                                 pos_separations, pos_vectors, FALSE);

//     if (DEBUG) printf("DB>> Energy after interacting with the pore %10.6f\n",interaction_energy[imol].non_bonded );

                }
              else if (strcmp(pot_info.combination, NONE) == 0 )
                {
                   interaction_energy[imol].non_bonded +=  
                          non_bond_a_b_none(p_pore, num_p_atoms, 
                                            guest_ptrs[index], p_demarc->end,
                                            p_need_grad, p_grad, 
                                            pos_separations, pos_vectors, FALSE);

//     if (DEBUG) printf("DB>> Energy after interacting with the pore %10.6f\n",interaction_energy[imol].non_bonded );

                }
              else
                {
                   printf ("\nERROR: Attempt to use unknown combination ");
                   printf ("rules for this potential type\n when calculating vdw energy.\n");
                   printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                   pot_info.type, pot_info.combination);
                   printf ("Program terminating, goodbye.\n");
                   exit(EXIT_FAILURE);
                }
           }
         else
           {
              printf ("\nERROR: Attempt to use unknown potential ");
              printf ("type when calculating vdw energy.\n");
              printf ("Potential type is currently : %s\n", pot_info.type);
              printf ("Program terminating, goodbye.\n");
              exit(EXIT_FAILURE);
           }
          p_demarc++;
          free(pos_separations);
          free(pos_vectors);
         }
      }
/**** guest - guest interactions, do entire list on first pass ******/

//  printf("Non-bond contributions after host...guest interaction = %10.6f\n",
//                                                   interaction_energy[0].non_bonded);

//  printf("Starting on molmol list for %d molecule molecule interactions\n", num_molmol_list);

  p_this_list=p_molmol_list;
  for (ilist=0; ilist < num_molmol_list; ilist++)
    {
       imol=p_this_list->imol;    
       jmol=p_this_list->jmol;    
       p_i_demarc=p_guest_demarc+imol;
       p_j_demarc=p_guest_demarc+jmol;    

       ind=p_this_list->ind;
       jnd=p_this_list->jnd;

       if (ind==jnd) self_term=TRUE; else  self_term=FALSE;

/**** Count number of periodic images required and malloc accordingly    *******/
/**** note that this is the inclusion of mol..mol indices in list rather *******/
/**** than atom...atom.                                                  *******/

       list_size= count_vdw_pbc(guest_ptrs[jnd], p_j_demarc->end,  
                                guest_ptrs[ind], p_i_demarc->end,
                                self_term);

       if ( !self_term )
         {
                  
//             printf("\n Taking guest %d with guest %d for inter-molecular energy periodic image list_size = %d\n\n", 
//                                                   ind, jnd, list_size);

           if(list_size < 1)      pairs_to_malloc=10;
                              else pairs_to_malloc=list_size;

           pos_separations=(double*)malloc(pairs_to_malloc*sizeof(double));
           pos_vectors=(double*)malloc(3*pairs_to_malloc*sizeof(double));

           if ( strcmp(pot_info.type, R_EPS) == 0 )
             {
               if (strcmp(pot_info.combination, SIXTH_POWER) == 0 )
                 {
//                  printf("using r_eps_sixth non_bond\n");
                    gg_inter = non_bond_r_eps_sixth(guest_ptrs[jnd], p_j_demarc->end, 
                                                    guest_ptrs[ind], p_i_demarc->end,
                                                    p_need_grad, p_grad, not_intra, 
                                                    dum_vdw_list, zero,
                                                    pos_separations, pos_vectors, 
                                                    self_term);
                 }
             else if (strcmp(pot_info.combination, ARITHMETIC) == 0 )
               {
                    gg_inter =
                          non_bond_r_eps_arithmetic(guest_ptrs[jnd], p_j_demarc->end, 
                                                    guest_ptrs[ind], p_i_demarc->end,
                                                    p_need_grad, p_grad, 
                                                    pos_separations, pos_vectors,
                                                    self_term);

               }
             else
               {
                  printf ("\nERROR: Attempt to use unknown combination ");
                  printf ("rules for this potential type\n when calculating guest-guest vdw energy.\n");
                  printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                       pot_info.type, pot_info.combination);
                  printf ("Program terminating, goodbye.\n");
                  exit(EXIT_FAILURE);
               }
             }
          else if ( strcmp(pot_info.type, A_B) == 0 )
             {    
               if (strcmp(pot_info.combination, GEOMETRIC) == 0 )
                 {
                   gg_inter =  
                          non_bond_a_b_geometric(guest_ptrs[jnd], p_j_demarc->end, 
                                                 guest_ptrs[ind], p_i_demarc->end,
                                                 p_need_grad, p_grad, 
                                                 pos_separations, pos_vectors, 
                                                 self_term);

               }
             else if (strcmp(pot_info.combination, NONE) == 0 )
               {
                   gg_inter =
                          non_bond_a_b_none(guest_ptrs[jnd], p_j_demarc->end, 
                                            guest_ptrs[ind], p_i_demarc->end,
                                            p_need_grad, p_grad, 
                                            pos_separations, pos_vectors, self_term);

               }
             else
               {
                  printf ("\nERROR: Attempt to use unknown combination ");
                  printf ("rules for this potential type\n when calculating vdw energy for guest-guest interaction.\n");
                  printf ("Potential type is currently: %s, quoted combining rules: %s\n",
                                                       pot_info.type, pot_info.combination);
                  printf ("Program terminating, goodbye.\n");
                  exit(EXIT_FAILURE);
               }
           }

          free(pos_separations);
          free(pos_vectors);
/**** Divide the molecule-molecule interaction in half to avoid double counting when summing total ***/
          interaction_energy[imol].non_bonded += gg_inter/2.0; 
          interaction_energy[jmol].non_bonded += gg_inter/2.0; 

//          printf("Adding %10.6f interaction energy from this pass\n", gg_inter);
       }
      p_this_list++;
    }
//  printf("Non-bond contributions after all interactions = %10.6f\n",
//                                                   interaction_energy[imol].non_bonded);
  }
/*=============================================================================*/
/******** end non_bonded *******************************************************/
/*=============================================================================*/
//printf("DEBUG>> Stopping for debugging\n");
//exit(0);
/*=============================================================================*/
/******** calculate coulomb energy *********************************************/
/*=============================================================================*/

if (charges) 
  {
    if (pbc)
      {
/*******************************************************************************/
/***** Template-> framework terms **********************************************/
/*******************************************************************************/

         if (!is_empty_pore)
           {
             zero_first=TRUE;
             self_term= FALSE;
             using_host= TRUE;

             p_demarc=p_guest_demarc;    
             for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
               {
                  if (symm_set) index=imol*(num_symm_ops+2);
                                                 else index=imol;

                 ewald_sum( guest_ptrs[index], p_demarc->end, p_pore, num_p_atoms, 
                            p_kvecs, p_kvec2, p_gvec2, num_kvecs, 
                            p_cos_sum, p_sin_sum, zero_first, 
                            self_term, using_host); 

                 p_demarc++;
              }
           }

/*******************************************************************************/
/***** Basic template self interaction energy **********************************/
/***** using molmol list, added March 2013, Dave Willock ***********************/
/*******************************************************************************/

         zero_first=FALSE;
         if (is_empty_pore) zero_first=TRUE;
         using_host= FALSE;

         p_this_list=p_molmol_list;
         for (ilist=0; ilist < num_molmol_list; ilist++)
           {
             imol=p_this_list->imol;    
             jmol=p_this_list->jmol;    
             p_i_demarc=p_guest_demarc+imol;
             p_j_demarc=p_guest_demarc+jmol;    

             ind=p_this_list->ind;
             jnd=p_this_list->jnd;

             self_term= FALSE;
             if (ind == jnd) self_term= TRUE;

             ewald_sum( guest_ptrs[ind], p_i_demarc->end, guest_ptrs[jnd], p_j_demarc->end, 
                        p_kvecs, p_kvec2, p_gvec2, num_kvecs, p_cos_sum, p_sin_sum, 
                        zero_first, self_term, using_host);

             p_this_list++;
           }

/*******************************************************************************/
/***** Remove self interaction term with own atoms in the same unit cell *******/
/*******************************************************************************/

         zero_first=FALSE;
         self_term= TRUE;

         p_demarc=p_guest_demarc;    
         for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
           {
              if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

              for (isymm=0; isymm <= num_symm_ops; isymm++)
                {
                    isol_coul( guest_ptrs[index], p_demarc->end,guest_ptrs[index], p_demarc->end, zero_first, self_term);
                    index++;
                }
              p_demarc++;
           }

      }
    else
      {
         zero_first=TRUE;
         self_term= FALSE;

         p_demarc=p_guest_demarc;
         for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
           {
              if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

              for (isymm=0; isymm <= num_symm_ops; isymm++)
                {
                   isol_coul( guest_ptrs[index], p_demarc->end, p_pore, num_p_atoms, 
                              zero_first, self_term); 
                   index++;
                }
              p_demarc++;
           }
        }
/*******************************************************************************/
/***** Sum potentials **********************************************************/
/*******************************************************************************/
       
         interaction_energy[imol].charges = 0.0;

         p_demarc=p_guest_demarc;    
         for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
           {
              if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

              for (isymm=0; isymm <= num_symm_ops; isymm++)
                {
                   p_atom=guest_ptrs[index];
                   for ( iatom=0; iatom <= p_demarc->end; iatom++)
                     {
                       p_atom->electrostatic_pot = coul_prefactor * (p_atom->electrostatic_pot);

                       interaction_energy[imol].charges += (p_atom->part_chge)* (p_atom->electrostatic_pot);
                       p_atom++;
                     }
                  index++;
               }
             p_demarc++;
           } 
//if (DEBUG) 
//   { 
//      printf(" pref: %10.6f\n",coul_prefactor);
//      printf("DB>>Simple summation for Coulomb interaction energy gives: %10.6f\n", 
//                                                                 interaction_energy[imol].charges);
//   }

  }


/**********************************************************************/
/**** Calculate restraint energy DJW July 1998 ************************/
/**** Only makes sense if non_bond energy is being calculated *********/
/**********************************************************************/

if (non_bonded)
  {
     if (line_restraints)
       {
           printf("ERROR: restraints energy not implemented in multi-molecule version\n");
           exit(0);
          if (line_hold)
            {
               for (irestr=0; irestr <= num_line_restraints; irestr++)
                 {
          /*        find_line_atoms(irestr, p_templ, num_t_atoms, &iA, &iB); */

          /*        p_A = p_templ+iA;  */
          /*        p_B = p_templ+iB;  */
 
                    line_now[0] = p_B->x-p_A->x;
                    line_now[1] = p_B->y-p_A->y;
                    line_now[2] = p_B->z-p_A->z;

                    unit_vector( &line_now[0], &size);

                    dot = line_atom[irestr].ref_x*line_now[0]
                         +line_atom[irestr].ref_y*line_now[1]
                         +line_atom[irestr].ref_z*line_now[2];

                    if (dot >  1.0 && dot <  1.001) dot = 1.0;
                    if (dot < -1.0 && dot > -1.001) dot = -1.0;

                    if (dot >  1.0 || dot < -1.0) 
                      {
                         printf("ERROR: dot is out of range in restraints routine\n");
                         exit(0);
                      }

                    theta= acos(dot);

/***** Shift obtuse angles back to the equivalent angular displacement from the origin **/

                    if (theta > pi_b2) theta= theta-pi;
                    if (theta < -pi_b2) theta= theta+pi;

                    interaction_energy[imol].restraint += 
                                   line_atom[irestr].k*(theta*RAD_TO_DEG)*(theta*RAD_TO_DEG);
  
                 }
            }
         if (have_tethers)
            {
           printf("ERROR: tether energy not implemented in multi-molecule version\n");
           exit(0);

               if (is_empty_pore)
                 {
                   printf("ERROR: Empty pore cannot be used to define tethers\n");
                   exit(0);
                 }

         /*    p_patom= p_pore+tether.index_a;  */
         /*    p_tatom= p_templ+tether.index_b; */

               dx = p_tatom->x - p_patom->x;
               dy = p_tatom->y - p_patom->y;
               dz = p_tatom->z - p_patom->z;
              
               if (pbc) min_image(&dx, &dy, &dz);
               r= dx*dx + dy*dy + dz*dz; 
               r= sqrt(r);

               interaction_energy[imol].restraint += tether.k * (r - tether.r0 ) * (r - tether.r0 );
            }
       }

/************* Try and get hbond energy *************************************************/
   
      if ( strcmp(pot_info.hbond, H_BOND_AMBER) == 0 ) 
        {

           printf("ERROR: hbond energy not implemented in multi-molecule version\n");
           exit(0);

/*         interaction_energy.hbond=                                                            */
/*                 hbond_amber(p_pore, num_p_atoms, p_templ, num_t_atoms, p_need_grad, p_grad); */
        }
  }


//printf("Returning from calculate_energy\n");

return;
}

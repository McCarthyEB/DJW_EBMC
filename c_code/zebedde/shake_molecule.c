#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "data.h"
#include "global_values.h"
#include "constants.h"

int cost_function(atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc,
                  atom *p_pore, double *p_box_limits, 
                  interaction_indices *p_molmol_list,
                  int num_molmol_list,
                  int which_test, double *p_cc_tot, 
                  int flags_set, int *p_flag_list,
                  double *p_kvecs, double *p_kvec2,
                  double *p_gvec2, int num_kvecs,
                  double *p_cos_sum, double *p_sin_sum,
                  int *p_need_grad, double *p_grad,
                  int have_comb_rules, int is_empty_pore,
                  angle_interact_list *p_angles_list_ptrs[],
                  int *p_num_angles_list, 
                  torsion_interact_list *p_torsions_list_ptrs[],
                  int *p_num_torsions_list, 
                  vdw_interact_list *p_vdw_list_ptrs[],
                  int *p_num_vdws_list,
                  int need_intra);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm );

double real_random(int done);

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *guest_ptrs[], int num_guests,
                      list_partition *p_guest_demarc,
                      interaction_indices *p_molmol_list,
                      int num_molmol_list,
                      double *p_kvecs, double *p_kvec2,
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad,
                      int have_comb_rules, int is_empty_pore);

void unit_vector(double *p_vector, double *p_size);

/******* random rigid body move of template aimed at  ***********/
/******* maximizing the total of the template -> pore ***********/
/******* close contact list                           ***********/

int  shake_molecule(atom *guest_ptrs[], int num_guests,
                    list_partition *p_guest_demarc,
                    atom *p_pore, double *p_box_limits,
                    interaction_indices *p_molmol_list, 
                    int num_molmol_list,
                    symm_ops *p_symm, int which_test, 
                    int num_attempts, double max_move, 
                    double *p_kvecs, double *p_kvec2,
                    double *p_gvec2, int num_kvecs,
                    double *p_cos_sum, double *p_sin_sum,
                    int *p_need_grad, double *p_grad,
                    int have_comb_rules, 
                    int is_empty_pore,
                    int which_mol)
{
#include "header.h"

int iatom, isign, reject_move, accepted, num_moves;
int imol, index, index_image, isymm, iii;
int dummy_int, allow, done_move;
double cc_tot=0.0, delta[3], dist, size, ediff; 
double before_temp[MAX_MOLS];

list_partition *p_demarc;

double last_cc_tot=0.0, last_energy=0.0, total_energy=0.0;
double rand;

angle_interact_list *p_angles_dummy_list[1];
torsion_interact_list *p_torsions_dummy_list[1];
vdw_interact_list *p_vdw_dummy_list[1];

/**************************************************************************/
/** Note that shake operation will not affect the intra-molecular *********/
/** energy so no need to include that in last_energy etc.         *********/
/**************************************************************************/

if (DEBUG) 
  {
    fprintf(output_fp,"DEBUG: Entered shake_molecule, which_test=%d\n", which_test);
    fprintf(output_fp,"DEBUG: Shaking molecule %d of %d guests.\n", which_mol, num_guests);
  } 

calculate_energy(p_pore, num_pore_atoms,
                 guest_ptrs, num_guests, p_guest_demarc,
                 p_molmol_list, num_molmol_list,
                 p_kvecs, p_kvec2, p_gvec2,
                 num_kvecs, p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                 have_comb_rules, is_empty_pore);

/**** No longer need intra energy for a rock as each molecule is in a separate     ***/
/**** list. Dave Willock March 2013                                                ***/

/**************************************************************************/
/**** Initialise the indexing to the required molecule ********************/
/**************************************************************************/
  p_demarc=p_guest_demarc+which_mol;
  if (symm_set) index=which_mol*(num_symm_ops+2);
                                         else index=which_mol;

  if (which_test == STERIC_COST)
    {
       last_cc_tot=0;
       for (iatom= 0; iatom <= p_demarc->end; iatom++)
                                    last_cc_tot+= closest_contact[iatom];

    }
  else if (which_test == NON_BOND_COST)
    {
       last_energy=0.0;
       for (imol=0; imol < num_guests; imol++)
                          last_energy+= interaction_energy[imol].non_bonded;

/*********************************************/
/**** For DEBUG remember initial energies ****/
           for (imol=0; imol < num_guests; imol++)
                          before_temp[imol]= interaction_energy[imol].non_bonded;
/**** For DEBUG remember initial energies ****/
/*********************************************/

    }
  else if (which_test == ENERGY_COST)
    {
       last_energy=0.0;
       for (imol=0; imol < num_guests; imol++)
                          last_energy+= interaction_energy[imol].non_bonded
                                      + interaction_energy[imol].charges;
    }

/************************************************************************/
/*** Attempt num_attempts trial moves along random vector. **************/
/************************************************************************/

  delta[0]= real_random(1);
  isign=  real_random(1) < 0.5;
  if (isign) delta[0]= -1.0*delta[0];

  delta[1]= real_random(1);
  isign=  real_random(1) < 0.5;
  if (isign) delta[1]= -1.0*delta[1];

  delta[2]= real_random(1);
  isign=  real_random(1) < 0.5;
  if (isign) delta[2]= -1.0*delta[2];

  unit_vector(&delta[0], &size);
  dist = real_random(1)*max_move;

  delta[0]= dist*delta[0];
  delta[1]= dist*delta[1];
  delta[2]= dist*delta[2];

  accepted= TRUE;

/************************************************************************/
/***** push in the trial direction till the cost function says stop *****/
/***** or we run out of patience                                    *****/
/************************************************************************/

//  printf("\nDB >> Trying a shake move for molecule %d\n", which_mol);

  num_moves=0;
  done_move=FALSE;
  while (accepted && num_moves < num_attempts)
    {
       num_moves++;

       move_molecule(guest_ptrs[which_mol], p_demarc->end, &delta[0], -1);

/****** Recalculate symmetry related atoms  ***************/


        if (symm_set)
          {
             for (isymm=0; isymm <= num_symm_ops; isymm++)
               {
                  index_image= index+isymm+1;
                  printf("        ....doing image %d with %d atoms to index %d\n",
                                                    isymm, p_demarc->num, index_image);

                  do_symmetry( guest_ptrs[index], p_demarc->end,
                               guest_ptrs[index_image], &symm[isymm] );
               }
          }


       reject_move = cost_function(guest_ptrs, num_guests, p_guest_demarc,
                                   p_pore, p_box_limits,
                                   p_molmol_list, num_molmol_list,
                                   which_test, &cc_tot,
                                   FALSE, &dummy_int,
                                   p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                   p_cos_sum, p_sin_sum,
                                   p_need_grad, p_grad,
                                   have_comb_rules,
                                   is_empty_pore, 
                                   p_angles_dummy_list, &dummy_int,
                                   p_torsions_dummy_list, &dummy_int, 
                                   p_vdw_dummy_list,&dummy_int, 
                                   FALSE);

//        printf("Back in shake_molecule\n");
/**** ignore bump check as other guests may have set this *****/
        reject_move = FALSE;

       if (which_test == STERIC_COST) 
         {
           reject_move= reject_move || cc_tot <= last_cc_tot;

if (DEBUG)
  {
     printf("DB >> shake_mol last_cc_tot = %10.6f cc_tot= %10.6f ave_cc %10.6f\n",
                                               last_cc_tot,cc_tot, last_cc_tot/p_demarc->num);
  }

         }
       else if (which_test == NON_BOND_COST) 
         {

/***** do real MC ********************/
/***** Added 1/8/06 Dave Willock *****/


           total_energy=0.0;
           for (imol=0; imol < num_guests; imol++)
                          total_energy+= interaction_energy[imol].non_bonded;

           ediff = total_energy-last_energy;
           rand  = real_random(1);

//           printf("DB>> In shake last_energy=%10.6f, total_energy=%10.6f\n",
//                                     last_energy, total_energy);

           allow = ediff <= 0.0 || rand < exp(-ediff/(BOLTZ*temperature));

           if (DEBUG)
             {
               for (imol=0; imol < num_guests; imol++)
                  fprintf(output_fp,"DEBUG: Energy for mol %d now: %10.6f before move: %10.6f\n",
                                                       imol,
                                                       interaction_energy[imol].non_bonded,
                                                       before_temp[imol]);

                  fprintf(output_fp,"DEBUG: total: %10.6f last: %10.6f\n", 
                                                       total_energy, last_energy);

                  fprintf(output_fp,"DEBUG: ediff = %10.6f, temp: %10.6f, BOLTZ*temp.: %10.6f\n",
                                                                ediff, temperature, BOLTZ*temperature);

                  fprintf(output_fp,"DEBUG: boltzmann factor: %10.6f random number %10.6f\n", 
                                                                exp(-ediff/(BOLTZ*temperature)), rand);

                  if (allow) fprintf(output_fp,"allow this shake\n");
                                   else fprintf(output_fp,"DO NOT allow this shake\n");
             }

           reject_move= reject_move || !allow;
                
         }
       else if (which_test == ENERGY_COST) 
         {
           total_energy= 0.0;
           for (imol=0; imol < num_guests; imol++)
                          total_energy+= interaction_energy[imol].non_bonded
                                       + interaction_energy[imol].charges;

           reject_move= reject_move || total_energy >= last_energy;
         }


if ( which_test == STERIC_COST && last_cc_tot/p_demarc->num > MAX_AVE_CC)
    {
    fprintf(output_fp,"ZEBEDDE ERROR: Close Contact Enormous in Shake last_cc_tot=%10.6f\n", last_cc_tot);
    fflush(output_fp);
    }

            if (!reject_move)
              {
                 if (which_test == STERIC_COST) 
                   {
                     last_cc_tot= cc_tot;
                   }
                 else if (which_test == NON_BOND_COST) 
                   {
                     last_energy= total_energy;
                   }
                 if (which_test == ENERGY_COST) 
                   {
                     last_energy= total_energy;
                   }
                 done_move=TRUE;
              }
            else
              {
                 delta[0] = -1.0 * delta[0];
                 delta[1] = -1.0 * delta[1];
                 delta[2] = -1.0 * delta[2];
                 move_molecule(guest_ptrs[which_mol], p_demarc->end, &delta[0], -1);
                 accepted= FALSE;

/****** Recalculate symmetry related atoms  ***************/

                 if (symm_set)
                   {
                     for (isymm=0; isymm <= num_symm_ops; isymm++)
                       {
                          index_image= index+isymm+1;
                          printf("        ....doing image %d with %d atoms to index %d\n",
                                                            isymm, p_demarc->num, index_image);

                          do_symmetry( guest_ptrs[index], p_demarc->end,
                                       guest_ptrs[index_image], &symm[isymm] );
                       }
                   }

/*** If this move is straight after a success this image will be accepted so redo energy **/
/*** Added Aug 09 Dave Willock ************************************************************/

                if (done_move)
                   {
                      calculate_energy(p_pore, num_pore_atoms,
                                       guest_ptrs, num_guests, p_guest_demarc,
                                       p_molmol_list, num_molmol_list,
                                       p_kvecs, p_kvec2, p_gvec2,
                                       num_kvecs, p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                                       have_comb_rules, is_empty_pore);
                   }
              }
          }

return done_move;
}

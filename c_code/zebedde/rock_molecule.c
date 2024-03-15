/***************************************************************/
/******* Rigid body rotation aimed at maximising the ***********/
/******* template -> pore close contact total distance *********/
/***************************************************************/

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

void do_symmetry( atom *p_molecule, int num_atoms, atom *p_image,
                  symm_ops *p_symm );

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_atoms, int which_mol);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

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

double real_random(int done);

void print_dist_mat( atom *p_molecule, int num_atoms, int use_pbc, FILE *fp);

void unit_vector(double *p_vector, double *p_size);

int rock_molecule(atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc,
                  atom *p_pore, double *p_box_limits,
                  interaction_indices *p_molmol_list, int num_molmol_list,
                  symm_ops *p_symm, int which_test, int num_attempts,
                  double max_rot, double *p_kvecs, double *p_kvec2,
                  double *p_gvec2, int num_kvecs,
                  double *p_cos_sum, double *p_sin_sum,
                  int *p_need_grad, double *p_grad, int have_comb_rules,
                  int is_empty_pore,
                  int which_mol)
{
#include "header.h"
   int  iatom, isign, reject_rot, imol, index, index_image; 
   int isymm;
   int num_symm_atoms,dummy_int, allow;
   int done_rot, accepted=FALSE, num_rots=0;

   list_partition *p_demarc;

   double origin[3], axis[3], cc_tot, theta;
   double last_cc_tot=0.0;
   double last_energy=0.0;
   double total_energy=0.0;
   double total_mass, size;
   double ediff, rand;

  angle_interact_list *p_angles_dummy_list[1];
  torsion_interact_list *p_torsions_dummy_list[1];
  vdw_interact_list *p_vdw_dummy_list[1];

   done_rot=FALSE;

/***************************************************************************/
/***** Start of executable lines *******************************************/
/***************************************************************************/

// printf("DEBUG>> Entered rock_molecule, which_test=%d\n", which_test);

calculate_energy(p_pore, num_pore_atoms, 
                 guest_ptrs, num_guests, p_guest_demarc,
                 p_molmol_list, num_molmol_list,
                 p_kvecs, p_kvec2, p_gvec2,
                 num_kvecs, p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                 have_comb_rules, is_empty_pore);

/**** No longer need intra energy for a rock as it each molecule is in a separate  ***/
/**** list. Dave Willock March 2013                                                ***/

//if (DEBUG)
//  {
//    printf("DB >> Trying a rock for molecule %d\n", which_mol);
//  }

/**************************************************************************/
/**** Initialise the required test ****************************************/
/**************************************************************************/
  p_demarc=p_guest_demarc+which_mol;
  if (symm_set) index=which_mol*(num_symm_ops+2);
                                         else index=which_mol;

  if (which_test == STERIC_COST)
    {
      last_cc_tot=0.0;

      for (iatom= 0; iatom <= p_demarc->end; iatom++)
                                    last_cc_tot+= closest_contact[iatom];
    }
  else if (which_test == NON_BOND_COST)
    {
       last_energy=0.0;
       for (imol=0; imol < num_guests; imol++)
                          last_energy+= interaction_energy[imol].non_bonded;
    }
  else if (which_test == ENERGY_COST)
    {
       last_energy=0.0;
       for (imol=0; imol < num_guests; imol++)
                          last_energy+= interaction_energy[imol].non_bonded
                                      + interaction_energy[imol].charges;
    }

/***************************************************************************/
/******** find new centre of mass ******************************************/
/***************************************************************************/

  centre_of_mass(&origin[0], &total_mass,
                            guest_ptrs[index], p_demarc->end, -1 );

  if (DEBUG) print_dist_mat(guest_ptrs[index], p_demarc->end,
                                    TRUE, stdout);

/********************************************************/
/*** Decide on angle to rotate by ***********************/
/********************************************************/

 theta = real_random(1)*max_rot;
 isign =  real_random(1) < 0.5;
 if (isign) theta= -1.0*theta;

/********************************************************/
/*** Random vector in space by assigning            *****/
/*** random components.                             *****/
/*** Must go to rotate with axis as a unit vector   *****/
/********************************************************/
         axis[0]= 1.0*real_random(1);
         axis[1]= 1.0*real_random(1);
         axis[2]= 1.0*real_random(1);

         unit_vector(&axis[0], &size);

//       if (DEBUG) 
//         {
//           printf("Selected axis %5.2f %5.2f %5.2f by angle : %10.6f degrees\n",
//                                   axis[0],axis[1],axis[2], theta*RAD_TO_DEG);
//
//           printf("DB >> Still trying a rock %d %d \n",accepted,num_rots);
//         }
 

         num_rots=0;
         accepted= TRUE;
         while (accepted && num_rots < num_attempts)
           { 
             num_rots++;

//           if (DEBUG) 
//             {
//               printf("Before this rotation step %d neighbour distances are:\n", num_rots);
//               print_dist_mat(guest_ptrs[index], p_demarc->end,
//                                                                 TRUE, stdout);
//             }

             rotate(guest_ptrs[index], &axis[0],
                          &origin[0], theta, p_demarc->end, -1);

//           if (DEBUG) 
//             {
//               printf("After this rotation step %d neighbour distances are:\n", num_rots);
//               print_dist_mat(guest_ptrs[index], p_demarc->end,
//                                                                 TRUE, stdout);
//             }

/****** Recalculate symmetry related atoms  ***************/

             if (symm_set)
               {
                  for (isymm=0; isymm <= num_symm_ops; isymm++)
                    {
                       index_image= index+isymm+1;
//                     printf("        ....doing image %d with %d atoms to index %d\n",
//                                                       isymm, p_demarc->num, index_image);

                       do_symmetry( guest_ptrs[index], p_demarc->end,
                                    guest_ptrs[index_image], &symm[isymm] );
                    }
               }

            reject_rot  = cost_function(guest_ptrs, num_guests, p_guest_demarc, 
                                        p_pore, p_box_limits, p_molmol_list, num_molmol_list,
                                        which_test, &cc_tot, FALSE, &dummy_int,
                                        p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                        p_cos_sum, p_sin_sum,
                                        p_need_grad, p_grad, have_comb_rules,
                                        is_empty_pore, 
                                        p_angles_dummy_list, &dummy_int,
                                        p_torsions_dummy_list, &dummy_int, 
                                        p_vdw_dummy_list,&dummy_int, 
                                        FALSE);

/**** ignore bump check as other guests may have set this *****/
        reject_rot = FALSE;


            if (which_test == STERIC_COST)  
              { 
                reject_rot= reject_rot  || cc_tot <= last_cc_tot;
              }
            else if (which_test == NON_BOND_COST)
              {

/***** do real MC ********************************************************************/
/***** Added 1/8/06 Dave Willock *****************************************************/
/***** Note that rock will not affect intra molecular energy so no need for that *****/

//              if (DEBUG) printf("DB rot: non bond Energy for molecule %d now: %10.6f last: %10.6f", 
//                                         which_mol, interaction_energy[which_mol].non_bonded, last_energy);

                total_energy=0.0;
                for (imol=0; imol < num_guests; imol++)
                                   total_energy+= interaction_energy[imol].non_bonded;

                ediff = total_energy-last_energy;
                rand  = real_random(1);

                allow = ediff <= 0.0 || rand < exp(-ediff/(BOLTZ*temperature));

//              if (DEBUG)
//                {
//                printf("DEBUG: ediff = %10.6f, temp: %10.6f, BOLTZ*temp.: %10.6f arg:%10.6f\n",
//                                ediff, temperature, BOLTZ*temperature, -ediff/(BOLTZ*temperature));
//
//                  if (allow) 
//                    {
//                      printf("allow this rock\n");
//                      if (ediff > 0.0)
//                        {
//                           printf("DEBUG: boltzmann factor: %10.6f random number %10.6f\n",
//                                                    exp(-ediff/(BOLTZ*temperature)), rand);
//                        }
//                    }
//                else printf("DO NOT allow this rock\n");
//                }

                reject_rot= reject_rot || !allow;
                
              }
            else if (which_test == ENERGY_COST)
              {
                total_energy=0.0;
                total_energy=   interaction_energy[which_mol].non_bonded 
                              + interaction_energy[which_mol].charges;
                reject_rot= reject_rot || total_energy >= last_energy;
              }


     if (which_test == STERIC_COST && last_cc_tot/p_demarc->num > MAX_AVE_CC) 
       {
          printf("Error: ave cc_tot greater than MAX_AVE_CC in rock_molecule\n");
          printf("DB >> Rock_mol last_cc_tot = %10.6f cc_tot= %10.6f\n",last_cc_tot,cc_tot);
          
          exit(-1);
       }

             if (!reject_rot)
               {
//               if (DEBUG) printf("Accepting rotation step\n");
                 if (which_test == STERIC_COST) last_cc_tot= cc_tot;
                 if (which_test == NON_BOND_COST) last_energy= interaction_energy[which_mol].non_bonded;
                 if (which_test == ENERGY_COST) last_energy= total_energy;
                 done_rot= TRUE;
               }
             else
               {
//                if (DEBUG) printf("Rejecting rotation step, winding back %d atoms\n",
//                                                     p_demarc->num);

/*** BUG fix March 07 Dave Willock ****/
/*** Reverse axis, for rejecting   ****/
/*** Used to only do one component ****/

                  axis[0] = -1.0 * axis[0];
                  axis[1] = -1.0 * axis[1];
                  axis[2] = -1.0 * axis[2];

                  rotate(guest_ptrs[index], &axis[0],
                              &origin[0], theta, p_demarc->end, -1);
                  accepted= FALSE;

/****** Recalculate symmetry related atoms  ***************/

                 if (symm_set)
                   {
                     for (isymm=0; isymm <= num_symm_ops; isymm++)
                       {
                          index_image= index+isymm+1;
//                        printf("        ....doing image %d with %d atoms to index %d\n",
//                                                          isymm, p_demarc->num, index_image);

                          do_symmetry( guest_ptrs[index], p_demarc->end,
                                       guest_ptrs[index_image], &symm[isymm] );
                       }
                   }

/**** If this is a step after a success we will allow the new configuration so need energy update ***/
/**** Dave Willock Aug 09 ***************************************************************************/
                 if (done_rot)
                   {
                      calculate_energy(p_pore, num_pore_atoms, 
                                       guest_ptrs, num_guests, p_guest_demarc,
                                       p_molmol_list, num_molmol_list,
                                       p_kvecs, p_kvec2, p_gvec2,
                                       num_kvecs, p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                                       have_comb_rules, is_empty_pore);
//intf("DEBUG>> Back in rock_molecule, 2\n");
                   }
               }
           }
return done_rot;
}


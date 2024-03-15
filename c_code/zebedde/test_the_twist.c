/*******************************************************/
/***** routine to rotate section of molecule   *********/
/***** assumes first atom in the list is fixed *********/
/***** Dave Willock April 21st 1995            *********/
/***** Adapted for intra-molecular energy terms*********/
/***** Oct. 2006 Dave Willock                  *********/
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

/* required functions ------------------------------------*/

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

/*--------------------------------------------------------*/

int test_the_twist(atom *p_pore, atom *guest_ptrs[], int num_guests,
                   list_partition *p_guest_demarc,
                   double *p_box_limits,
                   interaction_indices *p_molmol_list,
                   int num_molmol_list,
                   int which_test, 
                   double *p_kvecs, double *p_kvec2,
                   double *p_gvec2, int num_kvecs,
                   double *p_cos_sum, double *p_sin_sum,
                   int *p_need_grad, double *p_grad, int have_comb_rules,
                   int is_empty_pore, 
                   angle_interact_list *p_angles_list_ptrs[],
                   int *p_num_angles_list, 
                   torsion_interact_list *p_torsions_list_ptrs[],
                   int *p_num_torsions_list, 
                   vdw_interact_list *p_vdw_list_ptrs[],
                   int *p_num_vdws_list,
                   int which_mol)
{
#include "header.h"

double origin[3], axis[3];
int *p_flag_list, *p_this_flag;
int atom1,atom2,flags_set, do_twist;
double  rand, last_cc_tot, cc_tot;
int num_rots,isign, iatom, which_bit, index, index_image;
int allow, reject_twist;
int icount, iii, imol;
int just_count, num_vdws;
atom *p_atom1, *p_atom2;

list_partition *p_demarc;

double last_energy, total_energy, ediff;

/****************************************************************/
/**** Start of Executable lines *********************************/
/****************************************************************/
 p_demarc=p_guest_demarc+which_mol;

/****** Initialise the test variables *********************************/
   last_cc_tot=0.0;
   last_energy=0.0;
   total_energy=0.0;

   if (which_test == STERIC_COST)
     {
        cc_tot= 100.0;
        last_cc_tot= 0.0;
        for (iatom= 0; iatom <= p_demarc->end; iatom++)
                                         last_cc_tot+= closest_contact[iatom];
     }
   else if (which_test == NON_BOND_COST)
     {
       last_energy=0.0;
       for (imol=0; imol< num_guests; imol++)
          {
                          last_energy+=   interaction_energy[imol].non_bonded 
                                        + intra_energy[imol].total;
//       printf("DEBUG>> Last energy for twist is %10.6f vdw: %10.6f intra: %10.6f\n", last_energy,
//                                                       interaction_energy[imol].non_bonded,
//                                                       intra_energy[imol].total);
          }
     }
   else if (which_test == ENERGY_COST)
     {
       last_energy=   interaction_energy[which_mol].non_bonded 
                    + interaction_energy[which_mol].charges;
     }

if ( which_test == STERIC_COST && last_cc_tot/p_demarc->num > MAX_AVE_CC) 
    {
    fprintf(output_fp,"ZEBEDDE ERROR: Close Contact Enormous before Rotation\n");
    fflush(output_fp);
    exit(-1);
    }

    reject_twist = FALSE;
    allow = TRUE;

        reject_twist = cost_function(guest_ptrs, num_guests, p_guest_demarc, 
                                     p_pore, p_box_limits, p_molmol_list, num_molmol_list,
                                     which_test, &cc_tot, FALSE, p_flag_list,
                                     p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                     p_cos_sum, p_sin_sum,
                                     p_need_grad, p_grad,
                                     have_comb_rules, 
                                     is_empty_pore,
                                     p_angles_list_ptrs, p_num_angles_list,
                                     p_torsions_list_ptrs, p_num_torsions_list, 
                                     p_vdw_list_ptrs, p_num_vdws_list,
                                     TRUE);


if ( which_test == STERIC_COST && last_cc_tot/p_demarc->num > MAX_AVE_CC)
    {
    fprintf(output_fp,"ZEBEDDE ERROR: Close Contact Enormous during Rotation\n");
    fflush(output_fp);
    exit(-1);
    }

        if (which_test == STERIC_COST) 
          {
            reject_twist= reject_twist || cc_tot <= last_cc_tot;
            allow= !reject_twist;
          }
        else if (which_test == NON_BOND_COST)
           {
              if (reject_twist)
                {
//                  printf("Twist rejected by cost_function.....\n");
                  total_energy= last_energy;
                  allow=FALSE;
                }
              else
                {
/*****************************************************/
/** twist will also affect intra-molecular energy ****/
/*****************************************************************************************/
/** Note that num_guests is the real number of guests rather than an upper array limit ***/
/*****************************************************************************************/
                  total_energy=0.0;
//                  printf("Summing %d contributions\n", num_guests);
                  for (imol=0; imol< num_guests; imol++)
                     {
                          total_energy+=   interaction_energy[imol].non_bonded 
                                         + intra_energy[imol].total;

//                          printf("Contribution from molecule %d interaction: %10.6f intra: %10.6f\n",
//                                       imol, interaction_energy[imol].non_bonded, intra_energy[imol].total);
                     }

//                  printf("DEBUG>>.... test_twist last energy: %10.6f\n", last_energy);
//                  printf("DEBUG>>.... test_twist      energy: %10.6f\n", total_energy);
                  ediff = total_energy - last_energy;

//                  printf("DEBUG>> ediff twist energy: %10.6f\n", ediff);
                 rand  = real_random(1);

                 allow = ediff <= 0.0 || rand < exp(-ediff/(BOLTZ*temperature));
               }
           }

        else if (which_test == ENERGY_COST)
          {
             total_energy= interaction_energy[which_mol].non_bonded + interaction_energy[which_mol].charges;
 
             rand  = real_random(1);

             allow = ediff <= 0.0 || rand < exp(-ediff/(BOLTZ*temperature));
          }

        do_twist=FALSE;
        if (!reject_twist && allow) do_twist=TRUE;

//                   if (do_twist) printf("This will be accepted...\n");
return do_twist;
}

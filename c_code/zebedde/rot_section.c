/*******************************************************/
/***** routine to rotate section of molecule   *********/
/***** assumes first atom in the list is fixed *********/
/***** Dave Willock April 21st 1995            *********/
/***** Adapted for intra-molecular energy and  *********/
/***** Real MC testing, October 2006           *********/
/***** Dave Willock                            *********/
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"
#include "header.h"

void rotate(atom *p_molecule, double *p_axis, double *p_origin,
            double theta, int num_template_atoms);

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);

void unit_vector(double *p_vector, double *p_size);

double real_random(int done);

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

void rot_section(atom *guest_ptrs[], int num_guests,
                 list_partition *p_guest_demarc,
                 atom *p_pore, double *p_box_limits,
                 interaction_indices *p_molmol_list, 
                 int num_molmol_list,
                 symm_ops *p_symm,
                 int which_test, links *p_latest_link,
                 int num_to_rot, double max_rot,
                 double *p_kvecs, double *p_kvec2,
                 double *p_gvec2, int num_kvecs,
                 double *p_cos_sum, double *p_sin_sum,
                 int *p_need_grad, double *p_grad,
                 int have_comb_rules, 
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
atom *p_start, *p_temp_atom;

double origin[3], axis[3], theta ,cc_tot, last_cc_tot;
double size;
double last_energy, total_energy;

list_partition *p_demarc;

int reject_rot, num_rots, isign, iatom, start_atom, template_atom;
int index, dummy_int, isymm, index_image;

p_demarc=p_guest_demarc+which_mol;

if (DEBUG)
   {
printf("DB >> Trying rotate last fragment joined\n");
   }

   start_atom    = p_latest_link->end;
   template_atom = p_latest_link->start;

   p_start    = guest_ptrs[index]+ start_atom;
   p_temp_atom= guest_ptrs[index]+ template_atom;

if (DEBUG)
   {
   printf("DB >> Fixing atom No %d label %s\n", start_atom, p_start->label);
    }

   origin[0]= p_start->x;
   origin[1]= p_start->y;
   origin[2]= p_start->z;

   join_atoms( p_start, p_temp_atom, &axis[0]);
   unit_vector(&axis[0], &size);
  
   theta = real_random(1)*max_rot;
   isign =  real_random(1) < 0.5;
   if (isign == 1) theta= -1.0*theta;

/**********************************************************************/
/****** Initialise the relavent tests *********************************/
/**********************************************************************/
   last_cc_tot= 0.0;
   last_energy= 0.0;
   total_energy= 0.0;

   if (which_test == STERIC_COST)
     {
       for (iatom= 0; iatom <= p_demarc->end; iatom++)
                                         last_cc_tot+= closest_contact[iatom];
     }
   else if (which_test == NON_BOND_COST)
     {
       last_energy= interaction_energy[which_mol].non_bonded;
     }
   else if (which_test == ENERGY_COST)
     {
       last_energy=  interaction_energy[which_mol].non_bonded 
                   + interaction_energy[which_mol].charges;
     }


if ( which_test == STERIC_COST && last_cc_tot/p_demarc->num > MAX_AVE_CC)
    {
    fprintf(output_fp,"ZEBEDDE ERROR: Close Contact Enormous before Last Section Rotation\n");
    fflush(output_fp);
    exit(-1);
    }

   reject_rot = FALSE;
   num_rots= 0;

   while (!reject_rot && num_rots < 20)
     {
        num_rots++;
        rotate(p_start, &axis[0], &origin[0], theta, num_to_rot);

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

        reject_rot =  cost_function(guest_ptrs, num_guests, p_guest_demarc,
                                    p_pore, p_box_limits,
                                    p_molmol_list, num_molmol_list,
                                    which_test, &cc_tot,
                                    FALSE, &dummy_int, 
                                    p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                                    p_cos_sum, p_sin_sum,
                                    p_need_grad, p_grad,
                                    have_comb_rules, is_empty_pore, 
                                    p_angles_list_ptrs, p_num_angles_list,
                                    p_torsions_list_ptrs, p_num_torsions_list, 
                                    p_vdw_list_ptrs, p_num_vdws_list,
                                    TRUE);

          
     if (DEBUG)
       {
   printf("DB >> rot_section last_cc_tot = %10.6f cc_tot= %10.6f, theta = %10.6f\n"
                                ,last_cc_tot,cc_tot, theta);
       }

        if (which_test == STERIC_COST) 
          {
            reject_rot= reject_rot || cc_tot <= last_cc_tot;
          }
        else if (which_test == NON_BOND_COST)
          {
            reject_rot= reject_rot || 
                    interaction_energy[which_mol].non_bonded  >= last_energy;
                  
            if (reject_rot)
                    reject_rot = real_random(1) > 
                                exp(-(interaction_energy[which_mol].non_bonded -last_energy) 
                                                               /(BOLTZ*temperature));
          }
        else if (which_test == ENERGY_COST)
          {
            total_energy=   interaction_energy[which_mol].non_bonded 
                          + interaction_energy[which_mol].charges;

            reject_rot= reject_rot || total_energy >= last_energy;
          }

        if (!reject_rot)
          {
             if (which_test == STERIC_COST) 
               {
                  last_cc_tot= cc_tot;
               }
             else if (which_test == NON_BOND_COST)
               {
                  last_energy= interaction_energy[which_mol].non_bonded;
               }
             else if (which_test == ENERGY_COST)
               {
                  last_energy= total_energy;
               }
          }

     }


/********** Last rotation must have been rejected so reverse it ********/

  theta = -theta;
  rotate(p_start, &axis[0], &origin[0], theta, num_to_rot);

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


  return;
}

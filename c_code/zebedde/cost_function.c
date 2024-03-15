/******************************************************************************/
/***************** Cost function used to decide if a template *****************/    
/***************** is acceptable                              *****************/    
/*****************                                            *****************/    
/***************** Started Wednesday May 24th 1995            *****************/    
/***************** Last updated 6th July 1998                 *****************/    
/*****************                                            *****************/    
/***************** Update October 2006                        *****************/    
/***************** Make Zebedde able to do Monte Carlo        *****************/
/***************** type assessment so that this routine       *****************/
/***************** can use energy and temperature to make its *****************/
/***************** decision.                                  *****************/
/***************** Update May 2013 if called by shake or rock *****************/    
/***************** there can be no change of conformation so  *****************/    
/***************** added flag "need_intra".                   *****************/    
/***************** Dave Willock                               *****************/    
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

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

void calculate_intra_energy( atom *guest_ptrs[], int num_guests,
                             list_partition *p_guest_demarc,
                             angle_interact_list *p_angles_list_ptrs[],
                             int *p_num_angles_list,
                             torsion_interact_list *p_torsions_list_ptrs[],
                             int *p_num_torsions_list,
                             vdw_interact_list *p_vdw_list_ptrs[],
                             int *p_num_vdws_list);

int bump_check(atom *guest_ptrs[], int num_guests,
               list_partition *p_guest_demarc,
               interaction_indices *p_molmol_list, int num_molmol_list);

int bump_check_with_flags(atom *guest_ptrs[], int num_guests,
                          list_partition *p_guest_demarc,
                          interaction_indices *p_molmol_list, int num_molmol_list,
                          int *p_flag_list);

int box_test( double *p_box_limits, atom *guest_ptrs[], int num_guests,
                  list_partition *p_guest_demarc);

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
                  int need_intra)
{
#include "header.h"
int iatom, out_of_your_box, inter_guest_bump=FALSE;
int imol, index;

list_partition *p_demarc;

out_of_your_box = FALSE;

if (which_test == STERIC_COST)
   {

/********************************************************************************/
/************* Test only for steric hinderence: between template and ************/
/************* pore and internally within template                   ************/
/********************************************************************************/

     if (!pbc)
       {
         out_of_your_box= box_test( p_box_limits, guest_ptrs, num_guests,
                                    p_guest_demarc );
       }

if (DEBUG)
  {  
     if (out_of_your_box)
       {
         printf("DB >> Your out of your bleeding box\n");
       }
     else
       {
         printf("DB >> Your NOT out of your bleeding box\n");
       }
  }

     *p_cc_tot = 0.0;
     if (!out_of_your_box)
       {
          calculate_energy(p_pore, num_pore_atoms,
                           guest_ptrs, num_guests, p_guest_demarc,
                           p_molmol_list, num_molmol_list,
                           p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                           p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                           have_comb_rules, is_empty_pore);

//          printf("DEBUG>> Back in cost_function 1\n");

/**************************************************************/
/****** Calculate the total of the close contact list *********/
/**************************************************************/

          p_demarc= p_guest_demarc;
          for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
            {
              if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

                for (iatom=0; iatom<= p_demarc->end; iatom++)   /* loop over this guests atoms */
                                         *p_cc_tot+= closest_contact[iatom];

                p_demarc++;
            }
           

//          if (DEBUG) printf("DEBUG>> Calculating close_contact total as %10.6f\n",  *p_cc_tot);


          if (!flags_set)
            {
              inter_guest_bump = bump_check(guest_ptrs, num_guests, p_guest_demarc,
                                            p_molmol_list, num_molmol_list);
            }
          else 
            {
              inter_guest_bump = bump_check_with_flags(guest_ptrs, num_guests, p_guest_demarc, 
                                                       p_molmol_list, num_molmol_list, p_flag_list);
            }
       }

DEBUG=TRUE;
       if (DEBUG)
           {
             
             for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
               {
                 if (interaction_energy[imol].steric) printf("FAIL on steric\n");
               }
             if (inter_guest_bump) printf("...................FAIL on inter_guest_bump\n");
             if (out_of_your_box) printf("...................FAIL on out_of_your_box\n");
           }
DEBUG=FALSE;

    for (imol=0; imol< num_guests; imol++)  
                 if (interaction_energy[imol].steric) inter_guest_bump = TRUE;

    return (inter_guest_bump || out_of_your_box);
  }

/****************************************************************************************/
/*** Monte Carlo testing approach introduced for NON_BOND cost function *****************/
/*** Including intra-molecular terms in potential.                      *****************/
/****************************************************************************************/
 else if (which_test == NON_BOND_COST)
  {
//     printf("Cost function using non-bond\n");
     if (!pbc)
       {
         out_of_your_box= box_test( p_box_limits, guest_ptrs, num_guests,
                                    p_guest_demarc );
       }
 
     if (!out_of_your_box)
       {
//         printf("DEBUG>> checking for bumps....\n");
         if (!flags_set)
           {
//             printf("DEBUG>> no flags..............\n");
             inter_guest_bump = bump_check(guest_ptrs, num_guests, p_guest_demarc, 
                                           p_molmol_list, num_molmol_list);
           }
         else 
           {
//             printf("DEBUG>> with flags..............\n");
             inter_guest_bump = bump_check_with_flags(guest_ptrs, num_guests, p_guest_demarc,
                                                      p_molmol_list, num_molmol_list, p_flag_list);
           }

         if (!inter_guest_bump && symm_set)
           {
              inter_guest_bump = bump_check(guest_ptrs, num_guests, p_guest_demarc,
                                            p_molmol_list, num_molmol_list);
           }

//         printf("DEBUG>> Doing host...guest..energy....\n");

   /*    if (!inter_guest_bump ) */
   /*      { */
               calculate_energy(p_pore, num_pore_atoms, 
                                guest_ptrs, num_guests, p_guest_demarc,
                                p_molmol_list, num_molmol_list,
                                p_kvecs, p_kvec2, 
                                p_gvec2, num_kvecs, p_cos_sum, p_sin_sum,
                                p_need_grad, p_grad, have_comb_rules, is_empty_pore);

               if (need_intra) 
                 { 
//                     printf("..................lets get intra energy.........................\n");
                     calculate_intra_energy(guest_ptrs, num_guests, p_guest_demarc,
                                            p_angles_list_ptrs, p_num_angles_list, 
                                            p_torsions_list_ptrs, p_num_torsions_list, 
                                            p_vdw_list_ptrs, p_num_vdws_list);
                 }
 
//             if (DEBUG)
//               {
//                  for (imol=0; imol< num_guests; imol++)              /* loop over guest molecules */
//                    {
//                      printf("In cost_function intra terms for molecule %d are:\n", imol);
//                      printf("Bond Stretch : %10.6f\n", intra_energy[imol].stretch);
//                      printf("Angle bend   : %10.6f\n", intra_energy[imol].angle  );
//                      printf("Torsion      : %10.6f\n", intra_energy[imol].torsion);
//                      printf("Van der Waals: %10.6f\n", intra_energy[imol].vdw);
//                      printf("------------------------------------------\n");
//                      printf("Total        : %10.6f\n", intra_energy[imol].total);
//                      printf("------------------------------------------\n\n");
//                   }
//               }
               
       /*  } */
       }

//       if (inter_guest_bump) printf("rejecting because of bumps\n");
//       if (out_of_your_box ) printf("rejecting because of out of box\n");

       return (inter_guest_bump || out_of_your_box);
  }

 else if (which_test == ENERGY_COST)
  {
     if (!pbc)
       {
         out_of_your_box= box_test( p_box_limits, guest_ptrs, num_guests,
                                    p_guest_demarc );
         if (out_of_your_box) printf("DEBUG>> template out of box\n");
       }

      
     if (!out_of_your_box)
       {
         if (!flags_set)
           {
             inter_guest_bump = bump_check(guest_ptrs, num_guests, p_guest_demarc,
                                           p_molmol_list, num_molmol_list);

             if (DEBUG && inter_guest_bump) 
                {
                   printf("DEBUG>> old and new clash\n"); 
                   print_molecule( guest_ptrs[index], p_demarc->end, stdout, FALSE);
                }
           }
         else if (flags_set)
           {
/*              printf("Flags were set\n"); */
             inter_guest_bump = bump_check_with_flags(guest_ptrs, num_guests, p_guest_demarc,
                                                      p_molmol_list, num_molmol_list, p_flag_list);
           }
         else
           {
             inter_guest_bump = FALSE;
           }

         if (!inter_guest_bump ) calculate_energy(p_pore, num_pore_atoms,
                                                  guest_ptrs, num_guests, p_guest_demarc,
                                                  p_molmol_list, num_molmol_list,
                                                  p_kvecs, p_kvec2,
                                                  p_gvec2, num_kvecs, p_cos_sum, p_sin_sum,
                                                  p_need_grad, p_grad, have_comb_rules, is_empty_pore);
          printf("DEBUG>> Back in cost_function 3\n");
       }

       return (inter_guest_bump || out_of_your_box);

  }
 else
  {
    printf("ERROR : Non-existant option requested from cost_function: option %d\n"
                                                                       ,which_test);
    exit(1);
  }
return -1;
}

/*******************************************************************************/
/*** This routine connects a new fragment to the selected guest for a build ****/
/*** operation. At this point the program should have selected the fragment ****/
/*** to join and the structure placed at the top of the realloced array of  ****/
/*** the guest.                                                             ****/
/*** The routine should be passed:                                          ****/
/*** p_guest - a pointer to the current list of guest atoms for the one     ****/
/***           that is being used in the build.                             ****/
/*** p_old_demarc - The demarcation structure for the old guest             ****/
/*** p_new_demarc - The demarcation structure for the new guest             ****/
/*** place_for_one - The point at which the first H is removed              ****/
/*** place_for_two - The point at which the second H is removed             ****/
/*** p_guest_hyd_list, p_num_guest_hyds - H atom indices for old part of    ****/
/***                                      guest                             ****/
/*** p_frag_hyd_list, p_frag_hyd_partition - hydrogen list and partition    ****/
/***                                         for fragment list              ****/
/*** frag_picked - index of chosen fragment                                 ****/
/*******************************************************************************/
/*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

int join_frag_to_temp(atom *p_template, int *p_num_template_atoms,
                      atom *p_molecule, int num_molecule_atoms,
                      int *p_place_for_one, int *p_place_for_two,
                      int *p_template_hyd_list, int *p_num_template_hyds,
                      int *p_molecule_hyd_list, int num_molecule_hyds,
                      bond *p_forbidden_bond, int *p_temp_hyd_weights,
                      int *p_sum_temp_hyd_weights,
                      int increment_hydrogen_weight, int have_AB);

void connect_new_fragment(atom *p_guest, list_partition *p_old_demarc, list_partition *p_new_demarc,   
                          int *p_place_for_one, int *p_place_for_two,
                          int *p_guest_hyd_list, int *p_num_guest_hyds,
                          int *p_frag_hyd_list, list_partition *p_frag_hyd_partition,
                          atom_number *p_guest_types, 
                          int *p_num_guest_types,
                          int *p_guest_hyd_weights, int *p_sum_guest_hyd_weights,
                          int increment_hydrogen_weight, int *p_have_AB)
  {
#include "header.h"

    int iloop, index, end_now, num_new_frag;
    int *p_this_guest_H, *p_this_frag_H, *p_copy_frag_hyd_list;
    int *p_this_copy;

    atom_number *p_this_type;

//    printf("Arrived in connect_new_fragment\n");
//    printf("Old demarc: start %d end %d num %d\n", p_old_demarc->start,  p_old_demarc->end,  p_old_demarc->num);
//    printf("New demarc: start %d end %d num %d\n", p_new_demarc->start,  p_new_demarc->end,  p_new_demarc->num);

    num_new_frag=p_new_demarc->num - p_old_demarc->num;

//    printf("Number of atoms added : %d\n", num_new_frag);

//    printf("Guest has %d hydrogens:\n", *p_num_guest_hyds);

//    p_this_guest_H= p_guest_hyd_list;
//    for (iloop=0; iloop <= *p_num_guest_hyds; iloop++)
//      {
//         printf("%d -> %d  %s (elem %s)\n", iloop, *p_this_guest_H,
//                     (p_guest+*p_this_guest_H)->label,(p_guest+*p_this_guest_H)->elem);
//         p_this_guest_H++;
//      }

//    printf("Fragment has %d hydrogens:\n", p_frag_hyd_partition->num);

//    p_this_frag_H= p_frag_hyd_list;
//    for (iloop=0; iloop <=  p_frag_hyd_partition->num; iloop++)
//      {
//         index=p_old_demarc->num + *p_this_frag_H;
//         printf("%d -> %d index:%d %s (elem %s)\n", iloop, *p_this_frag_H, index,
//                                          (p_guest+index)->label, (p_guest+index)->elem);
//         p_this_frag_H++;
//      }

//    printf("Current form of guest has %d types:\n", *p_num_guest_types);
    
//    p_this_type= p_guest_types;
//     for (iloop=0; iloop < *p_num_guest_types; iloop++)
//      {
//        printf("%s  %d\n", p_this_type->atom_type, p_this_type->num);
//        p_this_type++;
//      } 

    end_now = p_old_demarc->end;

/***** p_old_demarc->num gives an address one after the end of the old guest list *****/
/***** this is the first address that the new fragment copy starts and so should  *****/
/***** provide the sub-routine with the correct start point for the fragment to   *****/
/***** be added                                                                   *****/

//printf("Off to join_frag_to_temp\n");

/** Send a copy of the frag_hyd_list so that the original does not get corrupted ***/

 p_copy_frag_hyd_list=malloc((p_frag_hyd_partition->num+1)*sizeof(int));

printf("malloced space for %d frag H atoms...\n",(p_frag_hyd_partition->num+1));

 p_this_frag_H= p_frag_hyd_list;
 p_this_copy= p_copy_frag_hyd_list;
 for (iloop=0; iloop <= p_frag_hyd_partition->num; iloop++)
    {
      *p_this_copy=*p_this_frag_H;
      p_this_frag_H++; p_this_copy++;
    }

 join_frag_to_temp(p_guest, &end_now, 
                   p_guest+p_old_demarc->num, num_new_frag, 
                   p_place_for_one, p_place_for_two,
                   p_guest_hyd_list, p_num_guest_hyds, 
                   p_copy_frag_hyd_list, p_frag_hyd_partition->num, 
                   &forbidden_bond[0], p_guest_hyd_weights, 
                   p_sum_guest_hyd_weights,
                   increment_hydrogen_weight, *p_have_AB);

/*** Finished with copy of fragment hydrogen list now ***/
  free(p_copy_frag_hyd_list);

  p_new_demarc->end = end_now; p_new_demarc->num = end_now+1;

//
// Adjust types for loss of H or HA/HB atoms
//
 printf("DEBUG>> On returning from join_frag_to_temp fag H-list is...\n");
 p_this_frag_H= p_frag_hyd_list;
 for (iloop=0; iloop < p_frag_hyd_partition->num; iloop++)
    {
      printf("%d => %d\n", iloop, *p_this_frag_H);
      p_this_frag_H++;
    }

//printf("Adjusting number of H-types....\n");

  p_this_type= p_guest_types;
  if (*p_have_AB)
   {
     for (iloop=0; iloop < *p_num_guest_types; iloop++)
      {
        if(strcmp(p_this_type->atom_type, "HA") == 0) (p_this_type->num)--;
        if(strcmp(p_this_type->atom_type, "HB") == 0) (p_this_type->num)--;
        p_this_type++;
      } 
   }
  else
   {
     for (iloop=0; iloop < *p_num_guest_types; iloop++)
      {
        if(strcmp(p_this_type->atom_type, "H") == 0) p_this_type->num = p_this_type->num-2;
        p_this_type++;
      } 
   }

return;
}


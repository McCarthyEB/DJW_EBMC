/*************************************************************************/
/* assemble_intra_lists.c                                                */
/* Assemble lists of all angles, torsions and vdws required for          */
/* intra-molecular potential calculation.                                */
/*                                                                       */
/* Started Dave Willock 19th September 2006                              */
/* Updated for guest pointer arrays Dave Willock Feb. 2013               */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "header.h"

int assemble_angle_list(atom *p_molecule, int num_atoms, angle_interact_list *p_angles_list,
                        int just_count);
           
int assemble_torsion_list(atom *p_molecule, int num_atoms, torsion_interact_list *p_torsions_list, 
                          links *p_links_list, int *p_num_links, int just_count);
           
int assemble_vdw_list(atom *p_molecule, int num_atoms, vdw_interact_list *p_vdw_list, int just_count);

void assemble_intra_lists(atom *guest_ptrs[], int num_guests, list_partition *p_guest_demarc, 
                          int num_symm_ops, angle_interact_list *p_angles_list[],
                          torsion_interact_list *p_torsions_list[], vdw_interact_list *p_vdw_list[],
                          links *p_links_list_ptrs[],
                          int *p_num_angles_list, int *p_num_torsions_list, int *p_num_vdws_list, 
                          int *p_num_links, int just_count) 
  {

int iii, ind, imol;
atom *p_molecule;

int idebug;

//printf("DEBUG>> Arrived in assemble_intra_lists with %d guests\n", num_guests);
//printf("DEBUG>> Arrived in assemble_intra_lists with %d symmetry operations\n", num_guests);

//if (just_count) printf("just counting this time\n");
//else 
//   {
//      printf("filling this time\n");
//      for (imol=0; imol<num_guests; imol++)
//        {
//          printf("For molecule %d:\n", imol);
//          printf("expecting %d angles\n", *(p_num_angles_list+imol));
//          printf("expecting %d torsions\n", *(p_num_torsions_list+imol));
//          printf("expecting %d vdws\n", *(p_num_vdws_list+imol));
//          printf("expecting %d allowed torsion links\n", *(p_num_links+imol));
//        }
//   }

for (imol=0; imol<num_guests; imol++)
   {
     ind=imol;
     if ( num_symm_ops >= 0 ) ind= imol*(num_symm_ops+2);
/*****************************************************************************/
/*** Angle list, note that lists for symmetry images will not be required ****/
/*****************************************************************************/

     if ( just_count || *p_num_angles_list >=0 )
           *p_num_angles_list = assemble_angle_list(guest_ptrs[ind], p_guest_demarc->end, 
                                                    p_angles_list[imol], just_count);

/*******************************************/
/*** Now do list for torsion potentials ****/
/*******************************************/
     if ( just_count || *p_num_torsions_list >=0 )
      {
       *p_num_links=-1;
       *p_num_torsions_list = assemble_torsion_list(guest_ptrs[ind], p_guest_demarc->end,
                                                    p_torsions_list[imol], 
                                                    p_links_list_ptrs[imol], p_num_links,
                                                    just_count);
      }

//       if (just_count) printf("Max_length of torsions list for molecule %d will be: %d\n", 
//                                                              imol, *p_num_torsions_list);
//       else if ( *p_num_torsions_list >=0 )  
//               printf("Compiled %d torsion potentials for intra energy of molecule %d.\n", 
//                                                              *p_num_torsions_list, imol);

/***************************************/
/*** Now do list for vdw potentials ****/
/***************************************/
     if ( just_count || *p_num_vdws_list >=0 )
      *p_num_vdws_list = assemble_vdw_list(guest_ptrs[ind], p_guest_demarc->end, 
                                                      p_vdw_list[imol], just_count);

//      if (just_count) printf("Max_length of vdws list will be: %d\n", *p_num_vdws_list);
//      else if (*p_num_vdws_list >=0 ) printf("Compiled %d vdw potentials for intra energy.\n", 
//                                                                                 *p_num_vdws_list);

      p_num_angles_list++;
      p_num_torsions_list++;
      p_num_vdws_list++;
      p_num_links++;
      p_guest_demarc++;
    }

   return;
 }

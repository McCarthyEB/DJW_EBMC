/************************************************************************************************/
/*** Routine to bring together atom list and related information copying for build, GCMC etc ****/
/*** arrays relative to molecule A are copied into those for molecule B.                     ****/
/*** All array sizes must be set prior to calling routine.                                   ****/
/*** Started August 2014, Dave Willock **********************************************************/
/************************************************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"

void copy_molecule_data( atom *p_B_molecule, atom *p_A_molecule,
                         list_partition *p_B_demarc,  list_partition *p_A_demarc, 
                         links *p_B_links, links *p_A_links,
                         int *p_num_B_links, int *p_num_A_links,
                         int *p_B_hyd_list, int *p_A_hyd_list, 
                         int *p_B_hyd_weights, int *p_A_hyd_weights, 
                         int *p_B_num_hyds, int *p_A_num_hyds, 
                         int *p_B_sum_hyd_weights, int *p_A_sum_hyd_weights, 
                         atom_number *p_B_types, atom_number *p_A_types,
                         int *p_B_num_types, int *p_A_num_types ) 
{
  int iatom, ilink, ihyd, itype;

/*** Copy over atom list details ****/
//  printf("Copying molecule data for build\n");

  *p_B_demarc= *p_A_demarc;

  for (iatom=0; iatom <=p_A_demarc->end; iatom++)
    {
       *p_B_molecule = *p_A_molecule;
       p_A_molecule++; p_B_molecule++;
    }

//  printf("Copied %d atoms\n", p_A_demarc->end);
                        
/*** Copy over links list details ****/

   *p_num_B_links= *p_num_A_links;
 
   for (ilink=0; ilink <=*p_num_A_links; ilink++)
     {
       *p_B_links= *p_A_links;
       p_B_links++; p_A_links++;
     }
//  printf("Copied %d links\n", *p_num_A_links);

/*** Copy over hydrogen index lists and weights ****/

   *p_B_num_hyds= *p_A_num_hyds;
   *p_B_sum_hyd_weights= *p_A_sum_hyd_weights; 

   for (ihyd=0; ihyd <=*p_A_num_hyds; ihyd++)
     {
       *p_B_hyd_list= *p_A_hyd_list;
       p_B_hyd_list++; p_A_hyd_list++;
       *p_B_hyd_weights= *p_A_hyd_weights;
       p_B_hyd_weights++; p_A_hyd_weights++;
     }

//  printf("Copied %d hyd indices and weights\n", *p_A_num_hyds);
/*** Copy over types list  ****/

   *p_B_num_types= *p_A_num_types;
   for (itype=0; itype <=*p_A_num_types; itype++)
     {
       *p_B_types= *p_A_types;
       p_B_types++; p_A_types++;
     }

//  printf("Copied %d types\n", *p_A_num_types);
return;
}

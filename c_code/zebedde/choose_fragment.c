/************************************************************/
/* This routine takes a single guest molecule and selects   */
/* the fragment to try adding next. Checks on:              */
/* concentration limits would not be exceeded               */
/* there are no forbidden bonds made                        */
/*                                                          */
/* Originally part of build.c                               */
/* Started DWL 27/11/94, DJW 4/95                           */
/* Code split                                               */
/* Oct. 2013 Dave Willock                                   */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void add_types( atom_number *p_types1, int *p_num_types1,
                atom_number *p_types2, int num_types2);

void copy_types( atom_number *p_to, int *num_to,
                 atom_number *p_from, int num_from );

int pick_frm_wgt_list( int sum_weights, int *p_weights, int *p_chosen_weight,
                       int num_in_list );

int test_concentrations( atom_number *p_types, int num_types,
                         atom_number *p_limits, int num_limits);

int bonds_possible( atom *p_guest,          int num_guest_atoms,
                    int  *p_guest_hyd_list, int num_guest_hyds,
                    atom *p_fragment,          int num_fragment_atoms,
                    int  *p_frag_hyd_list,     int num_frag_hyds,
                    bond *p_forbidden_bond,    int have_AB);

/************************************************************/
/** Note num_types for guest is now really number NOT *******/
/** upper array limit.                                *******/
/************************************************************/

int choose_fragment(atom *p_guest, int num_guest_atoms, int *frag_weights,
                    int *p_guest_hyd_list, int num_guest_hyds,
                    int *p_frag_hyd_list, list_partition *p_frag_hyd_partition,
                    atom *frag_lib, atom_number *p_guest_types, int num_guest_types,
                    int *p_num_new_guest_types,
                    atom_number *p_frag_types, list_partition *p_frag_types_list,
                    int have_AB)
{
#include "header.h"

int start_this, ilimit, iloop, frag_picked, num_frag_types;
int acceptable, ifrag, can_bond;
int start_partition, num_frag_hyds, chosen_frag_weight;
int local_frag_hyd_list[MAX_HYDS], local_frag_weights[MAXFRAGMENTS];
int *p_this_frag_hyds;

atom_number *p_this_frag_types, *p_types, new_guest_types[100];
atom_number *p_this_type;

list_partition *p_this_hyd_partition;

atom *p_frag_atom;

/* this will be passed the first (or more) atoms of a guest */

/* altered to use user defined weights DW April 19th 1995 */

/*******DEBUG printing **************************************************************/
//DEBUG=FALSE;
//if (have_conc_limits && DEBUG)
//  {
//    printf("Have to impose %d concentration limits\n\n", num_conc_limits);
//    
//    for (ilimit= 0; ilimit < num_conc_limits; ilimit++)
//       {
//         printf("%s => %d\n", atom_limit[ilimit].atom_type, atom_limit[ilimit].num);
//       }
//  }
//    printf("\n%d forbidden bonds:\n\n",num_forbidden_bonds);
//
//    for (ilimit= 0; ilimit < num_forbidden_bonds; ilimit++)
//       {
//         printf("%s => %s \n",forbidden_bond[ilimit].atom1,forbidden_bond[ilimit].atom2);
//       }
//    printf("\n");
//    printf("Current guest has %d types:\n\n",num_guest_types);
//
//    p_types= p_guest_types;
//    for (iloop=0; iloop < num_guest_types; iloop++)
//      {
//      printf("%s => %d\n", p_types->atom_type, p_types->num);
//      p_types++;
//    }
//  for ( iloop = 0; iloop < number_of_fragments; iloop++)
//    {
//      printf("Weight for fragment %d = %d\n", iloop, frag_weights[iloop]);
//    }
/*******END OF DEBUG printing *******************************************************/

/**** We should work with local copy of the fragment weights as fragments may go in and out ****/
/**** of favour as a build progresses. Dave Willock Feb. 2014.                              ****/

    for ( iloop = 0; iloop < number_of_fragments; iloop++)
      {
        local_frag_weights[iloop] =  frag_weights[iloop];
      }

/******************************************************************/
/* The weights of dis-allowed fragments are set to                */
/* zero then check on sum of weights if it is zero no further     */
/* additions are possible                                         */
/******************************************************************/
//printf("Structure check on types........DEBUG.....\n");
//for (ifrag= 1; ifrag <= number_of_fragments; ifrag++)
//  {
//     num_frag_types=  (p_frag_types_list+ ifrag)->num;
//   
//     printf("Fragment %d (weight %d) has %d types\n", ifrag,  frag_weights[ifrag-1], num_frag_types);
//  }

sum_frag_weights = 0;
for (ifrag= 1; ifrag <= number_of_fragments; ifrag++)
  {
    if (frag_weights[ifrag-1] > 0)
       {
          p_this_frag_types= p_frag_types+ (p_frag_types_list+ ifrag)->start;
          num_frag_types=  (p_frag_types_list+ ifrag)->num;

//        if (DEBUG)
//           {
//              printf("/****************************************************/\n");
//              printf("In choose_fragment fragment %d has %d types\n",ifrag,num_frag_types);
//              printf("types list indexing: start= %d, end= %d, num= %d\n\n",
//                                               (p_frag_types_list+ ifrag)->start,
//                                               (p_frag_types_list+ ifrag)->end,
//                                               (p_frag_types_list+ ifrag)->num);
//
//          }
/****************************************************************************/
/******Test that this fragment will not push the atom concentrations********/
/************************over the limits************************************/
/****************************************************************************/

          copy_types( &new_guest_types[0], p_num_new_guest_types,
                      p_guest_types, num_guest_types);

//          printf(".....Copied atom types:\n");
//          p_this_type= p_guest_types;
//          for (iloop=0; iloop < *p_num_new_guest_types; iloop++)
//            {
//              printf("%s => %d   current guest %s => %d\n", new_guest_types[iloop].atom_type,
//                                                            new_guest_types[iloop].num,
//                                                            p_this_type->atom_type,
//                                                            p_this_type->num );
//              p_this_type++;
//            }
//          printf("Now adding fragment types:\n");

          add_types( &new_guest_types[0], p_num_new_guest_types,
                     p_this_frag_types, num_frag_types);

          if (*p_num_new_guest_types >= 100)
            {
               printf("ERROR>> In choose_fragment.c the number of types exceeds the size of the temporary array\n");
               exit(0);
            }

//        if (DEBUG)
//          {
//               printf("\n");
//               printf("New guest would have %d types (before H deletion):\n\n",*p_num_new_guest_types);
//  
//              for (iloop=0; iloop < *p_num_new_guest_types; iloop++)
//                 {
//                   printf("%s => %d\n", new_guest_types[iloop].atom_type, 
//                                                       new_guest_types[iloop].num);
//                 }
//          }
                
          acceptable= test_concentrations(&new_guest_types[0], *p_num_new_guest_types,
                                          &atom_limit[0], num_conc_limits); 

/****************************************************************************/
/******Test that this fragment has suitable valances to bond with the ******/
/******************* current guest **************************************/
/****************************************************************************/

          start_this= member_start[ifrag];
          p_this_hyd_partition= p_frag_hyd_partition+ifrag;
          start_partition= p_this_hyd_partition->start;
          num_frag_hyds= p_this_hyd_partition->num;

          can_bond= bonds_possible( p_guest, num_guest_atoms,
                                    p_guest_hyd_list,  num_guest_hyds,
                                    &frag_lib[start_this], 
                                    number_of_members[ifrag],
                                    p_frag_hyd_list+start_partition, 
                                    num_frag_hyds,
                                    &forbidden_bond[0], have_AB);

          if (!can_bond && DEBUG)
            {
               printf("This fragment is FORBIDDEN to bind to this guest\n");
            }
  
         if (!acceptable && DEBUG)
            {
               printf("Fragment %d is unacceptable due to concentration limits.\n",ifrag);
            }
 
         if (!can_bond || !acceptable) local_frag_weights[ifrag-1]=0;

         sum_frag_weights += local_frag_weights[ifrag-1];
      }
  }

//if (DEBUG)
//  {
//     for ( iloop = 0; iloop < number_of_fragments; iloop++)
//       {
//         printf("Latest weight for fragment %d = %d\n", iloop, local_frag_weights[iloop]);
//       }
//     printf("New sum of weights            = %d\n",sum_frag_weights );
//  }

/**************************************************************************/
/***** Program drops out when sum of weights hits zero ********************/
/***** may wish to exit with guest write at this point *****************/
/**************************************************************************/

if (sum_frag_weights == 0) 
   {
     printf("None of the available fragments can be added to the current guest\n");
     printf("Building cannot continue: program exiting\n\n");
     fflush(stdout);
     exit(0);
   }

/*************************************************************************/
/*******************now pick fragment from those allowed *****************/
/*************************************************************************/

//   if (DEBUG) 
//     {
//       printf("SI DEBUG>> Choosing fragment using pick_frm_wgt_list\n");
//       printf("SI DEBUG>> sending sum_weights=%d num=%d \n", sum_frag_weights, num_frag_weights-1);
//     }
   frag_picked= 1+pick_frm_wgt_list(sum_frag_weights, &local_frag_weights[0], 
                                    &chosen_frag_weight, num_frag_weights-1);

//   printf("Choose_fragment selected fragment %d of %d returning\n", frag_picked, number_of_fragments);

   return frag_picked;
  }

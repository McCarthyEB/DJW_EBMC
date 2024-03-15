
/*************************************************************************/
/* Order_double_pot.c                                                    */
/* Orders a potential list that depends on two atom types so that they   */
/* can be referenced by the single indexes of each atom.                 */
/* re-orders potential list so that A B C D atom types occur in order:   */
/* A:A A:B A:C A:D B:B B:C B:D C:C C:D D:D                               */
/* any missing parameters will be set to zero!                           */
/*                                                                       */
/* Started Dave Willock 8/2/99                                           */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void order_double_pot(potential_list *p_list_potent, int *p_num_in_list, 
                                              types *p_pot_types, int *p_num_types, int is_ordered)
{
#include "header.h"
int ipot, pot_new, num_now;
int itype, jtype, icount;
types *p_this_type, *p_that_type;
potential_list *p_this_potent;
potential_list orig_pots[MAX_POTS];

/****** Loop through and count number of types of potentials ******/
printf("Entered order_double_pot with %d types\n", *p_num_types);

*p_num_types=-1;
p_this_potent= p_list_potent-1;
for (ipot=0; ipot <= *p_num_in_list; ipot++)
  {
     pot_new=TRUE;
     p_this_potent++;
     orig_pots[ipot]= *p_this_potent;

     p_this_type=p_pot_types;
     for (itype = 0; itype  <= *p_num_types; itype++)
       {
         if ( strcmp(p_this_type->name, p_this_potent->pot) == 0 ) pot_new=FALSE;
         p_this_type++;
       }
     if (pot_new)
       {
         (*p_num_types)++;
         strcpy(p_this_type->name, p_this_potent->pot);
       }

     pot_new=TRUE;
     p_this_type=p_pot_types;
     for (itype = 0; itype  <= *p_num_types; itype++)
       {
         if ( strcmp(p_this_type->name, p_this_potent->pot2) == 0 ) pot_new=FALSE;
         p_this_type++;
       }
     if (pot_new)
       {
         (*p_num_types)++;
         strcpy(p_this_type->name, p_this_potent->pot2);
       }
  }

/****** Check all is well ******************************************/

printf("Found the following atom types can be involved in double parameter bonds\n");

p_this_type=p_pot_types;
for (itype = 0; itype  <= *p_num_types; itype++)
  {
     printf("%d %s \n", itype, p_this_type->name);    
     p_this_type++;
  }

p_this_potent= p_list_potent-1;
num_now= *p_num_in_list;
*p_num_in_list=-1;
p_this_type=p_pot_types;

/**** Double loop over potential types with shortened inner loop to avoid double counts *******/

icount=0;
for (itype = 0; itype  <= *p_num_types; itype++)
  {
    p_that_type=p_pot_types+itype;
    for (jtype = itype; jtype  <= *p_num_types; jtype++)
      {
          printf("Will set potential %d  for %s  %s \n", icount, p_this_type->name,  p_that_type->name);
          icount++;
          pot_new=TRUE;

/**** Get available parameters from original list and fill others with zeroes ****/

          for (ipot=0; ipot <= num_now; ipot++)
            {
               if (   strcmp(p_this_type->name, orig_pots[ipot].pot2) == 0 &&
                      strcmp(p_that_type->name, orig_pots[ipot].pot) == 0   )
                  {
                     p_this_potent++;
                     ++(*p_num_in_list);
                     pot_new=FALSE;
                    
                     strcpy(p_this_potent->pot2,  p_this_type->name);
                     strcpy(p_this_potent->pot, p_that_type->name);

/**************************************************************************/
/*** If the original is ordered differently need to swap params over for **/
/*** bond increment case                                                 **/
/**************************************************************************/
                     if (is_ordered)
                       {
                         p_this_potent->a = orig_pots[ipot].b;
                         p_this_potent->b = orig_pots[ipot].a;
                         p_this_potent->c = orig_pots[ipot].c;
                       }
                     else
                       {
                         p_this_potent->a = orig_pots[ipot].a;
                         p_this_potent->b = orig_pots[ipot].b;
                         p_this_potent->c = orig_pots[ipot].c;
                       }
                  }
                else if  ( strcmp(p_this_type->name, orig_pots[ipot].pot) == 0 &&
                           strcmp(p_that_type->name, orig_pots[ipot].pot2) == 0  )
                  {
                     p_this_potent++;
                     ++(*p_num_in_list);
                     pot_new=FALSE;

                     strcpy(p_this_potent->pot,  p_this_type->name);
                     strcpy(p_this_potent->pot2, p_that_type->name);
                     p_this_potent->a = orig_pots[ipot].a;
                     p_this_potent->b = orig_pots[ipot].b;
                     p_this_potent->c = orig_pots[ipot].c;

                  }
            }
          if ( pot_new )
            {
               p_this_potent++;
               ++(*p_num_in_list);
               strcpy(p_this_potent->pot,  p_this_type->name);
               strcpy(p_this_potent->pot2, p_that_type->name);
               p_this_potent->a =0.0;
               p_this_potent->b =0.0;
               p_this_potent->c =0.0;
            }
         p_that_type++;
      }
    p_this_type++;
  }

p_this_potent= p_list_potent-1;
for (ipot=0; ipot <= *p_num_in_list; ipot++)
  {
    p_this_potent++;
    printf("%s %s %10.6f %10.6f %10.6f \n", p_this_potent->pot,  p_this_potent->pot2,
                                  p_this_potent->a, p_this_potent->b, p_this_potent->c);
  }
return;
}



/************************************************************/
/* shuffle_hyd_list                                         */
/* shuffle hydrogen list to cover up for lost hydrogens     */
/* Started DJW 14/7/95                                      */
/*                                                          */
/* have_deactivated flag indicates that some hydrogens have */
/* been marked with an atom index in the list of -1.        */
/* These also need to be removed and the list closed up.    */
/* Added Nov 06 Dave Willock                                */
/*                                                          */
/************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "header.h"

int shuffle_hyd_list(int *p_list_start, int place_for_one,
                     int place_for_two, int num_in_list,
                     int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights,
                     int have_deactivated)
{
int ihyd, jhyd, ientry, index;

int *p_list, *p_new_list, *p_temp_list;
int *p_new_weights, *p_this_weight, *p_temp_weight;
int *p_this_entry, new_num_hyds;
int iii, num;

/********************Remove the template's hydrogen from the hydrogens list********/
/********************list will index the current list and new_list the updated ****/

p_list= p_list_start-1;
p_new_list= p_list;

p_this_weight= p_temp_hyd_weights-1;
p_new_weights= p_this_weight;

new_num_hyds= -1;

for (ihyd= 0; ihyd <= num_in_list; ihyd++)
  {
    p_list++;
    p_this_weight++;

/**** If these are not the ones to be removed copy over into new array positions ***/

    if ( *p_list != place_for_one && *p_list != place_for_two)
      {
        p_new_list++;
        *p_new_list= *p_list;

        p_new_weights++;
        *p_new_weights= *p_this_weight;

        new_num_hyds++;
      }

/**** If these are the ones to be removed re-adjust the weight total ***/
    else
      {
        *p_sum_temp_hyd_weights -= *(p_temp_hyd_weights+ihyd);
      }
  }

/**** For the H atoms that are being removed the real atom list will also ***/
/**** be compacted. This means that indicies of all higher hydrogens need ***/
/**** adjusting, so the content of this list needs decrementing           ***/

 p_this_entry= p_list_start-1;

 for (ientry=0; ientry <= new_num_hyds; ientry++)
   {
     p_this_entry++;

     index= *p_this_entry;

     if (index > place_for_two)
       {
         (*p_this_entry) -= 2;
       }
     else if (index > place_for_one)
       {
         (*p_this_entry) --;
       }
   }

/***************************************************************/
/*** Deal with deactivation ************************************/
/*** now new_num_hyds will contain the number in the list ******/
/*** with the bonding H atoms deleted                     ******/
/*** The deactivated hydrogens are still present in the   ******/
/*** template list and so simple removal from this list   ******/
/*** is all that is required.                             ******/
/***************************************************************/

p_list= p_list_start-1;
p_new_list= p_list;

p_this_weight= p_temp_hyd_weights-1;
p_new_weights= p_this_weight;

num= new_num_hyds;

     if (DEBUG)               
       {
         printf("Starting deactivation list is %d:\n", new_num_hyds);
         printf("Weight total is %d:\n", *p_sum_temp_hyd_weights);
         p_temp_list= p_list_start;
         for (iii= 0; iii <= new_num_hyds; iii++)
           {
             printf("%d ",*p_temp_list); 
             p_temp_list++;
           }
         printf("\n");
       }

jhyd = -1;
for (ihyd= 0; ihyd <= num; ihyd++)
  {
    jhyd++;
    p_list++; 
    p_this_weight++;

    if (*p_list == -1)
      {
         if (DEBUG) printf("Removing\n");
/********************************************************/
/* collapse list on top of this one                     */
/* jhyd counts current position accouting for deletions */
/********************************************************/
         p_temp_list= p_list;
         p_temp_weight = p_this_weight;

/**** Adjust total weight ****/
         *p_sum_temp_hyd_weights -= *p_this_weight;

         for (iii= jhyd; iii < new_num_hyds; iii++)
           {
             *p_temp_list = *(p_temp_list+1);  
             *p_temp_weight = *(p_temp_weight+1);  
             p_temp_list++;
             p_temp_weight++;
           }
         new_num_hyds--;
         p_this_weight--;
         p_list--;
         jhyd--;

      if (DEBUG)
       {
         printf("After deactivation list is %d:\n", new_num_hyds);
         p_temp_list= p_list_start;
         for (iii= 0; iii <= new_num_hyds; iii++)
           {
             printf("%d ",*p_temp_list); 
             p_temp_list++;
           }
         printf("\n");
       }
      }
    else
      {
        p_new_list++;
        p_new_weights++;

        if (DEBUG) printf("Keeping\n");

        *p_new_list = *p_list;
        *p_new_weights = *p_this_weight;
      }
  }
 
return new_num_hyds; 
}


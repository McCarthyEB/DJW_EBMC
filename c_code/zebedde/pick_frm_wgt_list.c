#include <stdio.h>
#include <stdlib.h>

double real_random(int done);

/****************************************************************/
/******* choose an option from a weighted list          *********/
/******* using a random number                          *********/
/*******                                                *********/
/******* Dave and Dewi May 1995                         *********/
/****************************************************************/

int pick_frm_wgt_list( int sum_weights, int *p_weights, int *p_chosen_weight,
                       int num_in_list )
{

int item_picked, random_num, top_weight;
double rando;

rando=real_random(1);
random_num = rando*sum_weights;

item_picked=0;
top_weight= *p_weights;

/**********************************************************************/
/* num_in_list will be sent as the highest index of the weight array  */
/* Dave Willock, June 2012.                                           */
/**********************************************************************/

// printf("Picking from weighted list with %d entries\n", num_in_list); 
// printf("Selecting random_num = %d, rando=%10.6f, top_weight starts at %d\n", random_num, rando, top_weight); 

for (item_picked = 0; item_picked <= num_in_list; item_picked++)
   {
     if ( random_num < top_weight )  
       {
         *p_chosen_weight= *p_weights;
         return item_picked;
       }

     p_weights++;
     top_weight+= *p_weights; 

   }

*p_chosen_weight= *p_weights;
return item_picked;
}

/************************************************************/
/* shuffle_links                                            */
/* shuffle links to cover up for lost hydrogens             */
/*                                                          */
/* Started DJW 14/7/95                                      */
/************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void shuffle_links(links *p_link_atoms, int place_for_one,
                        int place_for_two, int num_links)
{
int ilink, index;

links *p_this_link;

 p_this_link= p_link_atoms-1;
 for (ilink=0; ilink <= num_links; ilink++)
   {
     p_this_link++;

     index= p_this_link->start;

     if (index > place_for_two)
       {
         (p_this_link->start) -= 2;
       }
     else if (index > place_for_one)
       {
         (p_this_link->start) --;
       }

     index= p_this_link->end;

     if (index > place_for_two)
       {
          (p_this_link->end) -= 2;
       }
     else if (index > place_for_one)
       {
          (p_this_link->end) --;
       }
   }
return; 
}


/***************************************************************/
/*** Routine to give the double pot list index from ************/
/*** Single list indecies                           ************/
/***************************************************************/

#include <stdio.h>
#include "global_values.h"

int index_double_pot(int index_1, int index_2, int num_types, int *p_1_then_2)
  {
    int i, index;

    if ( index_1 < index_2 )
        {
           index= index_1*(num_types+1);
           for (i=1; i < index_1; i++) index-=i;
           index += index_2-index_1;
           *p_1_then_2= TRUE;
        }
      else
        {
           index= index_2*(num_types+1);
           for (i=1; i < index_2; i++) index-=i;
           index += index_1-index_2;
           *p_1_then_2= FALSE;
        }

    return index;
  }

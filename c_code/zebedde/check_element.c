/**************************************************************************/
/* check_element : checks string against periodic table for valid element */
/* Dewi 11/7/95                                                           */
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "data.h"

int check_element(char *p_element)

{
  int i;
  int found_elem;	

printf("checking element - ");
  i=-1;
  found_elem = FALSE;
  while ( (found_elem == FALSE) && (i < NUM_ELEMENTS) )
      {
      i++;
      if (strcmp(p_element, period_table[i].elem) == 0)
          {
          found_elem = TRUE;
          }
      }
printf("i=%i, %s %s\n", i,p_element,period_table[i].elem);

return(found_elem);
}

/***************************************************/
/*                                                 */
/* Test if this template                           */
/* exceeds the set concentration limits.           */
/*                                                 */
/* Dave Willock August 1995                        */
/***************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

int test_concentrations( atom_number *p_types, int num_types,
                         atom_number *p_limits, int num_limits)
{
 int ilimits, got_type, itype;

 char *p_label_limit, *p_this_label;

 atom_number *p_this_type;

/**************************************************************************/
/************ Loop over the limits list    ********************************/
/**************************************************************************/

 p_limits--;
 for (ilimits = 0; ilimits < num_limits; ilimits++)
   {
     p_limits++;
     p_this_type= p_types-1;
     p_label_limit= &(p_limits->atom_type[0]);
     got_type= FALSE;

/**************************************************************************/
/************ Check to see if have any of these in the template ***********/
/**************************************************************************/

     itype=0;
     while ( itype <= num_types && !got_type)
        {
           itype++;
           p_this_type++;
           p_this_label= &(p_this_type->atom_type[0]);

           if ( compare_strings( p_label_limit , p_this_label )) got_type= TRUE; 
        }

/**************************************************************************/
/******** check concentrations of this type *******************************/
/**************************************************************************/

     if (got_type && p_this_type->num >= p_limits->num) return FALSE;
  }
return TRUE;
}


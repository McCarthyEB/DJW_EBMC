/***************************************************/
/*                                                 */
/* Compare two atom element type strings to see if */
/* their bonding is forbidden                      */
/*                                                 */
/* Dave Willock August 1995                        */
/***************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

int forbid_bond( char *p_type1,          char *p_type2,
                 bond *p_forbidden_bond)
{
#include "header.h"
int itest, forbidden, this_way, that_way;

 itest=0;
 forbidden= FALSE;

 while ( itest < num_forbidden_bonds && !forbidden)
   {
      this_way= compare_strings(p_forbidden_bond->atom1, p_type1)
              &&compare_strings(p_forbidden_bond->atom2, p_type2);

      that_way= compare_strings(p_forbidden_bond->atom2, p_type1)
              &&compare_strings(p_forbidden_bond->atom1, p_type2);

      if (this_way || that_way) forbidden= TRUE;

      p_forbidden_bond++;
      itest++;
   }  
 return forbidden;
}

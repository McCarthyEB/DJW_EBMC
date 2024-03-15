/***************************************************/
/*                                                 */
/* Make a list of the element types which neighbour*/
/* hydrogen in this molecule, return their number  */
/*                                                 */
/* Dave Willock August 1995                        */
/***************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

void copy_strings( char to[], char from[] );

int find_hyd_neighbours(atom *p_molecule, int num_atoms,
                        int  *p_hyd_list, int num_hyds,
                        types *p_types_list)
{
#include "header.h"
 atom *p_hyd, *p_neigh;

 types *p_this_type;

 int itype, num_neigh_types, got_type, ihyd;

/**************************************************************************/
/************ loop over the molecule's hydrogens  *************************/
/************ and compare with those found so far *************************/
/**************************************************************************/
 num_neigh_types=-1;

 for (ihyd=0; ihyd <= num_hyds; ihyd++)
   {
     p_hyd= p_molecule+ *(p_hyd_list +ihyd);
     p_neigh= p_molecule+ p_hyd->neighb[0];

     got_type= FALSE;
 
     p_this_type= p_types_list;
     for (itype=0; itype <= num_neigh_types; itype++)
        {
          if (compare_strings(&(p_this_type->name[0]), &(p_neigh->elem[0])))
                                                                     got_type= TRUE;
          p_this_type++;
        }


     if (!got_type) 
       {
         num_neigh_types++;
         copy_strings( &(p_this_type->name[0]), &(p_neigh->elem[0]));
       }
   }

return num_neigh_types;
}

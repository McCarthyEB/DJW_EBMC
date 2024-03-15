/***************************************************/
/*                                                 */
/* Test to see if a bond can be formed between this*/
/* template and fragment which is not forbidden    */
/*                                                 */
/* Dave Willock August 1995                        */
/*                                                 */
/* Updated to allow for have_AB cases.             */
/* Dave Willock Feb. 2014                          */
/***************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

int forbid_bond( char *p_type1,          char *p_type2,
                 bond *p_forbidden_bond);

void copy_strings( char to[], char from[] );

int find_hyd_neighbours(atom *p_molecule, int num_atoms,
                        int  *p_hyd_list, int num_hyds,
                        types *p_types_list);
 
int bonds_possible( atom *p_template,          int num_template_atoms,
                    int  *p_template_hyd_list, int num_template_hyds,
                    atom *p_fragment,          int num_fragment_atoms,
                    int  *p_frag_hyd_list,     int num_frag_hyds,
                    bond *p_forbidden_bond,    int have_AB)
{
#include "header.h"

 types temp_neigh_types[100], frag_neigh_types[100];

 int num_temp_neigh_types, num_frag_neigh_types;
 int itemp_neigh, ifrag_neigh, ok_now;
 int temp_has_A, temp_has_B, frag_has_A, frag_has_B;
 int h_index, ihyd;

 char lab_HA[3], lab_HB[3];

/************************************************************************/
/********** get list of hydrogen neighbour types for ********************/
/********** template and fragment                    ********************/
/************************************************************************/

num_temp_neigh_types= find_hyd_neighbours(p_template, num_template_atoms,
                                          p_template_hyd_list, num_template_hyds,
                                          &temp_neigh_types[0]);

num_frag_neigh_types= find_hyd_neighbours(p_fragment, num_fragment_atoms,
                                          p_frag_hyd_list, num_frag_hyds,
                                          &frag_neigh_types[0]);

/************************************************************************/
/********** loop over pairs of the neighbours testing *******************/
/********** against the forbidden bonds list          *******************/
/********** Once we find a pair that are allowed to   *******************/
/********** bond we are OK from this respect....      *******************/
/************************************************************************/

ok_now=FALSE;
for (itemp_neigh= 0; itemp_neigh <= num_temp_neigh_types; itemp_neigh++)
   {
      for (ifrag_neigh=0; ifrag_neigh <= num_frag_neigh_types; ifrag_neigh++)
         {
            if (
                  !forbid_bond(temp_neigh_types[itemp_neigh].name,
                               frag_neigh_types[ifrag_neigh].name,
                               p_forbidden_bond)
                                                           ) ok_now=TRUE;
         }
   }

/*** Test that there are complementary HA/HB pairs on the guest and fragment ****/

//if (ok_now) printf("OK so far....\n");

temp_has_A=FALSE; temp_has_B=FALSE; frag_has_A=FALSE; frag_has_B=FALSE; 
if (ok_now && have_AB)
  {
    sprintf(lab_HA,"HA");sprintf(lab_HB,"HB");
    for (ihyd=0; ihyd<= num_template_hyds; ihyd++)
      {
         h_index= *(p_template_hyd_list+ihyd);
         if (compare_strings( (p_template+h_index)->elem, lab_HA )) temp_has_A=TRUE;
         if (compare_strings( (p_template+h_index)->elem, lab_HB )) temp_has_B=TRUE;
      }
    for (ihyd=0; ihyd<= num_frag_hyds; ihyd++)
      {
         h_index= *(p_frag_hyd_list+ihyd);
         if (compare_strings( (p_fragment+h_index)->elem, lab_HA )) frag_has_A=TRUE;
         if (compare_strings( (p_fragment+h_index)->elem, lab_HB )) frag_has_B=TRUE;
      }

//    if (temp_has_A) printf("Template has HA\n");
//    if (temp_has_B) printf("Template has HB\n");
//    if (frag_has_A) printf("Fragment has HA\n");
//    if (frag_has_B) printf("Fragment has HB\n");

    if ( !(( temp_has_A && frag_has_B ) || ( temp_has_B && frag_has_A )) ) ok_now=FALSE;
  }

return ok_now;
}

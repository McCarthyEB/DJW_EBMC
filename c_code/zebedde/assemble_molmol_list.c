/***************************************************************/
/* Code to assemble the list of molecule-molecule interactions */
/* this assembles all unique molecule pairs including self     */
/* interaction terms.                                          */
/*  imol/jmol  indexes the molecules numbered as if no symm.   */
/*  ind /jnd      are the actual indexes to use in the real    */
/*       list that combines the molecules and their images.    */
/* Started Dave Willock March 2013                             */
/***************************************************************/
#include <stdio.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

int assemble_molmol_list(interaction_indices *p_molmol_list, 
                         int num_guests, int num_expected)
{
#include "header.h"
int ind,jnd, tot_num;
int imol, jmol, isymm, jsymm;
int count;

//printf("Assembling molmol list:\n");
count=0;
for (imol=0; imol < num_guests; imol++)
  {
     if (symm_set) ind=imol*(num_symm_ops+2);
                                             else ind=imol;

     for (isymm=0; isymm<= num_symm_ops+1; isymm++)              /* loop over symmetry operations */
       {
         for (jmol=imol; jmol< num_guests; jmol++)              /* loop over guest molecules */
          {
            if (symm_set) jnd=jmol*(num_symm_ops+2);
                                             else jnd=jmol;

              for (jsymm=0; jsymm<= num_symm_ops+1; jsymm++)              /* loop over symmetry operations */
                {
                  if (jnd >= ind)
                    {
//                    printf("%d: imol=%d, ind=%d, jmol=%d, jnd=%d\n",       
//                                              count, imol, ind, jmol, jnd); 
                      p_molmol_list->imol=imol;
                      p_molmol_list->jmol=jmol;
                      p_molmol_list->ind =ind;
                      p_molmol_list->jnd =jnd;
                      p_molmol_list++;
                      count++;
                    }
                  jnd++;
                }
          }
        ind++;
      }
  }
if (count==num_expected)
                  return count;
else
   return -1;
}

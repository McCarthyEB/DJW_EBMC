/****************************************************************************/
/* print_peek_file.c : Routine to write a biosym.car file of the current    */
/*                     template when requested                              */
/* Started 29/8/96 DWL                                                      */
/* Notes: Currently file_fp is set in make_a_template as a constant from    */
/*        global_values.h                                                   */
/* Updates for multiple molecules March 2013, Dave Willock                  */
/****************************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);
void print_biosym_header(FILE *file_fp, int pbc_flag);
void print_frame_header(char *p_title, int num_frame, FILE *file_fp);
void car_end(FILE *fp);

void print_peek_file(atom *guest_ptrs[], int num_guests, 
                     list_partition *p_guest_demarc, int counter, FILE *peek_fp) 

{
#include "header.h"

char title_line[80];
int imol, isymm, index;

print_biosym_header(peek_fp, FALSE);         /* not PBC */
sprintf(title_line, "ZEBEDDE PEEK File, Attempt Number %d", counter);
print_frame_header(&title_line[0], -1, peek_fp);  /* not an archive file */

for (imol=0; imol < num_guests; imol++)
  {
    if (symm_set) index=imol*(num_symm_ops+2);
                                             else index=imol;

    for (isymm=0; isymm<=num_symm_ops+1; isymm++)
      {
        print_molecule(guest_ptrs[index], p_guest_demarc->end, peek_fp, FALSE);
        fprintf(peek_fp,"end\n");
        index++;
      }
    p_guest_demarc++;
  }
fprintf(peek_fp,"end\n");
return;
}


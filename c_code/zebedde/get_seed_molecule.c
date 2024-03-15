/*****************************************************************************/
/** get_seed_molecule: read seed from supplied file (single molecule)  *******/
/***                   dwl/djw 28/2/96                                 *******/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void get_seed_molecule(atom *p_original_template,
                       int *p_num_original_template_atoms,
                       int frame_number, int just_count)

{
#include "header.h"

int ineigh, i;
int got_date;
int icount=0;

atom *p_atom;

if (just_count)
  {
    printf("DB>> start seed read, just counting\n");
  }
else
  {
    printf("DB>> start seed read, getting co-ordinates\n");
  }

  if (!(seed_fp = fopen(seed_file, "r")))
        {
                printf("error opening seed file %s\n", seed_file);
                exit(1);
        }
  else
        {
        /* READ IN THE fragment COORDINATES */

        fgets(dummy_head,BUFFER,seed_fp);
        fgets(dummy_head,BUFFER,seed_fp);
		for (i=0; i<=frame_number; i++)
			{
			got_date = FALSE;
                        icount=0;
			while (!got_date)
				{
        		fgets(dummy_head,BUFFER,seed_fp);
				if (strstr(dummy_head, "!date") != NULL) got_date = TRUE;
				if (strstr(dummy_head, "!DATE") != NULL) got_date = TRUE;
                                 if (++icount >10)
                                  {
                                    printf("ERROR: Date line missing from the seed file\n");
                                    printf("ERROR: Please add something like:\n");
                                    printf("\n!DATE Sat Aug 17 10:21:28 2013\n\n");
                                    printf("ERROR: Under the title line and try again....\n");

                                    exit(0);
                                  }
				}
			}
/***** alteration 6th April 1995 DW: make template atoms reference to zero*****/
    (*p_num_original_template_atoms) = -1;

    p_atom = p_original_template;
      
    while(fgets(buffer,BUFFER,seed_fp) != NULL)
      {
      if (strstr(buffer,"end") != NULL) break;
      else
       {
        if (just_count)
          {
             (*p_num_original_template_atoms)++;
          }
        else
          {
             sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf", &(p_atom->label[0]), 
	            		                     &p_atom->x, 
                                                     &p_atom->y, 
		      	                             &p_atom->z, 
                                                     &(p_atom->group[0]),
			                             &(p_atom->group_no[0]), 
                                                     &(p_atom->pot[0]), 
			                             &(p_atom->elem[0]), 
                                                     &p_atom->part_chge);

		/* set up neighbour arrays */
             p_atom->num_neigh=-1;
             for (ineigh=0; ineigh < 4; ineigh++)
                                    p_atom->neighb[ineigh]=-1;

             p_atom++;
             (*p_num_original_template_atoms)++;
          }
        }
      }
     
    fclose(seed_fp);
    }

return;
}

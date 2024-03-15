/*****************************************************************************/
/** get_pre_anneal_template : get pre annealed template from pre_anneal.car***/
/***  (copied from get_seed_molecule) KJ 28/11/06                      *******/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void get_pre_anneal_template(atom *p_template,
                       int *p_num_template_atoms,
					   int frame_number)

{
#include "header.h"

int ineigh, i;
int got_date;

atom *p_atom;

printf("DB>> starting reading pre anneal template file\n");
  if (!(pre_opt_fp = fopen("pre_anneal.car", "r")))
        {
                printf("error opening pre_anneal.car file \n");
                exit(1);
        }
  else
        {
        /* READ IN THE fragment COORDINATES */

/***** alteration 6th April 1995 DW: make template atoms reference to zero*****/
       (*p_num_template_atoms) = -1;
		p_atom = p_template;
      
    while(fgets(buffer,BUFFER,pre_opt_fp) != NULL)
      {
      if (strstr(buffer,"end") != NULL) break;
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
		(*p_num_template_atoms)++;
           }
      }
     
	fclose(pre_opt_fp);
    }

return;
}

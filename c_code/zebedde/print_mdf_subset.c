/****************************************************************/
/* print_mdf_subset.c : prints a mdf for a molecule             */
/*                                                              */
/****************************************************************/


#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_mdf_subset(char *p_molname, atom *p_molecule, int num_atoms, 
	   	      int start, FILE *fp, int D_to_H, int have_subset)
{

int iatom1,iatom2;

char this_group[10], that_group[10];
atom *p_atom1;

if (!have_subset)  fprintf(fp, "\n\n");

/****** If required print out atom subset for fixing ********/
/****** Also need to put on ending in this case      ********/

if (have_subset)
  {
    fprintf(fp, "\n#atomset                                                                     \n\n\n"); 
    fprintf(fp, "@list subset MDF_FIXED_ATOMS_SUBSET\n\n");

    p_atom1= p_molecule;
    sprintf(this_group, "%s_%s", p_atom1->group, p_atom1->group_no);
    fprintf(fp, "pore:%s:", this_group);

    iatom2=0;
    for (iatom1=0; iatom1<= num_atoms; iatom1++)
      {
        sprintf(that_group, "%s_%s", p_atom1->group, p_atom1->group_no);
        if (strcmp(this_group, that_group) == 0 )
          {
            iatom2++;
            if (iatom2 < 17)
              {
                 fprintf(fp, "%s ",p_atom1->label); 
              }
            else
              {
                 fprintf(fp, "\n%s ",p_atom1->label); 
                 iatom2=0;
              }
          }
        else
          {
            iatom2=0;
            strcpy( this_group, that_group );
            fprintf(fp, "\n%s:%s ",this_group, p_atom1->label);
          }
        p_atom1++;
      }

    fprintf(fp, "\n#end\n");
  }

return;
}


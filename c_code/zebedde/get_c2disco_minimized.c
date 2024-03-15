/**************************************************************************/
/*get_c2disco_minimized: gets the minimized structure frm a Discover run  */
/*                   and copies the minimized structure to a molecule     */
/*                   num_molecule = gets the num_molecule molecule from   */
/*                                  the file (starts at 1)                */
/* cobbled from the get_disco routine August 98 DJW                       */
/**************************************************************************/

#include <stdio.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"

int get_c2disco_minimized(char *root_file,int mol_number,
						 atom *p_molecule, int num_atoms)

{
#include "header.h"
char fileroot[FILELEN_MAX];
FILE *fp;
int i;
int at_molecule = 0;
atom *p_atom;

strcpy(fileroot,root_file);
strcat(fileroot,".cor");

if (!(fp = fopen(fileroot, "r")))
    {
    fprintf(output_fp,"ZEBEDDE ERROR: error opening Discover file %s in get_disco_minimized\n", fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
    return(-1);
    }


else
	{
	/* ditch the top bit */
	fgets(dummy_head,MAX_LINE_LEN,fp);
        fgets(dummy_head,MAX_LINE_LEN,fp);
        fgets(dummy_head,MAX_LINE_LEN,fp);
        fgets(dummy_head,MAX_LINE_LEN,fp);

	i=0;

/******************************************************************/
/**** Modified to allow first molecule to be ignored DJW May 99 ***/
/******************************************************************/

		while(fgets(buffer,MAX_LINE_LEN,fp) != NULL )
	  	{
     		    if ( at_molecule == mol_number && strlen(buffer) > 10)
        		{
				p_atom = p_molecule + i;
        		sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf", &(p_atom->label[0]), 
	 				                        &p_atom->x,
            		                                        &p_atom->y, 
                                                                &p_atom->z,
                                                                &(p_atom->group[0]),
				                        	&(p_atom->group_no[0]),
            		                                        &(p_atom->pot[0]),
                                                                &(p_atom->elem[0]),
                                                                &p_atom->part_chge);
	

       			i++;
       			}
                     else if (strstr(buffer,"end") != NULL) 
			{
			at_molecule++;
			}
		}

		fclose(fp);
	}
return(1);
}

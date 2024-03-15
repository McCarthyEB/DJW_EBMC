/****************************************************************
******** Prints out a car file for further input   **************
******** Currently takes GULP input filename       **************
******** Started AJWL Jan 2008                     **************
****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
void print_biosym_header(FILE *file_fp, int pbc_flag);

void print_pbc_header(FILE *file_fp, double *p_abc);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void car_molecule_output(atom *p_pore, atom *p_molecule,
                                         int num_atoms, char *root_name)

{
#include "header.h"
char filename[FILELEN_MAX];
char tempno[FILELEN_MAX];
car_file_count++;
atom *p_start_atom;
strcpy(filename,root_name);
sprintf(tempno,"%d",car_file_count);
strcat(filename, tempno);
strcat(filename, ".car");
                                                                                                            
if (!(docked_template_fp = fopen(filename, "w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening CAR output file \n",
                      filename);
        return;
    }

/*print the headers*/

if (car_with_pore) print_biosym_header(docked_template_fp, pbc);
else print_biosym_header(docked_template_fp, FALSE);
fprintf(docked_template_fp,"Materials Studio generated CAR file\n");
fprintf(docked_template_fp,"!DATE Mon Oct 12 12:17:26 1998\n");

/* DEBUGfprintf(docked_template_fp, "%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n", abc[0], abc[1], abc[2], abc[3], abc[4], abc[5]);*/

/*print the molecule*/
if (car_with_pore)
    {
    if (pbc) print_pbc_header(docked_template_fp, &abc[0]);
	print_molecule( p_pore, num_pore_atoms, docked_template_fp, TRUE);
    fprintf(docked_template_fp,"end\n");
	p_start_atom= p_molecule;	
	print_molecule( p_molecule, num_atoms, docked_template_fp, TRUE);
	fprintf(docked_template_fp,"end\nend\n");
	}
else
	{
	print_molecule( p_molecule, num_atoms, docked_template_fp, FALSE);
	fprintf(docked_template_fp,"end\nend\n");
	}
strcpy(temp_car_output,filename);
fclose(docked_template_fp);
return;
}

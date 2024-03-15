/**************************************************************************/
/*anal_disco.c (!) : gets the relevant stuff from a Discover output file  */
/*                                                                        */
/*started DWl 3/5/95                                                      */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include <errno.h>
#include "global_values.h"
#include "structures.h"

#define NORMAL "Normal"
#define ABNORMAL "ABNORMAL"
#define START_MINIMIZE "Status Before"
#define END_MINIMIZE "Status After"
#define END "DISCOVER completed"
#define NONBOND "kcal=nonbond"
#define TOTAL "kcal=total"
#define DISCOVER_ERROR "ERROR"

int anal_disco(char *root_file, atom *p_molecule, int num_atoms, energy *p_energy)

{
#include "header.h"
char fileroot[FILELEN_MAX];

FILE *disc_output_fp;


int l_count, l_norm, l_abnorm;
int  l_start_minimize, l_end_minimize;
int status;
int end_found = FALSE;
int error_found = FALSE;
int delay_loop;



energy *p_answer;

p_answer = p_energy;


strcpy(fileroot,root_file);
strcat(fileroot,".out");
if (!(disc_output_fp = fopen(fileroot, "r")))
	{
      printf("trying to hang on a bit....\n");
      for (delay_loop=0;delay_loop<100000;delay_loop++)
         {
         /***** pause a while to clear buffers****/
         }
      printf("finished pausing...........\n");
      if (!(disc_output_fp = fopen(fileroot, "r")))
         {
         printf("...still no output file\n");
         fprintf(output_fp, "ZEBEDDE ERROR: error opening Discover output file %s in anal_disco\n",
                                         fileroot);
         fprintf(output_fp,"Errno = %d\n",errno);
         fprintf(output_fp,"Abandoning Minimisation\n");
         return(-1);
         }
      printf("....ok now\n");
	}
l_count = 0;
l_norm = 0;
l_abnorm = 0;
l_start_minimize = 0;
l_end_minimize = 0;

p_answer->minimizer_init_nonbond = -99999999.0;
p_answer->minimizer_init_total   = -99999999.0;

while (fgets(buffer, MAX_LINE_LEN, disc_output_fp) != NULL 
												&& end_found !=TRUE)
	{
	l_count++;
	if (strstr(buffer,NORMAL) != NULL) l_norm = l_count;
	else if (strstr(buffer,ABNORMAL) != NULL) l_abnorm = l_count;
	else if (strstr(buffer,START_MINIMIZE) != NULL) 
					l_start_minimize = l_count;
	else if (strstr(buffer,END_MINIMIZE)   != NULL) 
					l_end_minimize = l_count;
	else if (strstr(buffer,END) != NULL) end_found = TRUE;

	else if (strstr(buffer,DISCOVER_ERROR) != NULL) error_found = TRUE;

	else if (strstr(buffer,NONBOND) != NULL)	
		{
		sscanf(buffer,"%lf", &p_answer->minimizer_end_nonbond);
		if ((p_answer->minimizer_init_nonbond == -99999999.0) )
			{
			/*init_nonbond_energy = nonbond_energy; */
			p_answer->minimizer_init_nonbond=p_answer->minimizer_end_nonbond;
			}
		}

	else if (strstr(buffer,TOTAL) != NULL)	
		{
		sscanf(buffer,"%lf", &p_answer->minimizer_end_total);

		if ((p_answer->minimizer_init_total == -99999999.0) )
			{
			p_answer->minimizer_init_total=p_answer->minimizer_end_total;
			}
		}
			
	}

if (end_found == FALSE) 
	{
        fprintf(output_fp,"ZEBEDDE ERROR: Discover run %s did not terminate correctly\n",
            fileroot);
/***debug - fall out***/
/***exit(-1); ***/
	status = -1;
	}

if (error_found == TRUE) 
	{
        fprintf(output_fp,"ZEBEDDE ERROR: An error occured in run %s\n", fileroot);
/***debug - fall out***/
/***exit(-1);***/
	status = -1;
	}

if (l_norm > l_abnorm) status = TRUE; /* Normal > ABNORMAL */
else status = FALSE;

fclose(disc_output_fp);
return(status);
}


/**************************************************************************/
/*get_disco_minimized: gets the minimized structure frm a Discover run      */
/*                   and copies the minimized structure to a molecule     */
/*                   num_molecule = gets the num_molecule molecule from   */
/*                                  the file (starts at 1)                */
/*started DWl 5/5/95                                                      */
/**************************************************************************/
int get_disco_minimized(char *root_file,int mol_number,
						 atom *p_molecule, int num_atoms)

{
#include "header.h"
char fileroot[FILELEN_MAX];
FILE *fp;
int i;
int at_molecule = 1;
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

		while(fgets(buffer,MAX_LINE_LEN,fp) != NULL && 
										at_molecule == mol_number)
			{
        	if (strstr(buffer,"end") != NULL) 
				{
				at_molecule++;
				}
     		else
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
			}

		fclose(fp);
	}
return(1);
}

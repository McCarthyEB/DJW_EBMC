/****************************************************************************/
/*anal_c2disco.c (!) : gets the relevant stuff from a Discover output file  */
/*                                                                          */
/*started DWl 3/5/95                                                        */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include <errno.h>
#include "global_values.h"
#include "structures.h"

#define NORMAL "Normal"
#define ABNORMAL "ABNORMAL"
#define MINIMIZE_REPORT "Minimize:"
#define END_MINIMIZE "Minimize"
#define END "Exiting Discover"
#define NONBOND "vdW    "
#define TOTAL "Total Potential energy"
#define DISCOVER_ERROR "ERROR"

int anal_c2disco(char *root_file, atom *p_molecule, int num_atoms, energy *p_energy)

{
#include "header.h"
char fileroot[FILELEN_MAX];

FILE *disc_output_fp;

int status = 0;
int end_found = FALSE;
int error_found = TRUE;
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
         fprintf(output_fp, "ZEBEDDE ERROR: error opening Discover output file %s in anal_c2disco\n",
                                         fileroot);
         fprintf(output_fp,"Errno = %d\n",errno);
         fprintf(output_fp,"Abandoning Minimisation\n");
         return(-1);
         }
      printf("....ok now\n");
	}

p_answer->minimizer_init_nonbond = -99999999.0;
p_answer->minimizer_init_total   = -99999999.0;


while (fgets(buffer, MAX_LINE_LEN, disc_output_fp) != NULL && end_found !=TRUE)
	{
          if (strstr(buffer,NONBOND) != NULL)
            {
              sscanf(buffer,"%*s%lf%lf", &(p_answer->minimizer_init_nonbond),
                                         &(p_answer->minimizer_end_nonbond));
            }
 
	  else if (strstr(buffer,TOTAL) != NULL)	
            {
	      sscanf(buffer,"%*s%*s%*s%lf%lf", &(p_answer->minimizer_init_total),
                                                 &(p_answer->minimizer_end_total));
              error_found = FALSE;
            }
          else if (strstr(buffer,MINIMIZE_REPORT))
            {
              fprintf(output_fp, "%s\n", buffer);
            }
          else if (strstr(buffer, END) != NULL)
            {
              end_found= TRUE;
            }
	}

if (end_found == FALSE) 
	{
        fprintf(output_fp,"ZEBEDDE ERROR: Discover run %s did not terminate correctly\n",
            fileroot);
	status = -1;
	}

if (error_found == TRUE) 
	{
        fprintf(output_fp,"ZEBEDDE ERROR: No total energy reported after c2discover run %s\n", fileroot);
	status = -1;
	}

fclose(disc_output_fp);
return(status);
}



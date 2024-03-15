/**************************************************************************/
/*anal_gulp.c : gets the relevant stuff from a gulp output file           */
/*                                                                        */
/*started KEJ 14/11/05                                                    */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

#define GULP_ERROR "ERROR"

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );
int int_from_string(char *p_string, int *p_int);

int anal_gulp(char *root_file, atom *p_molecule, int num_atoms, 
                                                    energy *p_energy)

{
#include "header.h"
char fileroot[FILELEN_MAX];

FILE *gulp_output_fp;


int status = 0;
int end_found;
int error_found;
char string4[80];

int found_energy;
int single_run, opti_run;

energy *p_answer;

p_answer = p_energy;

end_found = FALSE;
error_found = FALSE;



/**** start of main routine ****/

strcpy(fileroot,root_file);
strcat(fileroot,".gout");

if (!(gulp_output_fp = fopen(fileroot, "r")))
   		{
                fprintf(output_fp,"ZEBEDDE WARNING: ERROR opening gulp output file %s in anal_gulp\n",
                                                        fileroot);
                fprintf(output_fp,"Probable cause : Error in gulp input\n");
                if (DEBUG) fprintf(stdout,"Probable cause : Error in gulp input\n");
                fprintf(output_fp,"Abandoning Minimisation\n");
                fprintf(output_fp," Please Report to Authors\n");

   		fclose(gulp_output_fp);
  		return(-1);
   		}
else
		{
		found_energy=FALSE;

		
		opti_run = FALSE;
		single_run = FALSE;

		while (!feof(gulp_output_fp)){
		    fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
		   
		    /***check for error message on any line of output***/
		    if (strstr(buffer, "ERROR") != 0)
			{
			 fprintf(output_fp, "Warning gulp output contains an error\n");
			 fprintf(output_fp, "Probable cause: Error in gulp input\n");
			 fprintf(output_fp, "Abandoning minimisation\n");	
			 return(-1);
			}
		    
		    /* is it an optimisation or a single calc? */
		    if(strstr(buffer, " opti ") != 0)
			{
			opti_run = TRUE;
			printf("DB> gulp is opti\n");
			}
		    else if(strstr(buffer, " single ") != 0)
			{
			single_run = TRUE;
		  	printf("DB> gulp is single\n");
			}


		    /*** Optimisation run: find the final energy from output ***/
		    if(opti_run)
		    {
		    if (strstr(buffer, "Final energy") != 0)
			{
			 printf("Line with final energy found:\n %s\n", buffer);
			 sscanf(buffer, "%*s %*s %*s %s", string4);
			 
			 /* convert string to double number */		   	 
			 p_answer->minimizer_end_total = atof (string4);
			 printf("The final energy taken from gulp output is %12.8f\n", p_answer->minimizer_end_total);
			 found_energy = TRUE;
			}
		    }
		    /*** Single run: find the final energy from output ***/
                    if(single_run)
		    {
		    if ((strstr(buffer, "Total lattice energy") != 0) && found_energy != TRUE)
                        {
                         printf("Line with final energy found:\n %s\n", buffer);
                         sscanf(buffer, "%*s %*s %*s %*s %s", string4);
                         
			 /* convert string to double number */                                                         
                         p_answer->minimizer_end_total = atof (string4);
                         printf("The total lattice energy taken from gulp output is %12.8f\n",
					 p_answer->minimizer_end_total);
                         found_energy = TRUE;
                        }
		    }
		   
		   }
		}

	    /*** check that an energy has been found and assign relevant status***/
	    if (found_energy == TRUE)
		{
		status = 1;
		}
	    else
		{
		 status=0;
		}

fclose(gulp_output_fp);

return(status);
}

/**************************************************************************/
/*get_gulp_minimized : gets the minimized structure from a GULP   run     */
/*                   and copies the minimized structure to a molecule     */
/*                                                                        */
/*started KEJ 14/11/05                                                    */
/**************************************************************************/
int get_gulp_minimized(char *root_file,int num_skip_atoms,
							 atom *p_molecule, int num_atoms)
 
{
#include "header.h"
#include <stdlib.h>
char fileroot[FILELEN_MAX];
FILE *gulp_output_fp;
int i;
char *found_molecule;
atom *p_atom;

int found_coords;
char x_string[80], y_string[80], z_string[80];
int single_run, opti_run;



/** main part of routine **/

strcpy(fileroot,root_file);
strcat(fileroot,".gout");
                                                                                                                       
 
if (!(gulp_output_fp = fopen(fileroot, "r")))
    {
    	fprintf(output_fp, "ZEBEDDE WARNING: Error opening GULP file %s in get_minimized\n", fileroot);
    	fprintf(output_fp, "Abandoning Minimisation\n");

    	return(-1);
    }
 
 
else
    {
	opti_run = FALSE;
        single_run = FALSE;
                                                                                                  
        while (opti_run == FALSE && single_run == FALSE)
         { 
	        fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
                                                                                                  
                    /* is it an optimisation or a single calc? */
                    if(strstr(buffer, "opti") != 0)
                        {
                        opti_run = TRUE;
                        }
                    else if(strstr(buffer, "single") != 0)
                        {
                        single_run = TRUE;
                        }
                                                                                                  
	 }

      if(opti_run)
	{
	/***read until we get the final structure **/
	printf("Reading GULP minimized structure\n");
	
	found_coords =FALSE;

	while (!feof(gulp_output_fp) && found_coords == FALSE)
	{
	fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
	if(strstr(buffer, "Final fractional coordinates of atoms") != 0)
	  {
	  found_coords = TRUE;
	  fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
	  fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
	  fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
	  fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
	  fgets(buffer, MAX_LINE_LEN, gulp_output_fp);

	   /** copy these lines (label,x,y,z) into atom  - until reach end **/
	   for (i=0; i<=10000;i++)
			{
			p_atom = p_molecule+i;
			fgets(buffer, MAX_LINE_LEN, gulp_output_fp);
			if (strstr(buffer, "-----------") == 0)
			   {
			    sscanf(buffer,"%*s %s %*s %s %s %s", &(p_atom->label[0]), x_string, y_string, z_string);
			    p_atom->x = atof (x_string);
			    p_atom->y = atof (y_string);
			    p_atom->z = atof (z_string);
				
			    printf("%s %12.8f %12.8f %12.8f\n", &(p_atom->label[0]), p_atom->x,p_atom->y,p_atom->z);
			   }
			else break;
			}
          }
	
	}
	}
      else if (single_run)
	{
	printf("No need to get final structure - as only a single point calc\n");
	}
	fclose(gulp_output_fp);
    }
return(1);	
}

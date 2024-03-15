/**************************************************************************/
/*anal_mopac.c (!) : gets the relevant stuff from a mopac output file  */
/*                                                                        */
/*started DWl 3/5/95                                                      */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

#define CYCLE "CYCLE"
#define HEAT "HEAT"
#define   NORMAL "GRADIENT TEST PASSED"
#define ABNORMAL "GRADIENT TEST NOT PASSED"
#define TOTAL "TOTAL ENERGY"
#define MIN_GEOM "FINAL GEOMETRY OBTAINED"
#define MOPAC_ERROR "ERROR"
#define TOO_BIG "MAX. NUMBER OF ATOMS"
#define SATISFIED "SATISFIED"
#define SUMMARY "SUMMARY"
#define ALL_READY_MINIMIZED "GRADIENTS WERE INITIALLY ACCEPTABLY SMALL"

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );
int int_from_string(char *p_string, int *p_int);

int anal_mopac(char *root_file, atom *p_molecule, int num_atoms, 
                                                    energy *p_energy)

{
#include "header.h"
char fileroot[FILELEN_MAX];

FILE *mopac_output_fp;


int status = 0;
int end_found;
int error_found;
int buffer_len;
int int_buffer[256];
int flag_int;
int dummy_int;
double temp_energy;

char  *hpos;

energy *p_answer;

p_answer = p_energy;

end_found = FALSE;
error_found = FALSE;

strcpy(fileroot,root_file);
strcat(fileroot,".out");

if (!(mopac_output_fp = fopen(fileroot, "r")))
   		{
                fprintf(output_fp,"ZEBEDDE WARNING: ERROR opening mopac output file %s in anal_mopac\n",
                                                        fileroot);
                fprintf(output_fp,"Probable cause : Error in Mopac input\n");
                if (DEBUG) fprintf(stdout,"Probable cause : Error in Mopac input\n");
                fprintf(output_fp,"Abandoning Minimisation\n");
                fprintf(output_fp," Please Report to Authors\n");

   		fclose(mopac_output_fp);
  		return(-1);
   		}
else
		{
		/***check for error messages on third line***/
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		if (strstr(buffer, TOO_BIG) != NULL)
   			{
                        max_mopac_atoms= atoi(strstr(buffer, ":")+1);
                        fprintf(output_fp,"SERIOUS MOPAC Error: Maximum number of atoms exceeded\n");
                        if (DEBUG) fprintf(stdout,"SERIOUS MOPAC Error: Maximum number of atoms exceeded\n");
                        fprintf(output_fp,"Max number allowed = %i\n",max_mopac_atoms);
                        fprintf(output_fp,"MOPAC must be recompiled if this is to be solved\n");
                        fclose(mopac_output_fp);
                        return(-2);

   			}
		}
fclose(mopac_output_fp);

/*****************read in the log file for the answers*************/
strcpy(fileroot,root_file);
strcat(fileroot,".log");

/********initialise energies***********************/
p_answer->minimizer_init_nonbond = -99999999.0;
p_answer->minimizer_init_total   = -99999999.0;
p_answer->minimizer_end_nonbond = +99999999.0;
p_answer->minimizer_end_total   = +99999999.0;

/****get the start and end heat of formatio****/

/***skip the top***/

if (!(mopac_output_fp = fopen(fileroot, "r")))
    {
    if (DEBUG) fprintf(stdout,"ZEBEDDE WARNING: ERROR opening mopac log file %s in anal_mopac\n",
                                                   fileroot);
    fprintf(output_fp,"ZEBEDDE WARNING: ERROR opening mopac log file %s in anal_mopac\n",
                                                   fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
    fclose(mopac_output_fp);
    return(-1);
    }
else
	{
	/*****skip the top 3 line******/
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);

	/****now read in each cycle line until no more****/
	/****get the initial heat****/

	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	hpos=strstr(buffer,HEAT);
	if (hpos !=NULL)
		{
		p_answer->minimizer_init_total = atof(hpos+5);	
                if (DEBUG) 
		    printf("Initial Energy: %f\n",p_answer->minimizer_init_total);
		/****find the last one***/
		while (hpos !=NULL)
			{
			fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
			hpos=strstr(buffer,HEAT);
			if (hpos !=NULL) p_answer->minimizer_end_total = atof(hpos+5);
			}

		/***get the heat from this line***/
                if (DEBUG) 
		    printf("Final Energy: %f\n",p_answer->minimizer_end_total);

	    /******check the final status of the run *****/
	    while(strstr(buffer, MIN_GEOM) == NULL)
		    {
		    fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		    if (strstr(buffer,ABNORMAL) != NULL) status =0;
		    else if (strstr(buffer,NORMAL) != NULL) status =1;
		    }
		}

	else
		/****Either and error or no need to minimize!*******/
		{
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		if (strstr(buffer,SUMMARY) == NULL)
			{
                        fprintf(output_fp,"ZEBEDDE WARNING: SERIOUS MOPAC Error: No results in logfile\n");
                        if (DEBUG) fprintf(stdout,"ZEBEDDE WARNING: SERIOUS MOPAC Error: No results in logfile\n");
                        fprintf(output_fp,"Abandoning Minimisation\n");
                        fprintf(output_fp,"Please Report to Authors\n");

			fclose(mopac_output_fp);
   			return(-1);
			}
		else
			{
			 /* printf("No minimization required: Gradient test already satisfied\n");*/
			while(strstr(buffer, MIN_GEOM) == NULL)
            	{
            	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
            	if (strstr(buffer,ALL_READY_MINIMIZED) != NULL) status =1;
            	else if (strstr(buffer,HEAT) != NULL) 
					{
					/****get the heat (init=final) from this line******/
					buffer_len=int_from_string(&buffer[0], &int_buffer[0]);

					temp_energy = get_doub(&int_buffer[0], 
								buffer_len, &dummy_int, &flag_int);
					if (!flag_int)
						{
                                                 if (DEBUG) fprintf(stdout,"Serious MOPAC Error: No Heat of formation\n");
                                                 fprintf(output_fp,"Serious MOPAC Error: No Heat of formation\n");
                                                 fprintf(output_fp,"Abandoning Minimisation\n");
                                                 fprintf(output_fp,"Please Report to Authors\n");
            			                 fclose(mopac_output_fp);
            			                 return(-1);
						}
					p_answer->minimizer_end_total = temp_energy;
					p_answer->minimizer_init_total =temp_energy;
					status = 1;
					}
            	}
			}

		}
		

	if (DEBUG) printf("MOPAC Molecule Minimization Status %i\n",status);
	if (status ==  0)
		{
		/***check to see if the energy decreased***/
		if (p_answer->minimizer_end_total < p_answer->minimizer_init_total)
			{
			if (DEBUG) printf("Energy was lowered so accepting move\n");
			status =1;
			}
		}
	}
fclose(mopac_output_fp);

if (p_answer->minimizer_end_total == 99999999) exit(1);

return(status);
}

/**************************************************************************/
/*get_mopac_minimized : gets the minimized structure from a MOPAC run     */
/*                   and copies the minimized structure to a molecule     */
/*                   num_skip_atoms = num atoms to skip before            */
/*                   starting to read				  */
/*                   to allow skipping of the pore                        */
/*                                                                        */
/*started DWl 20/6/95                                                     */
/**************************************************************************/
int get_mopac_minimized(char *root_file,int num_skip_atoms,
							 atom *p_molecule, int num_atoms)
 
{
#include "header.h"
char fileroot[FILELEN_MAX];
FILE *mopac_output_fp;
int i;
atom *p_atom;
 
strcpy(fileroot,root_file);
strcat(fileroot,".log");

 
if (!(mopac_output_fp = fopen(fileroot, "r")))
    {
    fprintf(output_fp, "ZEBEDDE WARNING: Error opening MOPAC file %s in get_minimized\n", fileroot);
    fprintf(output_fp, "Abandoning Minimisation\n");

    return(-1);
    }
 
 
else
    {

	/***read until we get the final srtructure **/
	/*printf("Reading MOPAC minmized structure\n");*/
	/*printf("skipping to start of molecule...\n");*/

	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	while(strstr(buffer,MIN_GEOM) == NULL)
		{
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		}

	/***skip to new geometry***/
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
	
	/***skip num_atoms***/
	for (i=0;i<num_skip_atoms;i++)
		{
		fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
		}
	for (i=0; i<=num_atoms;i++)
			{
			p_atom = p_molecule+i;
			fgets(buffer, MAX_LINE_LEN, mopac_output_fp);
			sscanf(buffer,"%*s%lf%*i%lf%*i%lf", &p_atom->x, 
										&p_atom->y,&p_atom->z);
			}
	fclose(mopac_output_fp);
	}

return(1);	
}

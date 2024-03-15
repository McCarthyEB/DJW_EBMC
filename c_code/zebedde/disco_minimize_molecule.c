/***************************************************************************/
/* disco_minimize_molecule.c : minimizes a molecule with respect to itself */
/*                       currently writes out files for DISCOVER           */
/* started 2/5/95 DWL                                                      */
/*                                                                         */
/* 12/6 SERIOUS MODS! DWL                                                  */
/* Now works with Discover 94                                              */
/*                                                                         */
/* August 98 removed call to calculate energy so that this routine does not*/
/* need all the gubbins for the Ewald sum DJW                              */
/* Adapted for multi-molecules March 2013 DJW                              */
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void regroup(atom *p_molecule, int num_atoms);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void print_biosym_header(FILE *file_fp,int pbc_flag);

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void car_end(FILE *fp);

void delete_file(char *fileroot, char *extension);

int anal_disco(char *root_file, atom *p_molecule,int num_atoms,energy *p_energy);

int  get_disco_minimized(char *root_file,int num_molecule, 
						atom *p_molecule, int num_atoms);

void print_mdf(char *p_molname, atom *p_molecule, int num_atoms,
                    int start, FILE *fp, int D_to_H);

void print_mdf_header(char *p_assembly, FILE *fp);

void print_mdf_end(FILE *fp);

void print_neighbours( atom *p_molecule, int num_atoms,  FILE *fp);

void disco_minimize_molecule( atom *p_molecule, int num_atoms, 
                              char *root_name, int which_mol )
{
#include "header.h"

char fileroot[FILELEN_MAX];
char disc_shell[FILELEN_MAX];
FILE *disc_shell_fp;
char *fullstop;
char discover_cmdline[256];
int discover_status,get_status;
int D_to_H;

/****DEBUG*****/
/* printf("Entering disco_minimise_molecule with %d atoms\n", num_atoms);*/
/****ENDDEBUG*****/

/*************print out discover .car file***********************************/

strcpy(fileroot,root_name);
strcat(fileroot,".car");
if (!(discover_fp = fopen(fileroot, "w")))
    {
      fprintf(output_fp, "ZEBEDDE ERROR: Error opening Discover car file %s in disco_minimize_molecule\n", fileroot);
      fprintf(output_fp, "Abandoning Minimisation\n");
      return;
    }

D_to_H= TRUE;
print_biosym_header(discover_fp,0); /***no pbc in molecule minimization**/
print_frame_header("molecule_minimize",1, discover_fp);
print_molecule(p_molecule, num_atoms, discover_fp, D_to_H);
car_end(discover_fp);

/*************print out discover .mdf file***********************************/

strcpy(fileroot,root_name);
strcat(fileroot,".mdf");

if (!(discover_fp = fopen(fileroot, "w")))
    {
    fprintf(output_fp,
         "ZEBEDDE ERROR: Error opening Discover mdf file %s in disco_minimize_molecule\n", fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
    return;
    }

if (DEBUG)
   {
    print_neighbours( p_molecule, num_atoms, stdout); 
   }

print_mdf_header(NULL,discover_fp);
print_mdf("molecule_minimize", p_molecule, num_atoms, 0, discover_fp, D_to_H);
print_mdf_end(discover_fp);

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

/******NOW RUN IN DISCOVER AND RETRIEVE OUTPUT *****/
/***** call Discover script to runjob now at full on priority ****/
/****** syntax: discover file_name discover_forcefield_name nice start_yn      ****/
/****** Not anymore! Latest discover uses just discover file_name - AJWL 9/08 ******/

sprintf(disc_shell,"./%s.csh",fileroot);

if (!(disc_shell_fp = fopen(disc_shell,"r")))
	{
	/* .disc_shell does not exist */
	/* so generate but don't run : BUG 4/7/96 DWL, as it detaches from shell */
    /* but first delete any old files just in case we get prompted           */
	delete_file(fileroot,".rst");
	delete_file(fileroot,".prm");
//	delete_file(fileroot,".out");
	delete_file(fileroot,".cor");
	sprintf(discover_cmdline,"%s %s", discover_path,fileroot);
	system (discover_cmdline);
	}

/*disc_shell exists so just run it*/
//sprintf(discover_cmdline,"%s>/dev/null",disc_shell);
//fclose(disc_shell_fp);

/* strcat(discover_cmdline,">/dev/null\n"); */
/***new version ****/
/* sprintf(discover_cmdline,"disco3 %s",fileroot); */


//system (discover_cmdline);	

/*******get results from DISCOVER run *************************/

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';
discover_status = anal_disco(fileroot, p_molecule,
		   	num_atoms, &interaction_energy[which_mol]);
/****DEBUG*****/
/* printf("Returned status %i  ",discover_status); */
/* if (discover_status == TRUE) */
	/* { */
	/* printf(" VALID Minimization\n"); */
	/* } */
/* else */
	/* { */
	/* printf(" INVALID Minimization\n"); */
	/* } */
/****DEBUG*****/




/******* Was it a good move in DISCOVER ***********************/

/**** since we are minimising a molecule, test on the total energy *****/
/**** accept as long as energy lowered (wether converged or not   *****/
if (discover_status >= 0)
    {
	if (interaction_energy[which_mol].minimizer_end_total 
			< interaction_energy[which_mol].minimizer_init_total)
		{

		get_status = get_disco_minimized(fileroot,1, p_molecule, num_atoms);
		if (get_status ==-1)
			{
			fprintf(output_fp, "ZEBEDDE WARNING: Error reading in minimized structure\n");
			fprintf(output_fp, "ZEBEDDE WARNING: in disco_minimize_molecule\n");
                        fprintf(output_fp, "Abandoning Minimization\n");
                        discover_status = FALSE;
			}
		}
	else
		{
		discover_status = FALSE;
		}
	}
if (discover_status == -1)
    {
    /*****error or unfinished run***/
    fprintf(output_fp,"ZEBEDDE WARNING: Ignoring Discover Run - Error in Run\n");
    fprintf(output_fp,"ZEBEDDE WARNING: If problem persists report to Authors\n");
    }

/***tidy up a bit to minimize disk usage***/

delete_file(fileroot,".rst");
delete_file(fileroot,".prm");
//delete_file(fileroot,".out");
delete_file(fileroot,".cor");


return;
}

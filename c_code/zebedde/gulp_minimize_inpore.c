/*********************************************************************/
/* gulp_minimize_inpore.c : minimizes a molecule with respect        */	
/*                       to (fixed) pore using GULP                  */
/* started KEJ 09/11/05                                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

#define GULP_SUCCESS 1

#define AL_CHARGE 3.00
#define SI_CHARGE 4.00
#define OCORE_CHARGE 0.86902
#define OSHEL_CHARGE -2.86902
#define N_CHARGE -0.353
#define H_CHARGE 0.153893
#define C1_C2_CHARGE -0.178
#define C3_CHARGE -0.383
#define NZ_CHARGE 0.000
#define F_CHARGE 0.000
#define CL_CHARGE 0.000
#define AR_CHARGE 0.000
#define OO_CHARGE 0.000
#define C_CHARGE 0.000
#define OH_CHARGE 0.000
#define HO_CHARGE 0.000
#define CP_CHARGE 0.000
#define HS_CHARGE 0.000
#define SH_CHARGE 0.000


void delete_file(char *fileroot, char *extension);

int get_gulp_minimized(char *root_file,int num_skip_atoms,
                                                atom *p_molecule, int num_atoms);

void minimize_inpore(atom *p_pore, int num_pore_atoms,
                               atom *p_molecule, int num_atoms, char *root_name);

int anal_gulp(char *root_file, atom *p_molecule, int num_atoms,
                                                    energy *p_energy);

int calculate_charge(atom *p_molecule, int num_atoms);

void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name);

void gulp_minimize_inpore(atom *p_pore, atom *p_molecule, 
                                         int num_atoms, char *root_name)


{
#include "header.h"

char fileroot[FILELEN_MAX];
char gulp_cmd[FILELEN_MAX];
char call_gulp[FILELEN_MAX];
char catenate_pots[FILELEN_MAX];
FILE *gulp_fp;
FILE *pots_fp;
char *fullstop;
int gulp_status;
atom *p_atom;


int i;
int charge;
char dummy[80];
char line_pots[BUFFER];
char extraline[BUFFER], oldline[BUFFER];
char extrafile[BUFFER];

char *p_Al_elem;
char *p_Si_elem;
char *p_O_elem;
char *p_N_elem;
char *p_H_elem;
char *p_C1_elem;
char *p_C2_elem;
char *p_C3_elem;
char *p_NZ_elem;
char *p_F_elem;
char *p_CL_elem;
char *p_AR_elem;
char *p_OO_elem;
char *p_C_elem;
char *p_OH_elem;
char *p_HO_elem;
char *p_CP_elem;
char *p_HS_elem;
char *p_SH_elem;



p_Al_elem= "Al";
p_Si_elem= "Si";
p_O_elem= "O";
p_N_elem= "N";
p_H_elem= "H";
p_C1_elem= "C1";
p_C2_elem= "C2";
p_C3_elem= "C3";
p_NZ_elem= "nz";
p_F_elem= "f";
p_CL_elem= "cl";
p_AR_elem= "ar";
p_OO_elem= "o=";
p_C_elem= "c";
p_OH_elem= "oh";
p_HO_elem= "ho";
p_CP_elem= "cp";
p_HS_elem= "hs";
p_SH_elem= "sh";

/*************print out GULP input file***********************************/

printf("in gulp routine\n");

strcpy(fileroot,root_name);
strcat(fileroot,".gin");

if (!(gulp_fp = fopen(fileroot, "w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening GULP input file %s in minimize_molecule\n",
                      fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
	return;
    }

/** print the command lines **/
fprintf(gulp_fp, "# Keywords:\n");
fprintf(gulp_fp,"%s\n", gulp_cmdline_inpore);

/** print out any extra command lines **/
strcpy(extrafile, "extra_gulp.lines");
extra_gulp_fp = fopen(extrafile, "r");


while (!feof(extra_gulp_fp))
   {
    fgets(extraline,FILELEN_MAX,extra_gulp_fp);
    if( strcmp(extraline, oldline) != 0)
	{
	fprintf(gulp_fp,"%s", extraline);
    	}
	strcpy(oldline, extraline);
   }


fclose(extra_gulp_fp);


/** print the cell parameters **/
fprintf(gulp_fp, "cell\n");
fprintf(gulp_fp, "%12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n", abc[0], abc[1], abc[2], abc[3], abc[4], abc[5]);

/** print the start of coordinates **/
fprintf(gulp_fp, "cartesian %d\n", (num_pore_atoms + num_atoms));



/***write the pore fixed geometry  - take charges from above***/
for (i=0;i<=num_pore_atoms;i++)
	{
	 p_atom = p_pore +i;
	 if (strncmp(p_Al_elem, (p_atom->label), 2) == 0)
	 {
	 fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 0 0 0\n", p_atom->label,
									p_atom->x, p_atom->y, p_atom->z, AL_CHARGE);
	 }
 	 if (strncmp(p_Si_elem, (p_atom->label), 2) == 0)
         {
         fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 0 0 0\n", p_atom->label, 
                                                                        p_atom->x, p_atom->y, p_atom->z, SI_CHARGE);
         }
         if (strncmp(p_O_elem, (p_atom->label), 1) == 0)
         {
         fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 0 0 0\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, OCORE_CHARGE);
	 fprintf(gulp_fp,"%s shel %12.8f %12.8f %12.8f %12.8f 0 0 0\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, OSHEL_CHARGE);

         }

	}
/***write the template for optimisation  - take charges from above ***/
for (i=0;i<=num_atoms;i++)
	{
	p_atom = p_molecule +i;
	if (strncmp(p_N_elem, (p_atom->label), 1) == 0)
         {
         fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, N_CHARGE);
	 }
	if (strncmp(p_H_elem, (p_atom->label), 1) == 0)
         {
	 fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, H_CHARGE);
 	 }
	if (strncmp(p_C1_elem, (p_atom->label), 2) == 0 || strncmp(p_C2_elem, (p_atom->label), 2) == 0 )
         {
         fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, C1_C2_CHARGE);
	 }
	if (strncmp(p_C3_elem, (p_atom->label), 2) == 0)
         {
	 fprintf(gulp_fp,"%s core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", p_atom->label,
                                                                        p_atom->x, p_atom->y, p_atom->z, C3_CHARGE);
	 }
	if (strcmp(p_NZ_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"N core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",
                                                                        p_atom->x, p_atom->y, p_atom->z, NZ_CHARGE);
         }

        if (strcmp(p_F_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"F core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",
                                                                        p_atom->x, p_atom->y, p_atom->z, F_CHARGE);
         }

        if (strcmp(p_CL_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"Cl core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",
                                                                       p_atom->x, p_atom->y, p_atom->z, CL_CHARGE);
         }

        if (strcmp(p_AR_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"Ar core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, AR_CHARGE);
         }

        if (strcmp(p_OO_elem, (p_atom->pot)) == 0)
        { 
         fprintf(gulp_fp,"O core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, OO_CHARGE);                  
        }

        if (strcmp(p_C_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"C core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, C_CHARGE);                   
        }

	if (strcmp(p_OH_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"O core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, OH_CHARGE);
        }

        if (strcmp(p_HO_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"H core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, HO_CHARGE);
        }

        if (strcmp(p_CP_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"C core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, CP_CHARGE);
        }

        if (strcmp(p_HS_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"H core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, HS_CHARGE);
        }

        if (strcmp(p_SH_elem, (p_atom->pot)) == 0)
        {
        fprintf(gulp_fp,"S core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, SH_CHARGE);
        }

	}

/** add  potentials to end of file **/

gulp_pots_fp = fopen(gulp_pots_file, "r");

while (!feof(gulp_pots_fp))
   { 
    fgets(line_pots,FILELEN_MAX,gulp_pots_fp);
    fprintf(gulp_fp,line_pots);
   }

fclose(gulp_fp);
fclose(gulp_pots_fp);

/******Run In GULP**************************/

sprintf(call_gulp,"%s <%s> %s.gout",gulp_path,fileroot, root_name);


printf("Going to gulp\n");
system (call_gulp);	
printf("Back from gulp\n");



/*******get results from gulp run *************************/

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

// gulp_status = anal_gulp(fileroot,p_molecule,num_atoms, &interaction_energy);

printf("DB>> back from gulp, gulp status %d\n", gulp_status);

/** if gulp run was successful, then copy the new structure using get_gulp_minimized **/
if (gulp_status == GULP_SUCCESS)
	{
	get_gulp_minimized(fileroot, num_pore_atoms, p_molecule, num_atoms);
        
	}
else
	{
        fprintf(output_fp,"Invalid Result from GULP\n");
	exit(1);
	}

return;
}

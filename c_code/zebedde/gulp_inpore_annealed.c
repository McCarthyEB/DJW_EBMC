/*********************************************************************/
/* gulp_inpore_annealed.c : just produces a GULP input file for     */
/*			template and pore,                           */
/*                       does not carry out minimisation             */
/*         KEJ 29/11/06                                              */
/*  The gulp input file will be labeled _annealed_no.gin,           */
/*  and is set up                                                    */
/*  to perform minimisation with fixed pore (other commands as       */
/*  supplied in input)                                               */
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
#define OC_CHARGE 0.000
void delete_file(char *fileroot, char *extension);

int get_gulp_minimized(char *root_file,int num_skip_atoms,
                                                atom *p_molecule, int num_atoms);

void minimize_inpore(atom *p_pore, int num_pore_atoms,
                               atom *p_molecule, int num_atoms, char *root_name);

int anal_gulp(char *root_file, atom *p_molecule, int num_atoms,
                                                    energy *p_energy);

int calculate_charge(atom *p_molecule, int num_atoms);

void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name);

void gulp_inpore_annealed(atom *p_pore, atom *p_molecule, 
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
char line_pots[BUFFER];
char extraline[BUFFER], oldline[BUFFER];
char extrafile[BUFFER];
char tempno[80];

char *p_Al_elem;
char *p_Si_elem;
char *p_O_elem;
char *p_N_elem;
char *p_N4_elem;
char *p_H_elem;
char *p_HC_elem;
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
char *p_OC_elem;
char *p_NA_elem;



p_Al_elem= "al";
p_Si_elem= "sz";
p_O_elem= "oss";
p_N_elem= "n";
p_N4_elem = "n4";
p_H_elem= "h";
p_HC_elem= "hc";
p_C1_elem= "c1";
p_C2_elem= "c2";
p_C3_elem= "c3";
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
p_OC_elem= "oc";
p_NA_elem= "na";

/*************print out GULP input file***********************************/

printf("in gulp_inpore_annealed routine-just producing input- no minimisation!!\n");

strcpy(fileroot, root_name);
strcat(fileroot,"_annealed");

sprintf(tempno, "%d", no_annealed_temps);

printf("hopefully no: %d has been converted to string: %s\n", no_annealed_temps, tempno);

strcat(fileroot,tempno);

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
	 if (strcmp(p_Al_elem, (p_atom->pot)) == 0)
	 {
	 fprintf(gulp_fp,"Al core %12.8f %12.8f %12.8f %12.8f 0 0 0\n",
									p_atom->x, p_atom->y, p_atom->z, AL_CHARGE);
	 }
 	 if (strcmp(p_Si_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"Si core %12.8f %12.8f %12.8f %12.8f 0 0 0\n", 
                                                                        p_atom->x, p_atom->y, p_atom->z, SI_CHARGE);
         }
         if (strcmp(p_O_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"O core %12.8f %12.8f %12.8f %12.8f 0 0 0\n", 
                                                                        p_atom->x, p_atom->y, p_atom->z, OCORE_CHARGE);
	 fprintf(gulp_fp,"O shel %12.8f %12.8f %12.8f %12.8f 0 0 0\n", 
                                                                        p_atom->x, p_atom->y, p_atom->z, OSHEL_CHARGE);
	 if(strcmp(p_O_elem, (p_atom->pot)) != 0 && strcmp(p_Al_elem, (p_atom->pot)) != 0 && strcmp(p_Si_elem, (p_atom->pot)) != 0 && strcmp(p_N_elem, (p_atom->pot)) != 0 && strcmp(p_H_elem, (p_atom->pot)) != 0 && strcmp(p_C1_elem, (p_atom->pot)) != 0 && strcmp(p_C2_elem, (p_atom->pot)) != 0 && strcmp(p_C3_elem, (p_atom->pot)) != 0 && strcmp(p_NZ_elem, (p_atom->pot)) != 0 && strcmp(p_F_elem, (p_atom->pot)) != 0 && strcmp(p_CL_elem, (p_atom->pot)) != 0 && strcmp(p_AR_elem, (p_atom->pot)) != 0 && strcmp(p_OO_elem, (p_atom->pot)) != 0 && strcmp(p_C_elem, (p_atom->pot)) != 0 && strcmp(p_OH_elem, (p_atom->pot)) != 0 && strcmp(p_HO_elem, (p_atom->pot)) != 0 && strcmp(p_CP_elem, (p_atom->pot)) != 0 && strcmp(p_SH_elem, (p_atom->pot)) != 0 && strcmp(p_HS_elem, (p_atom->pot)) != 0)
	{
	printf("Error - a gulp file can not be created as there are incompatible atom/pots types!\n");
	fprintf(gulp_fp, "ERROR - a gulp file can not be created as there are incompatible atom/pots types!\n");
	}


         }

	}
/***write the template for optimisation  - take charges from above ***/
for (i=0;i<=num_atoms;i++)
	{
	p_atom = p_molecule +i;
	if (strcmp(p_N_elem, (p_atom->pot)) == 0 || strcmp(p_N4_elem, (p_atom->pot)) == 0)
         {
         fprintf(gulp_fp,"N core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", 
                                                                        p_atom->x, p_atom->y, p_atom->z, N_CHARGE);
	 }
	if (strcmp(p_H_elem, (p_atom->pot)) == 0 || strcmp(p_HC_elem, (p_atom->pot)) == 0)
         {
	 fprintf(gulp_fp,"H core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", 
                                                                        p_atom->x, p_atom->y, p_atom->z, H_CHARGE);
 	 }
	if (strcmp(p_C1_elem, (p_atom->pot)) == 0 )
         {
         fprintf(gulp_fp,"C1 core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",
                                                                        p_atom->x, p_atom->y, p_atom->z, C1_C2_CHARGE);
	 }

	if (strcmp(p_C2_elem, (p_atom->pot)) == 0  )
         {
         fprintf(gulp_fp,"C2 core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",
                                                                        p_atom->x, p_atom->y, p_atom->z, C1_C2_CHARGE);
	 }

	if (strcmp(p_C3_elem, (p_atom->pot)) == 0)
         {
	 fprintf(gulp_fp,"C3 core %12.8f %12.8f %12.8f %12.8f 1 1 1\n", 
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

	if (strcmp(p_OC_elem, (p_atom->pot)) == 0)        {
        fprintf(gulp_fp,"O core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, OC_CHARGE);
        }

        if (strcmp(p_NA_elem, (p_atom->pot)) == 0)        {
        fprintf(gulp_fp,"N core %12.8f %12.8f %12.8f %12.8f 1 1 1\n",                                                                                                                                              p_atom->x, p_atom->y, p_atom->z, OC_CHARGE);
        }


	if(strcmp(p_O_elem, (p_atom->pot)) != 0 && strcmp(p_Al_elem, (p_atom->pot)) != 0 && strcmp(p_Si_elem, (p_atom->pot)) != 0 && strcmp(p_N_elem, (p_atom->pot)) != 0 && strcmp(p_N4_elem, (p_atom->pot)) != 0 && strcmp(p_H_elem, (p_atom->pot)) != 0 && strcmp(p_HC_elem, (p_atom->pot)) != 0 && strcmp(p_C1_elem, (p_atom->pot)) != 0 && strcmp(p_C2_elem, (p_atom->pot)) != 0 && strcmp(p_C3_elem, (p_atom->pot)) != 0  && strcmp(p_NZ_elem, (p_atom->pot)) != 0 && strcmp(p_F_elem, (p_atom->pot)) != 0 && strcmp(p_CL_elem, (p_atom->pot)) != 0 && strcmp(p_AR_elem, (p_atom->pot)) != 0 && strcmp(p_OO_elem, (p_atom->pot)) != 0 && strcmp(p_C_elem, (p_atom->pot)) != 0 && strcmp(p_OH_elem, (p_atom->pot)) != 0 && strcmp(p_HO_elem, (p_atom->pot)) != 0 && strcmp(p_CP_elem, (p_atom->pot)) != 0 && strcmp(p_SH_elem, (p_atom->pot)) != 0 && strcmp(p_HS_elem, (p_atom->pot)) != 0 && strcmp(p_OC_elem, (p_atom->pot)) !=0 && strcmp(p_NA_elem, (p_atom->pot)) !=0 )
        {
        printf("Error - a gulp file can not be created as there are incompatible atom/pots types!\n");
        fprintf(gulp_fp, "ERROR - a gulp file can not be created as there are incompatible atom/pots types!\n");
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
printf("in gulp routine- have produced input\n");


return;
}

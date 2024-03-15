/*********************************************************************/
/* MOPAC_minimize_molecule.c : minimizes a molecule with respect 	 */	
/*                       to itself using MOPAC                       */
/* started 19/6/ DWL                                                 */
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

#define MOPAC_SUCCESS 1
#define MOPAC_TOO_BIG -2

void delete_file(char *fileroot, char *extension);

int anal_mopac(char *root_file, atom *p_molecule, int num_atoms,
                                                    energy *p_energy);

int calculate_charge(atom *p_molecule, int num_atoms);

void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name);

int get_mopac_minimized(char *root_file,int num_skip_atoms,
                                                atom *p_molecule, int num_atoms);

void mopac_minimize_molecule(atom *p_molecule, int num_atoms, char *root_name,
                             int which_mol)
{

#include "header.h"
char fileroot[FILELEN_MAX];
char mopac_cmd[FILELEN_MAX];
char call_mopac[FILELEN_MAX];
FILE *mopac_fp;
char *fullstop;
int mopac_status;
atom *p_atom;
char dummy[80];

int i;
int charge;


/*************check we don't have too many atoms**************************/
if (num_atoms > max_mopac_atoms)
    {
    fprintf(output_fp, "ZEBEDDE WARNING: Too many atoms for MOPAC: %i sent %i is the max\n",
                num_atoms, max_mopac_atoms);
    fprintf(output_fp, "Abandoning Minimisation\n");
    fprintf(output_fp, "Trying minimization in default minimizer (%s)\n",
                                        MINIMIZER_DEFAULT);

	strcpy(dummy,minimizer_name);
	strcpy(minimizer_name,MINIMIZER_DEFAULT);
	minimize_molecule(p_molecule, num_atoms, root_name);
	strcpy(minimizer_name,dummy);
    return;
    }

/*************print out MOPAC input file***********************************/

strcpy(fileroot,root_name);
strcat(fileroot,".dat");
if (!(mopac_fp = fopen(fileroot, "w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening MOPAC input file %s in minimize_molecule\n", fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
	return;
    }

/***get the charge for the command line***/
charge = calculate_charge(p_molecule, num_atoms);

sprintf(mopac_cmd,"%s CHARGE=%i",mopac_cmdline_molecule, charge);

fprintf(mopac_fp,"%s\n\n",mopac_cmd);
fprintf(mopac_fp,"ZEBEDDE mopac_minimize_molecule\n");

for (i=0;i<=num_atoms;i++)
	{
	p_atom = p_molecule +i;
	fprintf(mopac_fp,"%s %12.8f 1 %12.8f 1 %12.8f 1\n", p_atom->elem, 
									p_atom->x, p_atom->y, p_atom->z);
	}
fclose(mopac_fp);

/******Run In MOPAC**************************/


sprintf(call_mopac,"%s %s",mopac_path,root_name);


system (call_mopac);	


/*******get results from mopac run *************************/

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

mopac_status = anal_mopac(fileroot,p_molecule,
						num_atoms, &interaction_energy[which_mol]);


if (mopac_status == MOPAC_SUCCESS)
	{
	get_mopac_minimized(fileroot, 0, p_molecule, num_atoms);
	}

else if (mopac_status == MOPAC_TOO_BIG) /***too many atoms***/
	{
        fprintf(output_fp,"\n\nZEBEDDE WARNING: Maximum number of atoms exceeded in MOPAC\n");
        fprintf(output_fp,"\nZEBEDDE WARNING: Switching to default minimiser (%s)\n",
                                                MINIMIZER_DEFAULT);
        fprintf(output_fp,"\nZEBEDDE WARNING: for remainder of job\n");
        fprintf(output_fp,"\nZEBEDDE WARNING: MOPAC should be recompiled to handle this\n");
        strcpy(minimizer_name,MINIMIZER_DEFAULT);


	}
else	
	{
        fprintf(output_fp,"Invalid Result from MOPAC. Ignoring run\n");
	}
/***tidy up a bit to minimize disk usage***/

delete_file(fileroot,".end");
delete_file(fileroot,".log");
delete_file(fileroot,".out");
delete_file(fileroot,".temp");
delete_file(fileroot,".syb");
delete_file(fileroot,".arc"); 
delete_file(fileroot,".end");
 
return;
}

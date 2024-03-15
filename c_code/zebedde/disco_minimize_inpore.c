/************************************************************************/
/* disco_minimize_inpore.c : minimizes a molecule with respect to a pore*/
/*                          currently writes out files for DISCOVER     */
/* started 2/5/95 DWL                                                   */
/* Set up to compile but not thought through how it should work for     */
/* multiple molecules. DJW March 2013.                                  */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_biosym_header(FILE *file_fp, int pbc_flag);

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void print_pbc_header(FILE *file_fp, double *p_abc);


void car_end(FILE *fp);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

int  get_disco_minimized(char *root_file,int mol_number,
                         atom *p_molecule, int num_atoms);

int anal_disco(char *root_file, atom *p_molecule, int num_atoms,
                                                    energy *p_energy);

void calculate_energy(atom *p_pore, int num_p_atoms,
                      atom *guest_ptrs[], int num_guests,
                      list_partition *p_guest_demarc,
                      double *p_kvecs, double *p_kvec2,
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad,
                      int have_comb_rules, int is_empty_pore);

void print_mdf(char *p_molname, atom *p_molecule, int num_atoms,
                                        int start, FILE *fp, int D_to_H);

void print_mdf_end(FILE *fp);

void print_mdf_pbc(FILE *fp);

void print_mdf_subset(char *p_molname, atom *p_molecule, int num_atoms,
                      int start, FILE *fp, int D_to_H, int have_subset);

void print_mdf_header(char *p_assembly, FILE *fp);

void delete_file(char *fileroot, char *extension);

void disco_minimize_inpore(atom *guest_ptrs[], int num_guests,
                           list_partition *p_guest_demarc,
                           atom *p_pore, char *root_name,
                           double *p_kvecs, double *p_kvec2, double *p_gvec2, 
                           int num_kvecs, double *p_cos_sum, double *p_sin_sum,
                           int *p_need_grad, double *p_grad, int have_comb_rules,
                           int is_empty_pore)
{
#include "header.h"
char fileroot[FILELEN_MAX];
char disc_shell[FILELEN_MAX];
FILE *disc_shell_fp;
char *fullstop;
char discover_cmdline[256];
int discover_status,get_status;
int iloop, D_to_H, imol, index, isymm;
char dummy_name[50];

list_partition *p_demarc;

strcpy(fileroot,root_name);
strcat(fileroot,".car");

if (!(discover_fp = fopen(fileroot, "w")))
	{
    fprintf(output_fp, "ZEBEDDE ERROR: Error opening Discover car file %s\n", fileroot);
    fprintf(output_fp, "Abandoning Minimisation\n");
    return;
    }

/******************************************************************/
/**** Write car file **********************************************/
/******************************************************************/
D_to_H=TRUE;

print_biosym_header(discover_fp, pbc);
print_frame_header("Zebedde Inpore minimisation", -1, discover_fp);
if (pbc) print_pbc_header(discover_fp, &abc[0]);
print_molecule(p_pore, num_pore_atoms, discover_fp, D_to_H);

fprintf(discover_fp,"end\n");
p_demarc= p_guest_demarc;
for (imol=0; imol < num_guests; imol++)
  {
    if (symm_set) index=imol*(num_symm_ops+2);
                                           else index=imol;

    for (isymm=0; isymm <= num_symm_ops; isymm++)
      {
        print_molecule(guest_ptrs[index], p_demarc->end, discover_fp, D_to_H);
        index++;
      }
    p_demarc++;
  }

car_end(discover_fp);

/******************************************************************/
/**** Write mdf file **********************************************/
/******************************************************************/

strcpy(fileroot,root_name);
strcat(fileroot,".mdf");
 
if (!(discover_fp = fopen(fileroot, "w")))
    {
    printf("ZEBEDDE ERROR: Error opening Discover mdf file %s\n", fileroot);
    printf("Abandoning Minimisation\n");
    return;
    }

print_mdf_header("pore_template",discover_fp);

print_mdf("pore", p_pore, num_pore_atoms, 0, discover_fp, D_to_H);

/*****write mdf for first template*****/
//print_mdf("template0", p_molecule, num_atoms, 0, discover_fp, D_to_H);

/****************************************************************************/
/***** Write the template images                                     ********/
/****************************************************************************/

//if (pbc)
//	{
//	for (iloop=0;iloop<=num_symm_ops;iloop++)
//		{
//		sprintf(dummy_name,"template%d",iloop);
//		/***just write out the same one cos they're all the same ***/
//print_mdf(dummy_name, p_molecule, num_atoms, 0, discover_fp, D_to_H);
//	}
// print_mdf_pbc(discover_fp);
//	}

//print_mdf_subset("pore", p_pore, num_pore_atoms, 0, discover_fp, D_to_H, TRUE);

/*****************************************************************************/
/***** NOW RUN IN DISCOVER AND RETRIEVE OUTPUT                           *****/
/***** call Discover script to runjob now at full on priority            *****/
/***** syntax: discover file_name discover_forcefield_name nice start_yn *****/
/***** syntax now:discover file_name AJWL 9/08 *******************************/
/*****************************************************************************/

strcpy(fileroot,root_name);

if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

sprintf(disc_shell,"./%s.csh",fileroot);


if (!(disc_shell_fp = fopen(disc_shell,"r")))
	{
	/* .disc_shell does not exist */
	/* so write it... */
	/* but delete any old stuff first to avoid prompting */
	delete_file(fileroot,".rst");
	delete_file(fileroot,".prm");
	delete_file(fileroot,".out"); 
	delete_file(fileroot,".cor");
	printf("Writing inpore.csh file for discover\n");
	sprintf(discover_cmdline,"%s %s\n",discover_path, fileroot);
	system (discover_cmdline);	
	}
/* .disc_shell exists so just run it */
//sprintf(discover_cmdline,"%s>/dev/null\n",disc_shell);
//fclose(disc_shell_fp);


//system (discover_cmdline);	


/*******get results from DISCOVER run *************************/

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

discover_status = anal_disco(fileroot, guest_ptrs[0],
                        p_guest_demarc->end, &interaction_energy[0]);

/*******DEBUG*******/
/* printf("Returned status %i  ",discover_status); */
/* if (discover_status == TRUE) */
    /* { */
    /* printf(" VALID Minimization\n"); */
    /* } */
/* else */
    /* { */
    /* printf(" INVALID Minimization\n"); */
    /* } */


/******* Was it a good move in DISCOVER ***********************/

/**** since we are minimising in the host, measure the non-bonding energy ***/
/**** accept move as long as energy decreased (even if not minimised***/
if (discover_status >= 0)
	{
	if (interaction_energy[0].minimizer_end_nonbond 
				< interaction_energy[0].minimizer_init_nonbond)
    	{

        /* update energy values */

//        calculate_energy(p_pore, num_pore_atoms, p_molecule,
//                         num_atoms, p_kvecs, p_kvec2, p_gvec2, num_kvecs, 
//                         p_cos_sum, p_sin_sum, p_need_grad, p_grad, have_comb_rules,
//                         is_empty_pore);


    	/* printf("INPORE Discover Minimisation did reduce the energy\n"); */
    	get_status = get_disco_minimized(fileroot,2, guest_ptrs[0], p_guest_demarc->end);

		if (get_status == -1 )
			{
			fprintf(output_fp, "ZEBEDDE ERROR: Error reading in minimized structure\n");
			fprintf(output_fp, "Abandoning Minimization\n");
			discover_status = FALSE;
			}
		}

	else
    	{
    	/* printf("INPORE Discover Minimisation did not reduce the energy\n");*/
    	discover_status = FALSE;
		}
    }

/***tidy up a bit to minimize disk usage***/

delete_file(fileroot,".rst");
delete_file(fileroot,".prm");
delete_file(fileroot,".out"); 
delete_file(fileroot,".cor");

return;
}

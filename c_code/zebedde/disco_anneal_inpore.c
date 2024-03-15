/***********************************************************************/
/****** disco_anneal_inpore.c : anneals molecule inpore         ********/
/******				after docking run               ********/
/****** Started AJWL 17/09/08					********/
/***********************************************************************/

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
                      atom *p_templ, int num_t_atoms,
                      double *p_kvecs, double *p_kvec2,
                      double *p_gvec2, int num_kvecs,
                      double *p_cos_sum, double *p_sin_sum,
                      int *p_need_grad, double *p_grad,
                      int have_comb_rules, int is_empty_pore);

void print_mdf(char *p_molname, atom *p_molecule, int num_atoms,
                                        int start, FILE *fp, int D_to_H);

void print_mdf_end(FILE *fp);
void print_mdf_pbc(FILE *fp);

void print_mdf_header(char *p_assembly, FILE *fp);

void delete_file(char *fileroot, char *extension);

void disco_anneal_inpore(atom *p_pore,
                           atom *p_molecule, int num_atoms, char *root_name,
                           double *p_kvecs, double *p_kvec2, double *p_gvec2,
                           int num_kvecs, double *p_cos_sum, double *p_sin_sum,
                           int *p_need_grad, double *p_grad, int have_comb_rules,
                           int is_empty_pore)

{
#include "header.h"
char fileroot[FILELEN_MAX];
char *fullstop;
char discover_cmdline[256];
int discover_status,get_status;
int iloop, D_to_H;
char dummy_name[40];
int iatom1,iatom2;
char this_group[10], that_group[10];
atom *p_atom1, *p_atom2;
atom *p_start_atom;
/*clean up a bit*/

strcpy(fileroot,root_name);

    delete_file(fileroot,".rst");
    delete_file(fileroot,".prm");
    delete_file(fileroot,".out");
    delete_file(fileroot,".cor");
	delete_file(fileroot,".mdf");

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
p_start_atom= p_molecule;
print_molecule(p_molecule, num_atoms, discover_fp, D_to_H);


if (pbc)
    {

/******************************************************************/
/********* write the template images ******************************/
/******************************************************************/

    for (iloop=1; iloop<=num_symm_ops; iloop++)
        {
        fprintf(discover_fp,"end\n");
                p_start_atom += (num_atoms+1);
        print_molecule(p_start_atom, num_atoms, discover_fp, D_to_H);
        }

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
print_mdf("template0", p_molecule, num_atoms, 0, discover_fp, D_to_H);

/****************************************************************************/
/***** Write the template images                                     ********/
/****************************************************************************/

if (pbc)
        {
        for (iloop=0;iloop<=num_symm_ops;iloop++)
                {
                sprintf(dummy_name,"template%d",iloop);
print_mdf(dummy_name, p_molecule, num_atoms, 0, discover_fp, D_to_H);
        }
 print_mdf_pbc(discover_fp);
        }

print_mdf_subset("pore", p_pore, num_pore_atoms, 0, discover_fp, D_to_H, TRUE);
return;
}

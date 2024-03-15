/*******************************************************************/
/* minimize_inpore.c :   Does a minimization on a molecule         */
/*                       in a fixed pore                           */
/*                       using a specifed minmizer                 */
/*                       Currently handles Discover and MOPAC      */
/* GULP minimiser added by Kim 2005                                */
/*  started Dewi 19/6                                              */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void car_molecule_output(atom *p_pore, atom *p_molecule,
                                         int num_atoms, char *root_name);

void disco_minimize_inpore(atom *p_pore,
                           atom *p_molecule, int num_atoms, char *root_name,
                           double *p_kvecs, double *p_kvec2,
                           double *p_gvec2, int num_kvecs,
                           double *p_cos_sum, double *p_sin_sum,
                           int *p_need_grad, double *p_grad,
                           int have_comb_rules, int is_empty_pore);
 
void disco_anneal_inpore(atom *p_pore,
                           atom *p_molecule, int num_atoms, char *root_name,
                           double *p_kvecs, double *p_kvec2,
                           double *p_gvec2, int num_kvecs,
                           double *p_cos_sum, double *p_sin_sum,
                           int *p_need_grad, double *p_grad,
                           int have_comb_rules, int is_empty_pore);


void mopac_minimize_inpore(atom *p_pore, atom *p_molecule, 
                                       int num_atoms, char *root_name);

void c2disco_minimize_inpore(atom *p_pore, atom *p_molecule, int num_atoms, 
                                       char *root_name, int which_mol);

void internal_minimize_inpore(atom *p_pore, int num_pore_atoms, atom *p_molecule,
                              int num_t_atoms, double *p_kvecs, double *p_kvec2,
                              double *p_gvec2, int num_kvecs,
                              double *p_cos_sum, double *p_sin_sum,
                              int *p_need_grad, double *p_grad,
                              int have_comb_rules, int is_empty_pore);

void gulp_inpore_justinput(atom *p_pore, atom *p_molecule, 
                                         int num_atoms, char *root_name);

void gulp_inpore_annealed(atom *p_pore, atom *p_molecule,
                                         int num_atoms, char *root_name);

void gulp_minimize_inpore(atom *p_pore, atom *p_molecule,
                                         int num_atoms, char *root_name);


void minimize_inpore(atom *p_pore, atom *p_molecule, 
                     int num_atoms, char *root_name,
                     double *p_kvecs, double *p_kvec2,
                     double *p_gvec2, int num_kvecs,
                     double *p_cos_sum, double *p_sin_sum,
                     int *p_need_grad, double *p_grad,
                     int have_comb_rules, int is_empty_pore,
                     int which_mol)
{
#include "header.h"
if(final_gulp==FALSE && final_annealed_gulp==FALSE && car_output == FALSE)
{
if (strcmp (minimizer_name, DISCOVER_MINIMIZER) ==0 && anneal_run == FALSE)
        {
		disco_minimize_inpore(p_pore, p_molecule, num_atoms, root_name, 
                              p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                              p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                              have_comb_rules, is_empty_pore);
        }

else if (strcmp (minimizer_name, DISCOVER_MINIMIZER) ==0 && anneal_run == TRUE)
        {
        disco_anneal_inpore(p_pore, p_molecule, num_atoms, root_name,
                              p_kvecs, p_kvec2, p_gvec2, num_kvecs,
                              p_cos_sum, p_sin_sum, p_need_grad, p_grad,
                              have_comb_rules, is_empty_pore);
        }


else if (strcmp (minimizer_name, C2DISCOVER_MINIMIZER) ==0)
        {
          printf("Calling c2disco_minimize_inpore\n");
          c2disco_minimize_inpore(p_pore, p_molecule, num_atoms, root_name, which_mol);
        }

else if (strcmp (minimizer_name, MOPAC_MINIMIZER) ==0)
        {
		if (pbc) 
                  {
                     fprintf(output_fp, "ZEBEDDE WARNING: No PBC facilities available for inpore minimisations (yet)\n");	
                  }
		else
                  {
                     mopac_minimize_inpore(p_pore, p_molecule, num_atoms, mopac_root);
                  }	
        }

else if (strcmp (minimizer_name, INTERNAL_MINIMIZER) ==0)
        {
          internal_minimize_inpore(p_pore, num_pore_atoms, p_molecule,
                                   num_atoms, p_kvecs, p_kvec2, p_gvec2, 
                                   num_kvecs, p_cos_sum, p_sin_sum,
                                   p_need_grad, p_grad, have_comb_rules, is_empty_pore); 
        }

else if (strcmp (minimizer_name, GULP_MINIMIZER) ==0)
        {
	  /*printf("Calling gulp_minimize_inpore\n");*/
	  gulp_minimize_inpore(p_pore, p_molecule, num_atoms, gulp_root);
        }


else
		{
		fprintf(output_fp, "SERIOUS ZEBEDDE ERROR: Minimizer Name not Recognized\n");
                fflush(output_fp);
		exit(-1);
		}
} /*if(!final_gulp)*/


if (final_gulp) 
	{
	  printf("in minimize_inpore, producing final input file\n");
	  gulp_inpore_justinput(p_pore, p_molecule, num_atoms, gulp_root);
	}
if(final_annealed_gulp)
	{
	  printf("producing final annealed gulp for template no: %d\n", no_annealed_temps);
	  gulp_inpore_annealed(p_pore, p_molecule, num_atoms, gulp_root);
	}	  
if (car_output)
	{
	car_molecule_output(p_pore, p_molecule, num_atoms, gulp_root);
	}
return;
}

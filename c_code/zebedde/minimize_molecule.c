/*******************************************************************/
/* minimize_molecule.c : Does a minimization on a free molecule    */
/*                       using a specifed minmizer                 */
/*                       Currently handles Discover and MOPAC      */
/*  started Dewi 19/6                                              */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void c2disco_minimize_molecule(atom *p_molecule, int num_atoms, 
                                         char *root_name, int which_mol);

void disco_minimize_molecule(atom *p_molecule, int num_atoms, 
                                         char *root_name, int which_mol);

void mopac_minimize_molecule(atom *p_molecule, int num_atoms, 
                                         char *root_name, int which_mol);

void gulp_minimize_molecule(atom *p_molecule, int num_atoms, 
                                         char *root_name, int which_mol); 

/*void internal_minimize_molecule(p_molecule, num_atoms, mopac_root,
                                  int *p_need_grad, double *p_grad)*/

void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name,
                       int *p_need_grad, double *p_grad, int which_mol)

{
#include "header.h"
printf("Inside minimize_molecule\n");
printf("Minimiser is %s\n", minimizer_name);
if (strcmp (minimizer_name, DISCOVER_MINIMIZER) ==0 && anneal_run ==FALSE)
        {
        disco_minimize_molecule(p_molecule, num_atoms, root_name, which_mol);
        }

else if (strcmp (minimizer_name, C2DISCOVER_MINIMIZER) ==0)
        {
        c2disco_minimize_molecule(p_molecule, num_atoms, root_name, which_mol);
        }

else if (strcmp (minimizer_name, MOPAC_MINIMIZER) ==0)
        {
        mopac_minimize_molecule(p_molecule, num_atoms, mopac_root, which_mol);
        }

else if (strcmp (minimizer_name, GULP_MINIMIZER) ==0) 
        {        
	printf("Going to use gulp minimiser\n");
	gulp_minimize_molecule(p_molecule, num_atoms, gulp_root, which_mol);
        }

else if (strcmp (minimizer_name, INTERNAL_MINIMIZER) ==0)
        {
        /*internal_minimize_molecule(p_molecule, num_atoms, mopac_root, p_need_grad, p_grad);*/
        }

else
		{
                fprintf(output_fp,"ZEBEDDE ERROR: Minimizer Name not Recognized\n");
                fflush(output_fp);
		exit(-1);
		}
printf("about to leave minimise molecule\n");
return;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "header.h"
#include "global_values.h"

/* protype list for this routine */

/*---------------------------------------------------------------------------*/


/* Work out template charge and make up pore charge to balance */

double neutralise_slab( atom *guest_ptrs[], int num_guests,
                        list_partition *p_guest_demarc, atom *pore)
{
  double tot_pore_charge;
  double tot_template_charge;

  atom *p_atom;

  int imol, iloop;
  int num_Si, num_O;
  int is_silaceous = FALSE;

/*-- Start off by getting total charges on template and pore -----------------*/

  tot_template_charge=0.0;
  for (imol=0; imol<num_guests; imol++)
    {
      for (iloop = 0; iloop < p_guest_demarc->num; iloop++ )
        {
          tot_template_charge += guest_ptrs[imol][iloop].part_chge; 
        }
    }

  num_Si=0;
  num_O=0;
  tot_pore_charge=0.0;
  for (iloop = 0; iloop <= num_pore_atoms; iloop++ )
    {
      tot_pore_charge += pore[iloop].part_chge; 

      if ( strcmp(pore[iloop].elem, "Si") == 0 ) num_Si++;
      if ( strcmp(pore[iloop].elem, "O") == 0 ) num_O++;
    }

  fprintf(output_fp,"Neutralising template charge:\n");
  fprintf(output_fp,"Current total charge on pore     : %10.6f\n", 
                                                     tot_pore_charge);
  fprintf(output_fp,"Current total charge on template : %10.6f\n\n", 
                                                 tot_template_charge);

  is_silaceous = num_Si + num_O == num_pore_atoms+1;
  if (is_silaceous)
    {
      fprintf(output_fp,"This pore is a simple SiO2 polymorph\n");
      fprintf(output_fp,"assigning formal charges.\n");

      tot_pore_charge=0.0;
      for (iloop = 0; iloop <= num_pore_atoms; iloop++ )
        {
          if ( strcmp(pore[iloop].elem, "Si") == 0 )
            {
/*****************************************/
/** If the template is neutral use    ****/
/** simple formal charges, if not     ****/
/** spread compensating charge evenly ****/
/** over silicons.                    ****/
/*****************************************/
               if (tot_template_charge >  0.001 ||
                   tot_template_charge < -0.001 )
                  {
                     pore[iloop].part_chge=4.0 - tot_template_charge/num_Si; 
                  }
               else
                  {
                     pore[iloop].part_chge=4.0; 
                  }
            }
          if ( strcmp(pore[iloop].elem, "O") == 0 ) pore[iloop].part_chge=-2.0;

          tot_pore_charge += pore[iloop].part_chge; 
        }
    }
  else
    {
       fprintf(output_fp,"ERROR: Cannot neutralise none standard pore structure\n");
       exit(0);
    }
  fprintf(output_fp, "Total charge on pore after balancing with");
  fprintf(output_fp," template  : %10.6f\n", tot_pore_charge);

  return tot_pore_charge+tot_template_charge;
}

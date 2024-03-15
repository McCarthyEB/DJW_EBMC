/******************************************************************************/
/**** print_seed_info.c : Prints out information regarding the seed being use */
/****                     at the start and during multi-seed runs             */
/**** Two routines:                                                           */
/****      print_seed_info     : global information                           */
/****      print_this_seed_info:  information  on this seed in multi-seed run */
/**** dwl 29/2                                                                */
/******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"


void print_dashes(int ndashes,FILE *fp);
void print_types( atom_number *p_types, int num_types,
                  int *p_hyd_list, int num_hyds, FILE *fp);


void print_seed_info(int seed_type, int num_of_seeds, int *p_seeds_to_use,
           int read_all_frames_as_seeds, int random_seed,
           int number_of_frames_in_arc,
		   int initial_hydrogen_weight, int increment_hydrogen_weight)
{
#include "header.h"

int i, have_ten;
int *p_seed;

print_dashes(80,output_fp);
printf("\nSeed Information\n");
fprintf(output_fp, "Seed Information\n");
print_dashes(30,output_fp);

fprintf(output_fp, "\nMaximum Number of Template to Grow : %d\n\n", max_templates);
fprintf(output_fp, "Template seed point will be located \n");

if (centralise_template)
    {
    fprintf(output_fp, "at the centre of mass of the pore\n\n");
    }
else
    {
    fprintf(output_fp, "at its given cartesian coordinates\n\n");
    }

if (seed_type == MOLE)
	{
	   fprintf(output_fp, "A molecule has been supplied as a seed:\n\n");
           fprintf(output_fp, "\tTemplate Seed File      : %s\n",seed_file);
           fprintf(output_fp, "\tTemplate Seed Title     : %s\n\n",seed_title);
        }
else if  (seed_type == FRAG)
    {
        fprintf(output_fp, "Template seed will be selected from the fragment library\n");
	if (random_seed == TRUE)
          {
            fprintf(output_fp, 
                    "The following list of seeds has been selected at random from the fragment library:\n");
          }
	else
          {
            fprintf(output_fp, "The following fragments will be used as seeds\n");
          }


        have_ten = 0;
        p_seed= p_seeds_to_use;

        for ( i=0; i < num_of_seeds; i++)
          {
             fprintf(output_fp, "%4d ",*p_seed );
             if (have_ten == 9)
               {
                 have_ten = 0;
                 fprintf(output_fp,"\n");
               }
             else
               {
                 have_ten++;
               }
           }

        fprintf(output_fp, "\n\n");
    }

else /** seed_type == ARCH  **/
    {
	printf("Seed type ARCH\n");
	fprintf(output_fp, "Template seed will be selected from the archive file %s \n", seed_file);
	fprintf(output_fp, "The following %d frames will be used as seeds\n", 
						num_of_seeds);
	for (i=0;i<num_of_seeds;i++)
           {
              fprintf(output_fp, "%d ",*p_seeds_to_use);
	      (*p_seeds_to_use)++;
            }
    }

fprintf(output_fp, "\n\n");

fprintf(output_fp, "During initialisation:\n");
fprintf(output_fp, "Template will be energy minimised      : ");

if (initial_minimize_template)
      {
      fprintf(output_fp, " ON\n");
      }
else
      {
      fprintf(output_fp, " OFF\n");
      }

fprintf(output_fp, "Template will be orientated in the Pore: ");
if (initial_minimize_inpore)
      {
      fprintf(output_fp, " ON\n");
      }
else
      {
      fprintf(output_fp, " OFF\n");
      }


print_dashes(40,output_fp);

fprintf(output_fp, "\nHydrogens will be weighted during growth:\n");
fprintf(output_fp, "\tInitial Hydrogen Weight = %d\n",initial_hydrogen_weight);
fprintf(output_fp, "\tIncrement on growth     = %d\n\n", increment_hydrogen_weight);
fprintf(output_fp,"\tNote: Negative numbers indicate that \"older\" hydrogens\n");
fprintf(output_fp,"\t      are *less* favoured\n");

print_dashes(80,output_fp);
fprintf(output_fp,"\n\n\n");

return;
}



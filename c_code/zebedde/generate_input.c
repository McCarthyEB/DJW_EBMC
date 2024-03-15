#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void generate_input(int number_of_actions, int *p_action_weights, 
                    int *p_frag_weights, int n_frag_weights) 
{
#include "header.h"
int i;
char filename[FILELEN_MAX];
char extraline[FILELEN_MAX];
char extrafile[FILELEN_MAX];
char oldline[FILELEN_MAX];

sprintf(filename,"new_input%i.dat",car_file_count);

if (!(docked_zeb_input_fp= fopen(filename,"w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening new ZEBEDDE input file : %s\n",
                      filename);
        return;
    }

/*print header stuff*/

fprintf(docked_zeb_input_fp,"# anything after a hash is now ignored\n# minimisers no longer work (interface to Accelrys uncertain!)\n#\n");

/*print the title*/
fprintf(docked_zeb_input_fp, "title\n");
fprintf(docked_zeb_input_fp, "%s\n",title);
fprintf(docked_zeb_input_fp,"initial_minimize template off\ninitial_minimize in_pore off\n");

fprintf(docked_zeb_input_fp,"#\n# No longer need biosym logfiles!\n#\n");
fprintf(docked_zeb_input_fp,"#box limits\n#-1.0 1.0 -2.0 2.0 -3.0 3.0\n#box fraction 0.75\n\n");

fprintf(docked_zeb_input_fp,"#to centre or not to centre initial seed or can use dock\n");

fprintf(docked_zeb_input_fp,"centre off\n\n");

fprintf(docked_zeb_input_fp,"#testing probability NOW A REAL PERCENTAGE\nprob_test %4.1f\n\n", prob_test*100);

fprintf(docked_zeb_input_fp,"#Rejection energy for docking run\ndokenergy %f\n\n", dock_energy);

if (exxon) fprintf(docked_zeb_input_fp,"exxon\n");

if (have_conc_limits == TRUE)
        {
	fprintf(docked_zeb_input_fp,"#concentration limits on each element\n# element (int) max conc in template\n");
	fprintf(docked_zeb_input_fp,"concentration %i\n", num_conc_limits);
	for (i=0; i< num_conc_limits;i++)
                {
                fprintf(docked_zeb_input_fp,"%s %i\n",
                    atom_limit[i].atom_type, atom_limit[i].num);
                }
        }


if (have_forbidden_bonds == TRUE)
	{
	fprintf(docked_zeb_input_fp,"# list of forbidden bonds\n");
	fprintf(docked_zeb_input_fp,"forbidden %i\n",num_forbidden_bonds);
		for (i=0; i< num_forbidden_bonds;i++)
      		{
        	fprintf(docked_zeb_input_fp,"%s %s\n",
        	forbidden_bond[i].atom1, forbidden_bond[i].atom2);
        	}
	}


fprintf(docked_zeb_input_fp,"#initial increment\nweight hydrogen\n1  50\n\n");

fprintf(docked_zeb_input_fp,"#set the minimizer\n");
fprintf(docked_zeb_input_fp,"minimizer %s\n",minimizer_name);

/*Should do final gulp?*/

fprintf(docked_zeb_input_fp,"#if you also want a final gulp file for each template produced, include the word fgul\n");
fprintf(docked_zeb_input_fp,"fcar\n");
fprintf(docked_zeb_input_fp,"anneal\n");

fprintf(docked_zeb_input_fp,"minimizer discover\n");
fprintf(docked_zeb_input_fp,"strategy template %s\n",template_strategy_file);
fprintf(docked_zeb_input_fp,"strategy in_pore %s\n",inpore_strategy_file);
fprintf(docked_zeb_input_fp,"discover path %s\n",discover_path);
fprintf(docked_zeb_input_fp,"discover forcefield %s\n",discover_forcefield_name);

fprintf(docked_zeb_input_fp,"#fileroot for gulp files\n");
fprintf(docked_zeb_input_fp,"gulp name final_docked\n");

fprintf(docked_zeb_input_fp,"#first commandline is for template, second for inpore minimization\n");
fprintf(docked_zeb_input_fp,"gcommandline\n");
fprintf(docked_zeb_input_fp,"%s",gulp_cmdline_molecule);
fprintf(docked_zeb_input_fp,"%s",gulp_cmdline_inpore);

/*GULP extra line stuff*/

fprintf(docked_zeb_input_fp,"#extra commandlines for gulp - it reads lines until \"end extralines for gulp\"\n");
fprintf(docked_zeb_input_fp,"gextra lines for gulp\n");
strcpy(extrafile, "extra_gulp.lines");
extra_gulp_fp = fopen(extrafile, "r");

while (!feof(extra_gulp_fp))
   {
    fgets(extraline,FILELEN_MAX,extra_gulp_fp);
    if( strcmp(extraline, oldline) != 0)
        {
        fprintf(docked_zeb_input_fp,"%s", extraline);
        }
        strcpy(oldline, extraline);
   }

fclose(extra_gulp_fp);

fprintf(docked_zeb_input_fp,"#end extralines for gulp\n");


fprintf(docked_zeb_input_fp,"#potentials file for gulp\n");
fprintf(docked_zeb_input_fp,"g_pots %s\n", gulp_pots_file);
fprintf(docked_zeb_input_fp,"gulp path %s\n",gulp_path);

fprintf(docked_zeb_input_fp,"#growth criteria\n");
fprintf(docked_zeb_input_fp,"vdw_scale %f\n",vdw_scale);
fprintf(docked_zeb_input_fp,"stop %f\n",stop_ctf);

fprintf(docked_zeb_input_fp,"#files\n");
fprintf(docked_zeb_input_fp,"pore %s\n",pore_file);

fprintf(docked_zeb_input_fp,"\nconvention XYZ\n");


fprintf(docked_zeb_input_fp,"#\n# verbose [on/off/debug] debug monitors which actions are being carried out and is extremely detailed\n#");
fprintf(docked_zeb_input_fp,"verbose\n");


fprintf(docked_zeb_input_fp,"seed molecule %s\n", temp_car_output);

fprintf(docked_zeb_input_fp,"\n#modification parameters\n");

fprintf(docked_zeb_input_fp,"modify attempts %d\n", num_modify_attempts);
fprintf(docked_zeb_input_fp,"maxtemplates %d\n\n",max_templates);
/*
fprintf(docked_zeb_input_fp,"modify attempts %i\n", num_modify_attempts);
fprintf(docked_zeb_input_fp,"maxtemplates %i\n\n",max_templates);
*/
fprintf(docked_zeb_input_fp,"# Cost function directive\n# cost_function [steric/non_bond/nonbond/energy]\n#\n\
# ONLY use non_bond OR steric at the moment\n#\n\n");

if (steric)
	{
	fprintf(docked_zeb_input_fp,"cost_function steric\n");
	}

if (non_bonded)
	{
	fprintf(docked_zeb_input_fp,"cost_function non_bond\n");
	}


fprintf(docked_zeb_input_fp,"nb_cutoff %f\n",nb_ctf);


fprintf(docked_zeb_input_fp,"#\n#General weighting of items\n");
fprintf(docked_zeb_input_fp,"#Note need for end statement to delimit the weights lines.\n#\n\
#-------------------------------------------------------------------------\n\
# Build | Rotate_last_bond | Shake | Rock | Twist | Ring_formation | Minimise_gas_phase | Minimise_in_host\n");
fprintf(docked_zeb_input_fp,"weights action\n");
/*fprintf(docked_zeb_input_fp,"0\t0\t10\t50\t50\t0\t0\t0\n");*/
for (i=0;i<number_of_actions; i++)
fprintf(docked_zeb_input_fp,"%i             ",*(p_action_weights+i));
fprintf(docked_zeb_input_fp,"\nend\n");


fprintf(docked_zeb_input_fp,"#\n# Temperature used in Monte Carlo testing\n#\n");
fprintf(docked_zeb_input_fp,"temperature %f\n",temperature);


fprintf(docked_zeb_input_fp,"#if we are carrying out an optimisation after each mc run, need keyword mcop\n\
#mcop\n\
#will then increment temperature, by tstep (a double)\n\
#tstep 50.0\n\
#carrying out amodi actions at each temperature step\n\
#amodi  6000\n\n\
#Interval at which to produce a current template coordinate file (.pcar) file\n\
# peek now also controls the trajectory writing:\n\
# when the number of attempt cycles hits the value set the next\n\
# accepted structure is written. In the MC docking mode lots of\n\
# structures can appear so this is to limit the size of the trajectory\n\
# file. The energy is still output to the out file for every success.\n\
#\n");
fprintf(docked_zeb_input_fp,"peek 100\n");

fprintf(docked_zeb_input_fp,"#fragment library control\n");
fprintf(docked_zeb_input_fp,"library fragment %s\n",fragment_file);

fprintf(docked_zeb_input_fp,"#library for internal forcefield\n");

fprintf(docked_zeb_input_fp,"library forcefield %s\n\n",forcefield_library);

fprintf(docked_zeb_input_fp,"#methane  ethane  benzene  ammonia  ammonium propane pyrrole adamantane cyclohexane\n");
fprintf(docked_zeb_input_fp,"weights fragments\n");

for (i=0; i < n_frag_weights; i++)
   {
     sum_frag_weights += *p_frag_weights;
     fprintf(docked_zeb_input_fp,"%d         ",*p_frag_weights);
     p_frag_weights++;
   }
fprintf(docked_zeb_input_fp,"\nend\n");



if (num_allowed_torsions >=0)
        {
        fprintf(docked_zeb_input_fp,"\n#\n# Allowed flexible torsions\n#\n");
        fprintf(docked_zeb_input_fp,"torsions %i\n", num_allowed_torsions+1);
      for (i=0; i<=num_allowed_torsions; i++)
        {
          fprintf(docked_zeb_input_fp,"%s  %s  %s  %s\n",allowed_torsions[i].A,
                                                  allowed_torsions[i].B,
                                                  allowed_torsions[i].C,
                                                  allowed_torsions[i].D);
        }
        }

sprintf(zebedde_input_filename, filename);
fclose(docked_zeb_input_fp);
return;
}

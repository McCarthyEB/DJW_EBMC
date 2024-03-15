/***********************************************************************/
/**output_routines.c ; general output routines                         */
/*                                                                     */
/* started 6.4.95 dewi                                                 */
/***********************************************************************/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"
#include "header.h"

/***********************************************************************/
/*********writebackinput ; write out user input and start parameters   */
/***********************************************************************/

void print_dashes(int ndashes,FILE *fp);
void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);
void banner(FILE *output_fp);
void print_types( atom_number *p_types, int num_types,
                  int *p_hyd_list, int num_hyds, FILE *fp);

void write_back_input(atom_number *p_frag_types, 
                      atom *pore, atom *frag_lib,
                      int *frag_weights,
                      list_partition *p_frag_types_list,
                      list_partition *p_frag_hyd_partition,
                      int *p_frag_hyd_list, char *p_animation_file,
					  int test_symm,
                      char *p_analyse_name,
                      int randomisation_method, int is_empty_pore)

{

time_t idate;
int start,ifrag;
int i,j;

fprintf(output_fp, "\n\n\n");
fprintf(output_fp, "Job running from : %s\n\n", inputfile);
fprintf(output_fp, "Job Title        : \n");
fprintf(output_fp, "                 : %s\n\n", title);
idate = time(NULL);
//fprintf(output_fp, "Job Started      : %24.24s\n\n", ctime(&idate));

if (verbose) 
	{
	  fprintf(output_fp, "Verbose Output Requested\n\n");
	}
else
	{
	  fprintf(output_fp, "Terse Output Requested\n\n");
	}
print_dashes(80,output_fp);

if (test_symm)
	{
	print_dashes(47,output_fp);
	fprintf(output_fp, "|Input Check Run - Will terminate before Growth|\n");
	print_dashes(47,output_fp);
	}

/************************Geometry input******************************/

fprintf(output_fp,"\n\nPore file     : %s\n",  pore_file);
fprintf(output_fp,"Pore title    : %s\n",  pore_title);

if (verbose) 
/* print out the pore coordinates */
     {
       if (!is_empty_pore)
         {
           fprintf(output_fp,"Pore has %i atoms\n\n",num_pore_atoms+1);
           fprintf(output_fp, "Pore Coordinates :\n");
           print_dashes(80,output_fp);
           print_molecule(&pore[0],num_pore_atoms,output_fp, FALSE);
	}
      else
        {
           fprintf(output_fp, "There are no pore coordinates," );
           fprintf(output_fp, "the pore file simply defines periodic boundary conditions.\n");
        }
      fprintf(output_fp, "\n\n");
     }

/*********seed stuff commented out : new routine! *********/
/*

fprintf(output_fp, "Maximum Number of Template to Grow : %i\n\n",max_templates);

fprintf(output_fp, "Template seed point will be located: \n");
if (centralise_template)
	{
	fprintf(output_fp, "\tat the centre of mass of the pore\n\n");
	}
	else
	{
	fprintf(output_fp, "\tat its given cartesian coordinates\n\n");
	}

if (seed_fp != NULL) 
	{
	fprintf(output_fp, "Template Seed File      : %s\n",seed_file);
	fprintf(output_fp, "Template Seed Title     : %s\n\n",seed_title);
	if (verbose) 
		{
		fprintf(output_fp, "Seed Coordinates :\n");
		print_dashes(80,output_fp);
		print_molecule(p_template, num_template_atoms,output_fp, FALSE);
		print_types(p_template_types, num_template_types,
                    p_template_hyd_list, num_template_hyds, output_fp);
		}
	}
else 
	{
	fprintf(output_fp, "No Template Seed supplied : \n");
	fprintf(output_fp, "Template seed will be selected from the fragment library\n");
	}

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

*/
/*********seed stuff commented out : new routine! *********/

print_dashes(80,output_fp);

fprintf(output_fp,"\nFragment Library for build : %s\n",fragment_file);
fprintf(output_fp,"Fragment Library contains %i entries\n\n",
					number_of_fragments);
fprintf(output_fp,"Fragment Library contains %i atoms\n\n",
					total_frag_atoms);


for (ifrag=1;ifrag<=number_of_fragments;ifrag++)
	{
	fprintf(output_fp, "\nFragment %i has %i atoms.  Fragment Weight = %i\n",
						ifrag,number_of_members[ifrag],frag_weights[ifrag-1]);

	if (verbose)
		/* print fragment library coordinates */
		{
		fprintf(output_fp,"\nCoordinates\n");
		print_dashes(20,output_fp);
		start =  member_start[ifrag];
		print_molecule(&frag_lib[start], number_of_members[ifrag],output_fp, FALSE);

        print_types(p_frag_types+((p_frag_types_list+ifrag)->start),
                    (p_frag_types_list+ifrag)->num,
                    p_frag_hyd_list+((p_frag_hyd_partition+ifrag)->start),
                   (p_frag_hyd_partition+ifrag)->num,
                    output_fp);

		print_dashes(80,output_fp);
		}

	}


/*******************Other supplied files ****************************/

fprintf(output_fp,"\nMinimisation Control\n");
print_dashes(30,output_fp);

/* fprintf(output_fp,"\nMinimizer to be used : %s\n\n", minimizer_name); */


/*******************Minimizer files**********************************/

if (strcmp(minimizer_name,DISCOVER_MINIMIZER) == 0)
	{
	/*******************DISCOVER files***********************************/
	fprintf(output_fp,"\nDISCOVER being used for Minmizations\n");
	fprintf(output_fp,"\nDiscover executable is %s\n", discover_path);
	fprintf(output_fp,"\nForcefield library        : %s\n\n",discover_forcefield_name);
	fprintf(output_fp,"\n\tTemplate Minimizer Strategy file: %s\n",template_strategy_file);
	fprintf(output_fp,"\n\tIn Pore Minimizer Strategy file : %s\n\n",inpore_strategy_file);
	fprintf(output_fp,"Intermediate files to be used by DISCOVER\n\n");
	fprintf(output_fp,"Minimisation of Template : \n\t\t\t%s and \n\t\t\t%s\n\n",
						template_min_car, template_min_mdf);
	fprintf(output_fp,"In Pore Minimisation     : \n\t\t\t%s and \n\t\t\t%s\n\n" ,inpore_min_car, inpore_min_mdf);

	}

else if (strcmp(minimizer_name,MOPAC_MINIMIZER) == 0)
	{
	fprintf(output_fp,"\nMOPAC being used for Minmizations\n");
	fprintf(output_fp,"\nMOPAC executable is %s\n", mopac_path);
	fprintf(output_fp,"\n\tCommand line for Molecule Minimization:\n");
	fprintf(output_fp,"\t\t%s\n",mopac_cmdline_molecule);
	fprintf(output_fp,"\n\tCommand line for In_Pore Minimization:\n");
	fprintf(output_fp,"\t\t%s\n",mopac_cmdline_inpore);
	fprintf(output_fp,"MOPAC runs to be called: %s\n",mopac_root);
	}

else if (strcmp(minimizer_name,GULP_MINIMIZER) == 0)
        {
        fprintf(output_fp,"\nGULP being used for Minimizations\n");
        fprintf(output_fp,"\nGULP executable is %s\n", gulp_path);
        fprintf(output_fp,"\n\tCommand line for Molecule Minimization:\n");
        fprintf(output_fp,"\t\t%s\n",gulp_cmdline_molecule);
        fprintf(output_fp,"\n\tCommand line for In_Pore Minimization:\n");
        fprintf(output_fp,"\t\t%s\n",gulp_cmdline_inpore);
        fprintf(output_fp,"GULP runs to be called: %s\n",gulp_root);
        }

else if (strcmp(minimizer_name,INTERNAL_MINIMIZER) == 0)
	{
	fprintf(output_fp,"\nInternal Minimizer Will be Used\n\n");
	}

/*******************Output Options **********************************/
print_dashes(80,output_fp);
fprintf(output_fp,"\n***********Output Control***********\n");
fprintf(output_fp,"\n\nTemplates that are accepted will be written to: %s.arc\n\n",
			gooduns);

if (animate_flag != 0)
	{
	fprintf(output_fp,"\nAnimations to be Generated\n");
	print_dashes(40,output_fp);
	fprintf(output_fp,"\nAnimation file            : %s\n",p_animation_file);
	fprintf(output_fp,"\nAnimation format          :");

	/****** multiple animation formats 6.6.96 DWL *********/
	if (animate_flag == BIOSYM_ANIMATION) 
			{
			fprintf(output_fp," Biosym archive\n");
			}
	else if (animate_flag == XMOL_ANIMATION_WITHPORE) 
			{
			fprintf(output_fp," Xmol XYZ animation with Pore\n");
			}
	else if (animate_flag == XMOL_ANIMATION_NOPORE) 
			{
			fprintf(output_fp," Xmol XYZ animation with no Pore\n");
			}
	else
			{
			fprintf(output_fp,"ZEBEDDE ERROR: Unknown Animation type %d\n", animate_flag);
			fprintf(output_fp,"               Exiting\n");
			exit(EXIT_FAILURE);
			}
	}
fprintf(output_fp,"\n***********Output Control***********\n");
print_dashes(80,output_fp);

/*******************Random Number   control**************************/
fprintf(output_fp,"\nRandomisation of Number Generator\n");
print_dashes(40,output_fp);
fprintf(output_fp,"\nGenerator will ");
if (randomisation_method == 0)
	{
	fprintf(output_fp,"be seeded from clock (i.e. \"really\" random)\n");
	}
else
	{
	fprintf(output_fp,"not be seeded (i.e. same random sequence will repeat)\n");
	}
/*******************Energy function control**************************/

fprintf(output_fp,"\nGrowth Acceptance Evaluation Criteria\n");
print_dashes(40,output_fp);

fprintf(output_fp,"\nSTERIC");
if (steric) 
	{
	fprintf(output_fp,"\nVan der Waals scaling factor = %5.4f\n\n",vdw_scale);
	}
else
	{
	fprintf(output_fp," Not Included\n\n");
	}
fprintf(output_fp,"\nNON_BONDED");
if (non_bonded) 
	{
	fprintf(output_fp," with a cutoff of %5.4f Angstrom\n",nb_ctf);
	}
else
	{
	fprintf(output_fp,"  Not Included\n\n");
	}
fprintf(output_fp,"\nCHARGE INTERACTION");
if (charges) 
	{
	fprintf(output_fp," Using Ewald Sum\n");
	}
else
        {
	fprintf(output_fp,"  Not Included\n\n");
        }
if (have_tethers)
        {
          fprintf(output_fp,"\nAtoms labeled %s in pore and %s in seed will be restrained by a \n",
                  tether.A, tether.B);
          fprintf(output_fp,"harmonic tethering potential with\nspring constant %10.6f ", tether.k);
          fprintf(output_fp,
   "(kcal/mol per sq Angstrom)\nminimum energy distance %10.6f (Angstroms)\n", 
                       tether.r0);
        }



fprintf(output_fp,"\n");
print_dashes(80,output_fp);

fprintf(output_fp,"Modification Control Parameters\n");
print_dashes(30,output_fp);

fprintf(output_fp,"\nClose Contact Sum Cutoff Parameter = %6.4f Angstrom\n", stop_ctf);

fprintf(output_fp,"\nMaximum number of modification attempts per template = %i\n", num_modify_attempts);

fprintf(output_fp,"\nTest Probability  = %4.1f %%\n", 100*prob_test);

fprintf(output_fp,"\nSHAKE move:  max attempts = %4i step = %6.4f Angstrom\n", num_shake_attempts, max_shake_step);

fprintf(output_fp,"\nROCK move :  max attempts = %4i step = %6.4f degree\n",
			 num_rock_attempts, RAD_TO_DEG*max_rock_step);

fprintf(output_fp,"\nRING Formation Control:\n");
fprintf(output_fp,"Will join 5th order neighbours within a cutoff of %5.4f Angstrom, (%5.4f sqrd)\n",ring_ctf,ring_ctf_2);

if (have_conc_limits == TRUE)
	{
	fprintf(output_fp,"\n\nElemental Concentration will be Limited as follows\n");
	print_dashes(52,output_fp);
	fprintf(output_fp,"\nElement  Limit\n");
	fprintf(output_fp,"--------------\n");
	for (i=0; i< num_conc_limits;i++)
		{
		fprintf(output_fp,"%-4s  %4i\n", 
                    atom_limit[i].atom_type, atom_limit[i].num);
		}
	fprintf(output_fp,"--------------\n");
	fprintf(output_fp,"\n");
	}

if (have_forbidden_bonds == TRUE)
   {
      fprintf(output_fp,"\n\nBonds are Forbidden Between the Following Atoms\n");
      fprintf(output_fp,"\n--------------\n");

      for (i=0; i< num_forbidden_bonds;i++)
        {
           fprintf(output_fp,"  %4s and %-4s\n", 
                   forbidden_bond[i].atom1, forbidden_bond[i].atom2);
        }

      fprintf(output_fp,"--------------\n");
      fprintf(output_fp,"\n");
   }

if (force_dihedrals == TRUE)
   {
      fprintf(output_fp,"\n\nIf the following dihedrals are formed will force the dihedral angles given\n");
      fprintf(output_fp,"\n--------------\n");

      for (i=0; i< num_for_dihedrals;i++)
        {
           fprintf(output_fp,"  %4s %4s %4s %4s %10.6f degrees\n",
                   set_diherals[i].A, set_diherals[i].B, set_diherals[i].C, set_diherals[i].D,
                   set_diherals[i].phi); 
        }

      fprintf(output_fp,"--------------\n");
      fprintf(output_fp,"\n");
   }


fprintf(output_fp,"\n");
print_dashes(80,output_fp);

fprintf(output_fp,"Template Symmetry Parameters\n");
print_dashes(40,output_fp);

if (!symm_set)
	{
	fprintf(output_fp,"\nNo Symmetry constraints for the Template\n");
	}
else
	{
	fprintf(output_fp,"\n%i Symmetry Equivalent Templates", num_symm_ops+1);
	fprintf(output_fp," will be generated\n\n");
	fprintf(output_fp,"Symmetry Operators Are:\n\n");
	fprintf(output_fp,"          Rotation        Translation\n");
	print_dashes(80,output_fp);
	for(i=0;i<=num_symm_ops;i++)
		{
		fprintf(output_fp,"Operator %i\n",i+1);
		print_dashes(12,output_fp);
		for (j=0;j<=2;j++)
			{
			fprintf(output_fp,"        | %6.4f %6.4f %6.4f |     | %6.4f |\n",
				symm[i].matrix[(j*3)],
                symm[i].matrix[(j*3)+1],
                symm[i].matrix[(j*3)+2],
				symm[i].translation[j]);
			}
		fprintf(output_fp,"\n");
		}
	}
fprintf(output_fp,"\n");
print_dashes(80,output_fp);

fprintf(output_fp,"Growth Box Parameters\n");
print_dashes(30,output_fp);
if (pbc)
	{
	fprintf(output_fp,"\nPeriodic Pore Supplied");
	fprintf(output_fp,"\nCell Parameters are:\n");
	fprintf(output_fp,"a     = %10.6f b    = %10.6f c     = %10.6f\n",
			abc[0],abc[1],abc[2]);
	fprintf(output_fp,"alpha = %10.6f beta = %10.6f gamma = %10.6f\n",
			abc[3],abc[4],abc[5]);
	}
	else
	{
	if (user_box == TRUE) 
		{
		fprintf(output_fp,"\nUser Supplied box limits :\n");
		}
	else
		{
		fprintf(output_fp,"\nBox calculated from Pore Extents\n");
		fprintf(output_fp,"Using a box of %4.3f of the pore extents\n\n",
								box_fraction);
		
		}
	fprintf(output_fp,"      -x : %8.4f  +x : %8.4f\n",
								box_limits[0],box_limits[1]);
	fprintf(output_fp,"      -y : %8.4f  +y : %8.4f\n",
								box_limits[2],box_limits[3]);
	fprintf(output_fp,"      -z : %8.4f  +z : %8.4f\n",
								box_limits[4],box_limits[5]);
	
	fprintf(output_fp,"\nNOTE: If a user box and a box fraction were supplied, the user box is used.\n\n");
	}
print_dashes(30,output_fp);

/*******DWL 3/7/96 Added output for analyse run ********/
if (strcmp(p_analyse_name, NO_ANALYSE) != 0)
	{
	fprintf(output_fp,"\n\n    This is an Analysis Run\n");
	fprintf(output_fp,"    =======================\n");
        fprintf(output_fp,"\nArchive to be analysed: %s\n\n\n", p_analyse_name);
	}
print_dashes(80,output_fp);

return;
}

/**************************************************************************/
/* output non-bonding parameters in a structure                           */
/* 10.4.95 DWL                                                            */
/**************************************************************************/

void nb_print(atom *p_molecule, int num_atoms)

	{
#include "header.h"
	int i;
	atom *p_atom;

	fprintf(output_fp,"\nLabel  Elem  Read_pot Use_pot    VdW rad");
        if (non_bonded == TRUE)
             {
             fprintf(output_fp,"             A             B\n");
             }
        else
             {
             fprintf(output_fp,"\n");
             }
	print_dashes(80,output_fp);

	for (i=0;i< num_atoms;i++)
		{
		p_atom = p_molecule+i;
		fprintf(output_fp,"%-5s   %-4s   %-4s   %-4s       %6.4f",
						p_atom->label, p_atom->elem, p_atom->pot,
						potent[p_atom->nb_list].pot,
						p_atom->vdw);
		if (non_bonded == TRUE)
			{
			fprintf(output_fp,"     %12.3f %12.3f\n",
				potent[p_atom->nb_list].a,potent[p_atom->nb_list].b);
			}
		else 
			{
			fprintf(output_fp,"\n");
			}
						
		}
	print_dashes(80,output_fp);
	return;
	}

/**************************************************************************/
/* print_dashes                                                           */
/*  prints n_dashes dashes followed by a new line                         */
/**************************************************************************/

void print_dashes(int ndashes,FILE *fp)
	{
#include "header.h"
	int i;
	for (i=0;i<ndashes;i++)
		 fprintf(fp,"-");
    fprintf(fp,"\n");
	return;
}

/**************************************************************************/
/* epilogue                                                             */
/*   closes up the files and stuff                                        */
/**************************************************************************/
void print_dashes(int ndashes,FILE *fp);
void timer(int stop_start);

void epilogue(char *p_this_animation_file)
{

#include "header.h"
time_t idate;

/**************** close up the necesary files ***************/

/* if (anim_file_fp) fclose(anim_file_fp); */

if (animate_flag== BIOSYM_ANIMATION)
     {
        if (logfile_needed)
          {
            fprintf(anim_read_fp,"# End of read\n");
            fclose(anim_read_fp);
            fprintf(anim_show_fp,"display on anim*\n");
            fprintf(anim_show_fp,"# End of show\n");
            fclose(anim_show_fp);
          }

	print_dashes(80,output_fp);
	fprintf(output_fp,"\nAnimation file %s contains %i frames\n\n",
			p_this_animation_file,num_anime_frames);
        if (logfile_needed)
          {
	    fprintf(output_fp,"Animations can be read into Insight by sourcing files:\n");
       	    fprintf(output_fp,"%s_read.log and %s_show.log\n", 
				p_this_animation_file,p_this_animation_file);
          }
	}
else if (animate_flag== XMOL_ANIMATION_WITHPORE || 
		 animate_flag== XMOL_ANIMATION_NOPORE)
	{
	fprintf(output_fp,"\nAnimation file %s contains %i frames\n\n",
            p_this_animation_file,num_anime_frames);
	fprintf(output_fp,"Animations can be read into XMOL as XYZ files\n");
	}

/******lob out the end time*****/
idate = time(NULL);
fprintf(output_fp, "\n\n");
print_dashes(80,output_fp);
printf( "               Cycle Finished     : %24.24s\n", ctime(&idate));
fprintf(output_fp, "               Cycle Finished     : %24.24s\n", ctime(&idate));
print_dashes(80,output_fp);
}


void print_energy(FILE *fp, energy *p_energy, internal_energy *p_intra, atom *p_molecule, 
					int num_atoms, int madalung, int have_built,
                                        int need_mini_ene, int guest_no, int build_no)
{

#include "header.h"
int iatom;
atom *p_atom;

print_dashes(60,fp);
fprintf(fp, "Internal Energies for Template (per molecule)\n");
fprintf(fp, "Bond Stretch : %10.6f\n", p_intra->stretch);
fprintf(fp, "Angle bend   : %10.6f\n", p_intra->angle  );
fprintf(fp, "Torsion      : %10.6f\n", p_intra->torsion);
fprintf(fp, "Van der Waals: %10.6f\n", p_intra->vdw);
fprintf(fp, "===================================\n");
fprintf(fp, "Total intra  : %10.6f\n", p_intra->total);


fprintf(fp, "Template host interaction energy\n");
fprintf(fp,"Steric Energy Flag : %i\n",p_energy->steric);
fprintf(fp,"Non-bonded Energy  : %lf\n",p_energy->non_bonded);
fprintf(fp,"Coulombic Energy   : %lf\n",p_energy->charges);
fprintf(fp,"Hbond Energy       : %lf\n",p_energy->hbond);
fprintf(fp,"Restraint Energy   : %lf\n",p_energy->restraint);
fprintf(fp,"Total interaction  : %lf for guest molecule %d for build step %d\n", 
                             p_energy->non_bonded+p_energy->charges+
                             p_energy->restraint+p_energy->hbond,
                             guest_no, build_no);
fprintf(fp, "===================================\n");

if (have_built)
 {
   fprintf(fp,"Total intra+inter  : %lf after successful build. for guest molecule %d for build step %d\n", 
                             p_energy->non_bonded+p_energy->charges+
                             p_energy->restraint+p_energy->hbond+
                             p_intra->total, 
                             guest_no, build_no);
 }
else
 {
   fprintf(fp,"Total intra+inter  : %lf\n for guest molecule %d for build step %d\n", 
                             p_energy->non_bonded+p_energy->charges+
                             p_energy->restraint+p_energy->hbond+
                             p_intra->total, 
                             guest_no, build_no);
 }

/* Should flag if not needed */
if (need_mini_ene)
  { 
    fprintf(fp,"\nEnergies from Minimizer (last run)\n");
    fprintf(fp,"\t\t\t\t\t\t\tinitial\t\tend\n");
    fprintf(fp,"Total Energy      : \t%10.5f\t%10.5f\n",
					p_energy->minimizer_init_total,
					p_energy->minimizer_end_total);
    fprintf(fp,"Non-bonded Energy : \t%10.5f\t%10.5f\n",
					p_energy->minimizer_init_nonbond,
					p_energy->minimizer_end_nonbond);
    print_dashes(60,fp);
  }

if (madalung)
  {

      fprintf(fp,"Atomic contributions to non-bond potential\n\n");

      fprintf(fp,"Atom No  type  vdw energy\n\n");

      p_atom= p_molecule-1;
      for (iatom=0; iatom <= num_atoms; iatom++)
         {
           p_atom++;
           fprintf(fp,"  %-4d    %-4s    %10.6f\n",
                                 iatom+1, p_atom->label, p_atom->vdw_energy); 
         }

      print_dashes(60,fp);
   }

return;
}

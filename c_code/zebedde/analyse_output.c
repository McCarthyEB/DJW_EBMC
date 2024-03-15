/*********************************************************************/
/* analyse_output.c : Takes an archive file, minimises the molecules */
/*                    in the host and in the gas phase, calculates   */
/*                    binding energy, sorts and writes new archive   */
/*                    currently uses DISCOVER                        */
/*                    but should use whatever MINIMISER is set       */
/* started 3/7/96 DWL                                                */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "maxima.h"
#include "structures.h"
#include "data.h"
#include "global_values.h"

void do_symmetry( atom *p_molecule, int num_atoms, int num_sym_ops,
                  symm_ops *p_symm, int *p_num_with_symm);

void index_sort(int n, float *p_arrin  ,int *p_indx);

void car_end(FILE *fp);

void generate_neighbours( atom *p_molecule, int num_atoms,
                          atom_number *p_types, int *p_num_types,
                          int use_pbc);

void strip_spaces(char *string);

void print_biosym_header(FILE *file_fp, int pbc_flag);

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void print_dashes(int ndashes,FILE *fp);

void print_neighbours( atom *p_molecule, int num_atoms,  FILE *fp);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void relabel(atom *p_molecule, int num_to_label, int *p_num_elems_used, int reset);

void minimize_molecule(atom *p_molecule, int num_atoms, char *root_name);

void minimize_inpore(atom *p_pore, atom *p_molecule,
                     int num_atoms, char *root_name,
                     double *p_kvecs, double *p_kvec2,
                     double *p_gvec2, int num_kvecs,
                     double *p_cos_sum, double *p_sin_sum,
                     int *p_need_grad, double *p_grad);

void analyse_output(char *p_analyse_file, double *p_kvecs, double *p_kvec2,
                     double *p_gvec2, int num_kvecs,
                     double *p_cos_sum, double *p_sin_sum,
                     int *p_need_grad, double *p_grad, atom *pore)
{
#include "header.h"
#define DATELINE "date"

char filename_root[FILELEN_MAX];
char discover_root[FILELEN_MAX];
char gas_min_file[FILELEN_MAX];
char sorted_file[FILELEN_MAX];
char inpore_min_file[FILELEN_MAX];

FILE *analyse_file_fp;
FILE *gas_min_fp;
FILE *inpore_min_fp;
FILE *sorted_fp;

char *fullstop;

int this_molecule;
int number_of_molecules;
int num_molecule_atoms;

char molecule_title[BUFFER];
char new_title[BUFFER];
char temp_title[BUFFER];

atom molecule[MAXTEMPLATE];
atom_number molecule_types[100];
int frame_counter;
int read_this_frame;
int molecule_index[MAXTEMPLATE]; 
							/* index of molecule energy : too many at moment */
float sort_this[MAXTEMPLATE]; 
energy gas_phase_energy[100];
energy inpore_energy[100];
energy binding_energy[100];

int num_molecule_types;
int num_elems_used[NUM_ELEMENTS];
int ineigh;

int num_temp_atoms_with_symm;

/****DEBUG*****/
/* printf("DBANAL>> Entering analyse_output of %s\n", p_analyse_file);*/
/****ENDDEBUG*****/

/******************************************************************************/
/****                  Set to use discover only                         *******/
/******************************************************************************/
strcpy(minimizer_name, DISCOVER_MINIMIZER);

/************************** derive output file names **************************/

/********************** get the root of the analyse_name   ********************/

    if (strstr(p_analyse_file, "./"))
       {
         strcpy(filename_root, strstr(p_analyse_file, "./")+2);
       }
    else
       {
         strcpy(filename_root,p_analyse_file);
       }

if ((fullstop = strchr(filename_root, '.')) != NULL)
                                        *(++fullstop) = '\0' ;
*(--fullstop) = '\0' ; /* now take off the stop */
strcat(filename_root,"_");

/******************** gas phase minimisation filename *************************/
strcpy(gas_min_file, filename_root);
strcat(gas_min_file,"gasphase.arc");
if ((gas_min_fp = fopen(gas_min_file, "w")) == NULL)
   {
   fprintf(output_fp," ZEBEDDE ERROR: Unable to open file %s for writing.\n",gas_min_file);
   exit(0);
   }


/******************** in_pore minimisation filename ***************************/
/*****in_pore  minimisation filename*****/
strcpy(inpore_min_file, filename_root);
strcat(inpore_min_file,"inpore.arc");
/* printf("DBANAL>> In pore results to be called: %s\n", inpore_min_file);*/
if ((inpore_min_fp = fopen(inpore_min_file, "w")) == NULL)
   {
   fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for writing.\n",inpore_min_file);
   exit(0);
   }

/**************** Read the arhive to be analysed to count molecules ***********/
number_of_molecules = 0;

if ((analyse_file_fp = fopen(p_analyse_file, "r")) == NULL)
   {
   fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for input.\n",p_analyse_file);
   exit(0);
   }
else
   {
   while (fgets(buffer,BUFFER,analyse_file_fp) != NULL)
  	{
  	if (strstr(buffer, DATELINE) != NULL) number_of_molecules++;
  	}
   }

fclose(analyse_file_fp);


/******************************************************************************/
/**************** Energy minimise each molecule in the gas phase **************/
/******************************************************************************/

/******************** get the root of the gas phase strategy ******************/
if (strstr(template_strategy_file,"./"))
   {
    strcpy(discover_root, strstr(template_strategy_file,"./")+2);
   }
else
   {
    strcpy(discover_root, template_strategy_file);
   }
if ((fullstop = strchr(discover_root,'.')) !=NULL) *(++fullstop) = '\0' ;
*(--fullstop) = '\0' ; /* now take off the stop */

/******************** re-open archive for reading again ***********************/
if ((analyse_file_fp = fopen(p_analyse_file, "r")) == NULL)
   {
   fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for input.\n",p_analyse_file);
   exit(0);
   }

/***throw away the top bits of the molecule.arc****/
fgets(molecule_title,BUFFER,analyse_file_fp);
fgets(dummy_head,BUFFER,analyse_file_fp);

/*** start writing results ****/
print_biosym_header(gas_min_fp, FALSE);

/**************** Loop doing Energy minimisations *****************************/
for (this_molecule=1; this_molecule <=number_of_molecules; this_molecule++)
    {
    /****** read in archive frame ***************/
    fgets(molecule_title, BUFFER,analyse_file_fp);
    fgets(dummy_head,BUFFER,analyse_file_fp);

    num_molecule_atoms=0;

    while(fgets(buffer,BUFFER,analyse_file_fp) != NULL)
      {
      if (strstr(buffer,"end") != NULL)
       { /*end of car */
       /****jettison the next end as well ******/
       fgets(buffer,BUFFER,analyse_file_fp);
       break;
       }
      else
       {
       /*read the atom*/
       sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf",
                        &(molecule[num_molecule_atoms].label[0]), 
			&molecule[num_molecule_atoms].x,
                        &molecule[num_molecule_atoms].y, 
			&molecule[num_molecule_atoms].z,
                        &(molecule[num_molecule_atoms].group[0]),
			&(molecule[num_molecule_atoms].group_no[0]),
                        &(molecule[num_molecule_atoms].pot[0]),
			&(molecule[num_molecule_atoms].elem[0]),
                        &molecule[num_molecule_atoms].part_chge);
       /* initialise neighbour information */
       molecule[num_molecule_atoms].num_neigh=-1;
       for (ineigh=0; ineigh < 4; ineigh++)
                       molecule[num_molecule_atoms].neighb[ineigh]=-1;

       num_molecule_atoms++;
       }
      }
    num_molecule_atoms--;

    /* Not sure yet if it is a good idea to relabel or not.........          */
    /* relabel( &molecule[0], num_molecule_atoms, &num_elems_used[0], TRUE); */
    /* Not sure yet if it is a good idea to relabel or not.........          */

/***************** generate neighbours for this molecule **********************/
    generate_neighbours( &molecule[0], num_molecule_atoms, &molecule_types[0],
                        &num_molecule_types, FALSE);

/***************** minimise as a molecule *************************************/
    minimize_molecule(&molecule[0], num_molecule_atoms, discover_root);
   
/**************** get the total energy and add to title line ******************/
//  gas_phase_energy[this_molecule].minimizer_init_total = 
//					interaction_energy.minimizer_init_total;
//  gas_phase_energy[this_molecule].minimizer_end_total = 
//                              interaction_energy.minimizer_end_total;
//  gas_phase_energy[this_molecule].minimizer_init_nonbond = 
//                              interaction_energy.minimizer_init_nonbond;
//  gas_phase_energy[this_molecule].minimizer_end_nonbond = 
//                              interaction_energy.minimizer_end_nonbond;

/**** this bit doesn't quite work! WHY NOT???? *************/
    strip_spaces(&molecule_title[0]);
    sprintf(new_title, "%s E(gp)= %10.6f", molecule_title, 
			gas_phase_energy[this_molecule].minimizer_end_total);
    
/***************** print out the molecule to the gas_phase archive ************/
    print_frame_header(&new_title[0], -1, gas_min_fp);
    print_molecule(&molecule[0], num_molecule_atoms, gas_min_fp, FALSE);
    fprintf(gas_min_fp,"end\n");
    fprintf(gas_min_fp,"end\n");
    
    }
fclose(gas_min_fp);


/******************************************************************************/
/**************** Energy minimise each molecule in the HOST      **************/
/******************************************************************************/

/******************** get the root of the inpore  strategy ******************/
strcpy(discover_root, inpore_strategy_file);
if (strstr(inpore_strategy_file,"./"))
   {
    strcpy(discover_root, strstr(inpore_strategy_file,"./")+2);
   }
else
   {
    strcpy(discover_root, inpore_strategy_file);
   }

if ((fullstop = strchr(discover_root,'.')) !=NULL) *(++fullstop) = '\0' ;
*(--fullstop) = '\0' ; /* now take off the stop */

/******************** re-open archive for reading again ***********************/
if ((analyse_file_fp = fopen(p_analyse_file, "r")) == NULL)
   {
   fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for input.\n",p_analyse_file);
   exit(0);
   }

/***throw away the top bits of the molecule.arc****/
fgets(molecule_title,BUFFER,analyse_file_fp);
fgets(dummy_head,BUFFER,analyse_file_fp);

/*** start writing results ****/
print_biosym_header(inpore_min_fp, FALSE);


/**************** Looping Energy minimising in the host ***********************/

for (this_molecule=1; this_molecule <=number_of_molecules; this_molecule++)
    {
    /****** read in archive frame *****/
    fgets(molecule_title, BUFFER,analyse_file_fp);
    fgets(dummy_head,BUFFER,analyse_file_fp);

    num_molecule_atoms=0;

    while(fgets(buffer,BUFFER,analyse_file_fp) != NULL)
      {
      if (strstr(buffer,"end") != NULL)
       { /*end of car */
       /****jettison the next end as well ******/
       fgets(buffer,BUFFER,analyse_file_fp);
       break;
       }
      else
       {
       /*read the atom*/
       sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf",
                        &(molecule[num_molecule_atoms].label[0]), 
			&molecule[num_molecule_atoms].x,
                        &molecule[num_molecule_atoms].y, 
			&molecule[num_molecule_atoms].z,
                        &(molecule[num_molecule_atoms].group[0]),
			&(molecule[num_molecule_atoms].group_no[0]),
                        &(molecule[num_molecule_atoms].pot[0]),
			&(molecule[num_molecule_atoms].elem[0]),
                        &molecule[num_molecule_atoms].part_chge);
       /* initialise neighbour information */
       molecule[num_molecule_atoms].num_neigh=-1;
       for (ineigh=0; ineigh < 4; ineigh++)
                       molecule[num_molecule_atoms].neighb[ineigh]=-1;

       num_molecule_atoms++;
       }
      }
    num_molecule_atoms--;

    printf("analyse_output 1: Calling relabel with num_to_label = %d\n", num_molecule_atoms);
    relabel( &molecule[0], num_molecule_atoms, &num_elems_used[0], TRUE);
    /* Not sure yet if it is a good idea to relabel or not.........          */

/***************** generate neighbours for this molecule **********************/
    generate_neighbours( &molecule[0], num_molecule_atoms, &molecule_types[0],
                        &num_molecule_types, FALSE);

/***************** Add symmetry related images if set    **********************/
if (symm_set)
   {
   do_symmetry(&molecule[0], num_molecule_atoms, num_symm_ops,
                                     &symm[0], &num_temp_atoms_with_symm);
   }

/****************** minimise inpore *******************************************/
    minimize_inpore(&pore[0], &molecule[0], num_molecule_atoms, discover_root,
                    p_kvecs, p_kvec2, p_gvec2, num_kvecs, p_cos_sum, p_sin_sum,
                    p_need_grad, p_grad);
   
    /******* get the total energy and make a title line out of it *********/
//  inpore_energy[this_molecule].minimizer_init_total = 
//                                              interaction_energy.minimizer_init_total;
//  inpore_energy[this_molecule].minimizer_end_total = 
//                                              interaction_energy.minimizer_end_total;
//  inpore_energy[this_molecule].minimizer_init_nonbond = 
//                                              interaction_energy.minimizer_init_nonbond;
//  inpore_energy[this_molecule].minimizer_end_nonbond = 
//                                            	interaction_energy.minimizer_end_nonbond;
    
    /****normalise if symmetry related images****/
    if (symm_set)
     {
      inpore_energy[this_molecule].minimizer_init_total =
           inpore_energy[this_molecule].minimizer_init_total / num_symm_ops+1;
      inpore_energy[this_molecule].minimizer_end_total =
           inpore_energy[this_molecule].minimizer_end_total / num_symm_ops+1;
      inpore_energy[this_molecule].minimizer_init_nonbond =
           inpore_energy[this_molecule].minimizer_init_nonbond / num_symm_ops+1;
      inpore_energy[this_molecule].minimizer_end_nonbond =
           inpore_energy[this_molecule].minimizer_end_nonbond / num_symm_ops+1;

     }

    strip_spaces(&molecule_title[0]);
    sprintf(new_title, "%s E(inhost)= %10.6f", molecule_title, 
			inpore_energy[this_molecule].minimizer_end_total);
    
/***************** print out the molecule to the in_pore archive **************/
    print_frame_header(&new_title[0], -1, inpore_min_fp);
    print_molecule(&molecule[0], num_molecule_atoms, inpore_min_fp, FALSE);
    fprintf(inpore_min_fp,"end\n");
    fprintf(inpore_min_fp,"end\n");
    
    }
fclose(inpore_min_fp);

/**************** Now calculate the binding energies **************************/
for (this_molecule=1; this_molecule <=number_of_molecules; this_molecule++)
  {
  binding_energy[this_molecule].minimizer_init_total =
	inpore_energy[this_molecule].minimizer_init_total -
		gas_phase_energy[this_molecule].minimizer_init_total;
  binding_energy[this_molecule].minimizer_end_total =
	inpore_energy[this_molecule].minimizer_end_total -
		gas_phase_energy[this_molecule].minimizer_end_total;
  binding_energy[this_molecule].minimizer_init_nonbond =
	inpore_energy[this_molecule].minimizer_end_total -
		gas_phase_energy[this_molecule].minimizer_end_total;
  binding_energy[this_molecule].minimizer_end_nonbond =
		inpore_energy[this_molecule].minimizer_end_nonbond - 
			gas_phase_energy[this_molecule].minimizer_end_nonbond;
  }


/**************** Sort them ***************************************************/
/** molecule_index is a blank array of integers which will be filled by the  **/
/** sort routine. It is then sorted accoring to the passed array and so will **/
/** come back with the order of molecules                                    **/

/*** doesn't seem to like being passed a structure so copy to new array ***/
for (this_molecule=1; this_molecule <=number_of_molecules; this_molecule++)
  {
  sort_this[this_molecule]=binding_energy[this_molecule].minimizer_end_total;
  /* printf("sort this %d %10.4f\n", this_molecule, sort_this[this_molecule]);*/
  }

index_sort(number_of_molecules, &sort_this[0],
				&molecule_index[0]);

/**** DEBUG ****/
/* for (this_molecule=1; this_molecule<=number_of_molecules; this_molecule++) */
/*  { */
/*      printf("%i %i %10.4f\n", this_molecule, molecule_index[this_molecule], */
/*            binding_energy[molecule_index[this_molecule]].minimizer_end_total); */
/*  } */
/**** DEBUG ****/



/******************** Write out sorted archive ********************************/
strcpy(sorted_file, filename_root);
strcat(sorted_file, "sort.arc");

if ((sorted_fp = fopen(sorted_file, "w")) == NULL)
   {
   fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for writing.\n",sorted_file);
   exit(0);
   }

print_biosym_header(sorted_fp, FALSE);

/* get the filename of the minimised in_pore molecules */
strcpy(inpore_min_file, filename_root);
strcat(inpore_min_file,"inpore.arc");

/* loop over all the template indexes */
for (this_molecule=1; this_molecule<= number_of_molecules; this_molecule++)
    {
	if ((inpore_min_fp = fopen(inpore_min_file, "r")) == NULL)
   		{
   		fprintf(output_fp, "ZEBEDDE ERROR: Unable to open file %s for reading\n",inpore_min_file);
   		exit(0);
   		}
	/*** which frame do I want? ***/
	read_this_frame = 1;
	while (this_molecule != molecule_index[read_this_frame]) read_this_frame++;

        /**** Wind to frame frame_counter in the archive ****/
	frame_counter = 0;
	while (frame_counter != read_this_frame)
		{
		if (fgets(buffer,BUFFER,inpore_min_fp) != NULL)
		   {
			if (strstr(buffer, DATELINE) != NULL) 
				{
				frame_counter++;
				}
			}
		else
			{
			fprintf(output_fp, "SERIOUS ZEBEDDE ERROR: End of archive file in analyse_output\n");
			fprintf(output_fp, "Report to Authors\n");
			exit(0);
			}
		}


	num_molecule_atoms=0;
	while(fgets(buffer,BUFFER,inpore_min_fp) != NULL)
      {
      if (strstr(buffer,"end") != NULL)
         { /*end of car */
          break;
         }
      else
         {
         /*read the atom*/
         sscanf(buffer,"%s%lf%lf%lf%s%s%s%s%lf",
                        &(molecule[num_molecule_atoms].label[0]), 
                        &molecule[num_molecule_atoms].x,
                        &molecule[num_molecule_atoms].y, 
                        &molecule[num_molecule_atoms].z,
                        &(molecule[num_molecule_atoms].group[0]),
                        &(molecule[num_molecule_atoms].group_no[0]),
                        &(molecule[num_molecule_atoms].pot[0]),
                        &(molecule[num_molecule_atoms].elem[0]),
                        &molecule[num_molecule_atoms].part_chge);
		  num_molecule_atoms++;
		  }
	   } /*end while fgets*/
	num_molecule_atoms--;

	/*** found this frame so close the file again ****/
	fclose(inpore_min_fp);

	/*** make up a title line ****/

    sprintf(temp_title, 
		"Previous frame number %d E(gp) = %8.4f E(h) = %8.4f E(BE) = %8.4f",
		molecule_index[this_molecule],
		gas_phase_energy[molecule_index[this_molecule]].minimizer_end_total,
		inpore_energy[molecule_index[this_molecule]].minimizer_end_total,
		binding_energy[molecule_index[this_molecule]].minimizer_end_total);

	/*** write to sorted archive ***/
	print_frame_header(&temp_title[0], -1, sorted_fp);
	print_molecule(&molecule[0], num_molecule_atoms, sorted_fp, FALSE);
	fprintf(sorted_fp,"end\n");
	fprintf(sorted_fp,"end\n");

	}
fclose(sorted_fp);


/***************** Write out a nice report! ***********************************/
print_dashes(80,output_fp);
fprintf(output_fp,"ANALYSIS REPORT\n\n");
fprintf(output_fp,"Output file %s has been energy minimised\n\n",
							 p_analyse_file);
fprintf(output_fp,"The following archive files have been created:\n\n");
fprintf(output_fp,"Gas phase minimised structure: %s_gasphase.arc\n", 
								filename_root);
fprintf(output_fp,"In host  minimised structures: %s_inpore.arc\n", 
								filename_root);
fprintf(output_fp,"\n\nNote: These are in the same order as in the original file\n");

fprintf(output_fp,"\nSorted by Binding Energy     : %s_sort.arc\n\n\n", 
								filename_root);

print_dashes(80,output_fp);
fprintf(output_fp,"Frame#\tOrig#\t\tE(gas phase)\tE(in_host)\tE(bind)\n");
print_dashes(80,output_fp);
for (this_molecule=1; this_molecule <=number_of_molecules; this_molecule++)
    {
    fprintf(output_fp, "%-5d\t\t%-5d\t\t%10.3f\t%10.3f\t%10.3f\n", 
		this_molecule, molecule_index[this_molecule],
		gas_phase_energy[molecule_index[this_molecule]].minimizer_end_total,
		inpore_energy[molecule_index[this_molecule]].minimizer_end_total,
		binding_energy[molecule_index[this_molecule]].minimizer_end_total);
    }
print_dashes(80,output_fp);
if (symm_set)
  {
    fprintf(output_fp, "NOTE: There are %d symmetry related templates\n", num_symm_ops+1);
    fprintf(output_fp, "      E(in_host) is normalised to give enegy per template\n");

  }
return;
}

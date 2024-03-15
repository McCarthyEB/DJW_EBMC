/*****************************************************************************/
/* c2disco_minimize_molecule.c : minimizes a molecule with respect to itself */
/* version adapted to run c2 version of discover                             */
/*                                                                           */
/* August 98 DJW                                                             */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void regroup(atom *p_molecule, int num_atoms);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void print_biosym_header(FILE *file_fp,int pbc_flag);

void print_pbc_header(FILE *file_fp, double *p_abc);

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void car_end(FILE *fp);

void delete_file(char *fileroot, char *extension);

int anal_c2disco(char *root_file, atom *p_molecule,int num_atoms,energy *p_energy);

int  get_c2disco_minimized(char *root_file,int num_molecule, 
						atom *p_molecule, int num_atoms);

void print_mdf(char *p_molname, atom *p_molecule, int num_atoms,
                int start, FILE *fp, int D_to_H, int have_subset);

void print_mdf_subset(char *p_molname, atom *p_molecule, int num_atoms,
                      int start, FILE *fp, int D_to_H, int have_subset);

void print_mdf_header(char *p_assembly, FILE *fp);

void print_mdf_pbc(FILE *fp);

void print_mdf_end(FILE *fp);

void print_neighbours( atom *p_molecule, int num_atoms,  FILE *fp);

void c2disco_minimize_inpore(atom *p_pore, atom *p_molecule, int num_atoms, char *root_name,
                             int which_mol)
{
#include "header.h"

char fileroot[FILELEN_MAX];
char disc_new[FILELEN_MAX];
char disc_stand[FILELEN_MAX];
FILE *disc_new_fp;
FILE *disc_stand_fp;
char *fullstop;
char dollar1[BUFFER];
char dummy_name[BUFFER];
char fix_label[5] = "AAAA";
char discover_cmdline[BUFFER];
char setup_cmdline[BUFFER];
int discover_status;
int iatom, ifix, fix_this, D_to_H;
int ineigh, neighs_expected, used_neigh[10];
atom *p_atom;
atom *p_neigh;
int *p_fixed_atoms, *p_this_fix;
int done, num_to_fix, ref_neigh, icheck, checked;
int iloop;

/*************write out discover .car file***********************************/
p_fixed_atoms=malloc(num_atoms*sizeof(int));

strcpy(fileroot,root_name);
strcat(fileroot,"_tmp.car");
if (!(discover_fp = fopen(fileroot, "w")))
    {
      fprintf(output_fp, "ZEBEDDE ERROR: Error opening Discover car file %s in disco_minimize_molecule\n", fileroot);
      fprintf(output_fp, "Abandoning Minimisation\n");
      return;
    }

D_to_H = TRUE;
print_biosym_header(discover_fp, pbc); 
print_frame_header("molecule_minimize",1, discover_fp);

if (pbc)
   {
     print_pbc_header(discover_fp, &abc[0]);

     print_molecule(p_pore, num_pore_atoms, discover_fp, FALSE );
     fprintf(discover_fp,"end\n");
     print_molecule(p_molecule, num_atoms, discover_fp, D_to_H);

/******************************************************************/
/********* write the template images ******************************/
/******************************************************************/

        p_atom= p_molecule;
        for (iloop=0; iloop<=num_symm_ops; iloop++)
          {
            p_atom += (num_atoms+1);
            fprintf(discover_fp,"end\n");
            print_molecule(p_atom, num_atoms, discover_fp, D_to_H);
          }
    }
else
   {
      print_molecule(p_pore, num_pore_atoms, discover_fp, FALSE );
      fprintf(discover_fp,"end\n");
      print_molecule(p_molecule, num_atoms, discover_fp, D_to_H);
   }

car_end(discover_fp);

strcpy(fileroot,root_name);
strcat(fileroot,"_tmp.mdf");

 if (!(discover_fp = fopen(fileroot, "w")))
     {
     fprintf(output_fp,
          "ZEBEDDE ERROR: Error opening temporary Discover mdf file %s in disco_minimize_molecule\n", fileroot);
     fprintf(output_fp,"Abandoning Minimisation\n");
     return;
     }

 print_mdf_header("pore_template",discover_fp);

 print_mdf("pore", p_pore, num_pore_atoms, 0, discover_fp, D_to_H, FALSE); 

/****************************************************************************/
/***** Write the template images                                     ********/
/***** just write out the same one cos they're all the same          ********/
/****************************************************************************/

      if (pbc)
        {
          for (iloop=0;iloop<=num_symm_ops+1;iloop++)
           {
              sprintf(dummy_name,"template%d",iloop);
              print_mdf(dummy_name, p_molecule, num_atoms, 0, discover_fp, D_to_H, FALSE);
           }
    
          print_mdf_pbc(discover_fp);
       }
      else
       {
          sprintf(dummy_name,"template");
          print_mdf(dummy_name, p_molecule, num_atoms, 0, discover_fp, D_to_H, FALSE);
       }

strcpy(fileroot,root_name);
if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

/******NOW RUN IN DISCOVER AND RETRIEVE OUTPUT *****/
/***** call Discover script to runjob now at full on priority ****/
/****** syntax: discover file_name discover_forcefield_name nice start_yn      ****/

sprintf(disc_new,"./%s_tmp.run",fileroot);
sprintf(disc_stand,"./%s.run",fileroot);

if (!(disc_stand_fp = fopen(disc_stand,"r")))
  {
/* Standard run file does not exist report and abandon */

     fprintf(output_fp, "ERROR >> Standard c2 discover run file: %s",disc_stand);
     fprintf(output_fp, " not supplied create this file or amend dat file\n");
     exit(0);
   }

/*******************************************************/
/*** Adapt the standard run file to our needs **********/
/*******************************************************/

disc_new_fp = fopen(disc_new,"w");

done=FALSE;
while (fgets(buffer,BUFFER, disc_stand_fp) != NULL)
  {
    sscanf(buffer,"%s", dollar1);

/******* Retain the setenvs ****************************/

    if (strcmp(dollar1, "setenv") == 0 || strcmp(dollar1, "#!/bin/csh") == 0 || strcmp(dollar1, "source") == 0 )
      {
        fputs(buffer, disc_new_fp);
      }    

/******* Substitute our molecule name for the discovery run part *****/
    
    if (strstr(dollar1, "discovery") != NULL && !done)
      {
        done=TRUE;
        fprintf(disc_new_fp, "%s %s_tmp\n", dollar1, root_name);  
      }
  }

fclose(disc_stand_fp);
fclose(disc_new_fp);

/******************************************************************/
/****Sort out the Discover .inp file                           ****/
/****Fix the pore atoms and                                    ****/
/****Fix the atoms requested by the user, detect by neighbours ****/
/******************************************************************/

num_to_fix=-1;
for (ifix=0; ifix <=num_mini_fix; ifix++)
  {
    p_atom = p_molecule-1;
    neighs_expected= fix_atoms[ifix].num_neigh;

    for (iatom=0; iatom <= num_atoms; iatom++)
      {
         p_atom++;
/******************************************************************/
/***** Check for Central atom element types ***********************/
/******************************************************************/

         if (strcmp(p_atom->elem, &(fix_atoms[ifix].central_elem[0])) == 0 &&
                                                    p_atom->num_neigh == neighs_expected )
           {
/******************************************************************/
/***** Check we have the right neighbours *************************/
/******************************************************************/

             for (ref_neigh=0; ref_neigh <= neighs_expected; ref_neigh++) used_neigh[ref_neigh]= FALSE;

             for (ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
               {
                  p_neigh= p_molecule+ p_atom->neighb[ineigh];

                  done= FALSE;
                  ref_neigh= 0;
                  while (ref_neigh <= neighs_expected && !done)
                    {
                      if (strcmp(p_neigh->elem, &(fix_atoms[ifix].neigh[ref_neigh][0]) ) == 0 &&
                                                                           !used_neigh[ref_neigh] )
                        {
                           used_neigh[ref_neigh]= TRUE;
                           done= TRUE;
                        }
                      ref_neigh++;
                    }
               }

/**************************************************************************/
/***** Fix those atoms with the right neighbours and the neighbours too ***/
/**************************************************************************/

             fix_this= TRUE;
/**************************************************************************/
/***** For the inpore case let everything rip  DJW Cludge June 99 *********/
/**************************************************************************/
             fix_this= FALSE;
/**************************************************************************/
             for (ref_neigh=0; ref_neigh <= neighs_expected; ref_neigh++)
               {
                 if (!used_neigh[ref_neigh]) fix_this= FALSE;
               }

             if (fix_this) 
               {
                 checked=TRUE;
                 p_this_fix=p_fixed_atoms;
                 for (icheck=0; icheck <= num_to_fix; icheck++) 
                   {
                      if (*p_this_fix == iatom) checked=FALSE;    
                      p_this_fix++;
                   }
  
                 if (checked) 
                   {
                      num_to_fix++;
                      *(p_fixed_atoms+num_to_fix)= iatom;
                   }
               
                 for (ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
                   {
                      checked=TRUE;
                      p_this_fix=p_fixed_atoms;
                      for (icheck=0; icheck <= num_to_fix; icheck++)  
                        {
                          if (*p_this_fix == p_atom->neighb[ineigh]) checked=FALSE;
                          p_this_fix++;
                        }
 
                      if (checked)
                        {  
                           num_to_fix++;
                           *(p_fixed_atoms+num_to_fix)= p_atom->neighb[ineigh]; 
                        }
                   }
               }
           }
      }  
  }

/*************write out discover .mdf file***********************************/

sprintf(disc_stand,"./%s.inp",fileroot);

if (!(disc_stand_fp = fopen(disc_stand,"r")))
  {
/* Standard run file does not exist report and abandon */

     fprintf(output_fp, "ERROR >> Standard c2 discover inp file for stem: %s",disc_stand);
     fprintf(output_fp, " not supplied create this file or amend dat file\n");
     exit(0);
   }

strcpy(fileroot,root_name);
strcat(fileroot,"_tmp.inp");

if (!(disc_new_fp = fopen(fileroot, "w")))
    {
    fprintf(output_fp,
         "ZEBEDDE ERROR: Error opening temporary Discover inp file %s in disco_minimize_molecule\n", fileroot);
    fprintf(output_fp,"Abandoning Minimisation\n");
    return;
    }

done=FALSE;
while (fgets(buffer,BUFFER, disc_stand_fp) != NULL)
  {
    sscanf(buffer,"%s", dollar1);

/******* Look for begin line in input file *************/

    if (strcmp(dollar1, "begin") != 0 || done)
      {
        fputs(buffer, disc_new_fp);
      }   
    else
      {
        done=TRUE;
        fputs(buffer, disc_new_fp);

/****** Fix Pore atoms using a subset ********************/

print_mdf_subset("pore", p_pore, num_pore_atoms, 0, discover_fp, D_to_H, TRUE);

        if (num_to_fix >= 0)
          {
             fprintf(disc_new_fp, "#\n# Fix constrained Atom subset added by zebedde\n" );
             p_this_fix=p_fixed_atoms;
             for (ifix=0; ifix <=num_to_fix; ifix++)
               {
                  fix_label[0]= 'T'; 
                  fix_label[1]= 'A'+ifix/26;
                  fix_label[2]= 'A'+ifix-26*(ifix/26);

                  p_atom= p_molecule+*p_this_fix;
                  fprintf(disc_new_fp, "  setAtomMovability  fixed %s *:*:%s\n", &fix_label[0], p_atom->label); 
                  p_this_fix++;
               }
          }
      }
  }

free(p_fixed_atoms);
fclose(discover_fp);
fclose(disc_stand_fp);
fclose(disc_new_fp);

/**************************************************************/
/***** Run Discover *******************************************/
/**************************************************************/

sprintf(setup_cmdline, "chmod 755 %s_tmp.run", root_name);
system (setup_cmdline);	

sprintf(discover_cmdline, "%s_tmp.run", root_name);
system (discover_cmdline);	

/*******get results from DISCOVER run *************************/

strcpy(fileroot,root_name);
strcat(fileroot,"_tmp");

if ((fullstop = strchr(fileroot,'.')) !=NULL) *fullstop = '\0';

discover_status = anal_c2disco(fileroot, p_molecule,
                                                num_atoms, &interaction_energy[which_mol]);

if (discover_status == 0)
  {

/********************************************************/
/****** Read back structure if energy lowered ***********/
/****** Because this is an inpore minimisation **********/
/****** the new template will be molecule 1    **********/
/****** i.e. the second molecule in the car file ********/
/********************************************************/

     printf("DEBUG>> Energies before and after: %10.6f %10.6f \n",
                                        interaction_energy[which_mol].minimizer_init_total,
                                        interaction_energy[which_mol].minimizer_end_total);

     if ( interaction_energy[which_mol].minimizer_end_total 
                           < interaction_energy[which_mol].minimizer_init_total )
       {
          get_c2disco_minimized( fileroot, 1, p_molecule, num_atoms);

/********************************************************/
/****** Protect the D atoms *****************************/
/********************************************************/

          p_atom=p_molecule;
          for (iatom=0; iatom <= num_atoms; iatom++)
            {
              if (strstr(p_atom->label, "D") != NULL && 
                  strcmp(p_atom->elem, "H") == 0)        strcpy(p_atom->elem, "D");
              p_atom++;
            }
       }
  }

/***Dont ! tidy up a bit to minimize disk usage***/
/*
delete_file(root_name,"_tmp.prm");
delete_file(root_name,"_tmp.out");
delete_file(root_name,"_tmp.cor");
*/

return;
}

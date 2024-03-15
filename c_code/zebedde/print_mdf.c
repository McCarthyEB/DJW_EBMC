/****************************************************************/
/* print_mdf.c : prints a mdf for a molecule                    */
/*               Biosym archive version 3                       */
/*three routines, print_mdf which lobs out the main part of the */
/*                mdf file                                      */
/*                print_mdf_header which lobs out the preamble  */
/*            and print_mdf_end which finishes it off           */
/*                                                              */
/* started Dewi 2/5/95                                          */
/* 12/6 SERIOUS MODS! DWL                                       */
/* Now works Properly and does Discover 94 mdfs                 */
/* and does pbc if required from abc[] and pbc flags            */
/*                                                              */
/* NB. Still need to sort oop for benzenes                      */
/* Adapted for D_to_H DJW August 98                             */
/* 
 * AJWL Dec 08  Added routine to print connectivities 		
 * into MDF files.     
 * AJWL Feb 08 Major mods. Removed subset printing till later
 * Finally it works!! 						*/
/****************************************************************/


#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void print_mdf(char *p_molname, atom *p_molecule, int num_atoms, 
		int start, FILE *fp, int D_to_H, int have_subset)
{
#include "header.h"
int iatom1,iatom2;
int idummy1;
int is_pore = FALSE;
char dummy[40];
char dummy2[40];
char elem[6];
char this_group[10], that_group[10];
double shift;
atom *p_atom1, *p_atom2;
int xshift, yshift, zshift;
fprintf(fp, "@molecule %s\n\n",p_molname);

if (strcmp(p_molname,"pore") == 0)
{
is_pore = TRUE;
}
/* loop over pairs */

if (D_to_H)
  {
    p_atom1= p_molecule-1;
    for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
      {
          p_atom1++;

/*** Do D to H conversion *****/
          strcpy(elem, p_atom1->elem);
          if (strcmp(elem,"D")==0) strcpy(elem, "H");

/*** Print out details of this atom *****/
		  sprintf(dummy,"%s_%s:%s",     p_atom1->group,
                                      p_atom1->group_no,p_atom1->label);
		  fprintf(fp,"%-20s",dummy);
          fprintf(fp, "%-2s %-3s",elem,p_atom1->pot);
          fprintf(fp, "     %-4s",p_atom1->group);
          fprintf(fp, " 0  0");
          fprintf(fp, "  %7.4f",p_atom1->part_chge);
          fprintf(fp, " 0 0 0");
          fprintf(fp, "  1.0000  0.0000");

/*  new mdf needs this according to manual but it doesn't!! */
/*  fprintf(fp, " %i ",p_atom1->num_neigh); */

          sprintf(dummy,"%s_%s:", p_atom1->group, p_atom1->group_no);

/*** Print neighbours for connectivity ****/
if (is_pore == FALSE)
{
          for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
            {
              iatom2= (p_atom1->neighb[idummy1]) - start ;
              p_atom2= p_molecule+iatom2;

              sprintf(dummy2,"%s_%s:", p_atom2->group, p_atom2->group_no);
              if (strcmp(dummy,dummy2) == 0)
                {
                  strcpy(dummy2,"");
                }
              fprintf(fp, " %s%s",dummy2, p_atom2->label);
           }
}
/*** Print neighbours for connectivity ****/
if (is_pore == TRUE)
          {
	for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
            {
              iatom2= (p_atom1->neighb[idummy1]) - start ;
              p_atom2= p_molecule+iatom2;

              sprintf(dummy2,"%s_%s:", p_atom2->group, p_atom2->group_no);
              if (strcmp(dummy,dummy2) == 0)
                {
                  strcpy(dummy2,"");
                }
	      /******CALC IF IN ANOTHER BOX - kludge for now********/
              xshift = 0;
              yshift = 0;
              zshift = 0;
              if((p_atom1->x/p_atom2->x) >=0.0);
              {
              /*in same box*/
              shift = (p_atom1->x - p_atom2->x);
              if (shift> abc[0]/2.0) xshift = 1;
              if (shift< -1.0*abc[0]/2.0) xshift = -1;
              }
            /*  printf("DB>>xshift %f %f %f %f %d\n", p_atom1->x , p_atom2->x,shift, abc[0]/2.0, xshift);*/
              shift = (p_atom1->y - p_atom2->y);
              if (shift> abc[1]/2.0) yshift = 1;
              if (shift< -1.0*abc[1]/2.0) yshift = -1;
/*              printf("DB>>yshift %f %f %f %f %d\n", p_atom1->y , p_atom2->y,shift, abc[0]/2.0, yshift);*/
              shift = (p_atom1->z - p_atom2->z);
              if (shift> abc[2]/2.0) zshift = 1;
              if (shift< -1.0*abc[2]/2.0) zshift = -1;
  /*            printf("DB>>zshift %f %f %f %f %d\n", p_atom1->z , p_atom2->z,shift, abc[0]/2.0, zshift);*/
              if ((xshift==0) && (yshift==0) && (zshift==0))
                 {
                  /*nothing to add*/
                  fprintf(fp, " %s%s",dummy2, p_atom2->label);
                 }
              else
                 {
                  fprintf(fp, " %s%s%%%d%d%d#1",dummy2, p_atom2->label,
                                                xshift,yshift,zshift);
                 }

}
}


/***** The extra space is required for the last neighbour to be recoginised by cerius DJW March 99 ***/
      fprintf(fp," \n");
     }
  }
else
  {
    for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
      {
          p_atom1= p_molecule+iatom1;

          sprintf(dummy,"%s_%s:%s",	p_atom1->group,
							p_atom1->group_no,p_atom1->label);
          fprintf(fp,"%-20s",dummy);
          fprintf(fp, "%-2s %-3s",p_atom1->elem,p_atom1->pot);
          fprintf(fp, "     %-4s",p_atom1->group);
          fprintf(fp, " 0  0");
          fprintf(fp, "  %7.4f",p_atom1->part_chge);
          fprintf(fp, " 0 0 0");
          fprintf(fp, "  1.0000  0.0000");
/*  new mdf needs this according to manual but it doesn't!! */
/*  fprintf(fp, " %i ",p_atom1->num_neigh); */

          sprintf(dummy,"%s_%s:", p_atom1->group, p_atom1->group_no);
          for (idummy1=0; idummy1 <= (p_atom1->num_neigh); idummy1++)
            {
              iatom2= (p_atom1->neighb[idummy1]) - start ;
	            	sprintf(dummy2,"%s_%s:", (p_molecule+iatom2)->group,
                                                    (p_molecule+iatom2)->group_no);
              if (strcmp(dummy,dummy2) == 0)
		{
                  strcpy(dummy2,"");
		}
              fprintf(fp, " %s%s",dummy2, (p_molecule+iatom2)->label);
           }
/***** The extra space is required for the last neighbour to be recoginised by cerius DJW March 99 ***/
      fprintf(fp," \n");
     }
  }


//if (anneal_run) have_subset = FALSE;
if (is_pore == TRUE) have_subset = TRUE;
if (!have_subset)  fprintf(fp, "\n\n");

/****** If required print out atom subset for fixing ********/
/****** Also need to put on ending in this case      ********/
/*
if (have_subset)
  {
	fprintf(fp, "\n#atomset                                                                     \n\n\n"); 
    fprintf(fp, "@list subset CERIUS2_FIXED_ATOMS\n");

    p_atom1= p_molecule;
    sprintf(this_group, "%s_%s", p_atom1->group, p_atom1->group_no);
    fprintf(fp, "pore:%s:", this_group);

    iatom2=0;
    for (iatom1=0; iatom1<= num_atoms; iatom1++)
      {
        sprintf(that_group, "%s_%s", p_atom1->group, p_atom1->group_no);
        if (strcmp(this_group, that_group) == 0 )
          {
            iatom2++;
            if (iatom2 < 17)
              {
                 fprintf(fp, "%s ",p_atom1->label); 
              }
            else
              {
                 fprintf(fp, "\n%s ",p_atom1->label); 
                 iatom2=0;
              }
          }
        else
          {
            iatom2=0;
            strcpy( this_group, that_group );
            fprintf(fp, "\n%s:%s ",this_group, p_atom1->label);
          }
        p_atom1++;
      }

    fprintf(fp, "\n#end\n");
  }
printf("Have closed fp\n");
*/
return;
}

void print_mdf_header(char *p_assembly, FILE *fp)
{
time_t idate;

fprintf(fp, "!BIOSYM molecular_data 4\n\n");
idate = time(NULL);

fprintf(fp, "!DATE:      %24.24s",ctime(&idate));
fprintf(fp, "     INSIGHT generated molecular data file\n\n");
fprintf(fp, "#topology\n\n");
fprintf(fp, "@column 1 element\n");
fprintf(fp, "@column 2 atom_type cff91\n");
fprintf(fp, "@column 3 charge_group cff91\n");
fprintf(fp, "@column 4 isotope\n");
fprintf(fp, "@column 5 formal_charge\n");
fprintf(fp, "@column 6 charge cff91\n");
fprintf(fp, "@column 7 switching_atom cff91\n");
fprintf(fp, "@column 8 oop_flag cff91\n");
fprintf(fp, "@column 9 chirality_flag\n");
fprintf(fp, "@column 10 occupancy\n");
fprintf(fp, "@column 11 xray_temp_factor\n");
fprintf(fp, "@column 12 connections\n\n");

if (p_assembly != NULL)  /* writing assembly not a molecule */
	{
	fprintf(fp, "@assembly %s\n\n\n",p_assembly);
	}

return;
}

void print_mdf_pbc(FILE *fp)
{
fprintf(fp, "!\n");
fprintf(fp, "#symmetry\n");
fprintf(fp, "@periodicity 3 xyz\n");
fprintf(fp, "@group (P1)\n");

return;
}
void print_mdf_end(FILE *fp)
{
fprintf(fp, "#end\n");
fclose(fp);
return;
}

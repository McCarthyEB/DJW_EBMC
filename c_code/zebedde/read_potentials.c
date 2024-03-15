/****************************************************/
/********* Routines for reading .frc files      *****/
/********* begun 27/11/95 Dave Willock          *****/
/*********                                      *****/
/********* Updated for Oie potentials with      *****/
/********* possibility of multiple non_bond     *****/
/********* potential types                      *****/
/********* Dave Willock, Jan 06                 *****/
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "data.h"

int open_file(FILE **p_file, char *p_filename, char *p_status);

int read_line(FILE *fp, int *p_ichar);

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars );
 
int string_from_int(int *p_int, char *p_string);

void get_nonbond_info(FILE *fp, 
                      int *p_line, char *p_line_string,
                      int *p_num_of_chars);

void get_nonbond_params(FILE *fp, int *p_line, char *p_line_string,
                        int *p_num_of_chars, int have_comb_rules, int type);

void get_hbond_params(FILE *fp, int *p_line, char *p_line_string,
                        int *p_num_of_chars);

void get_stretch_params(FILE *fp, int *p_line, char *p_line_string,
                                int *p_num_of_chars, int which);

void get_equivalences(FILE *fp, int *p_line, char *p_line_string,
                                int *p_num_of_chars);

void get_angle_params(FILE *fp, int *p_line, char *p_line_string,
                        int *p_num_of_chars, int which);

void get_torsion_params(FILE *fp, int *p_line, char *p_line_string, 
                        int *p_num_of_chars, int which);

void get_bond_incs(FILE *fp, int *p_line, char *p_line_string,
                                                 int *p_num_of_chars);

int read_potentials(int *p_have_comb_rules)
{

#include "header.h"

int num_of_chars,iloop;
int line[BUFFER],num_lines;
int idummy, version_found;
int which_stretch, found_stretch;
int which_angle, found_angle;
int which_torsion, found_torsion;
int found_vdw, this_type;
int worked;

double a2;

char line_string[BUFFER];
char h_bond[BUFFER], equivalence[BUFFER];
char bond_inc[BUFFER];

FILE *file_fp;

/**************************************************/
/** Set up counters for multi-entry lists *********/
/**************************************************/
num_potential_types=-1;
/**************************************************/
worked = open_file(&file_fp, &forcefield_library[0] , "r");

if (worked == EXIT_FAILURE) return EXIT_FAILURE;

/**************************************************/
/*** Work out what type of forcefield it is *******/
/**************************************************/

num_of_chars=1;
num_lines=0;
version_found=FALSE;

while (num_of_chars != END_OF_FILE && !version_found )
  {
    num_lines++;
    num_of_chars=read_line(file_fp, &line[0]);
    idummy= string_from_int(&line[0], &line_string[0]);

    if ( locate_string_in_string(VERSION, &line_string[0], num_of_chars) )
       {
          if ( locate_string_in_string(PCFF_STRING, &line_string[0], num_of_chars) )  
             {
                fprintf(output_fp,"The potential file contains terms for the pcff potential type.\n");
                pot_info.version= PCFF;
                version_found=TRUE;
                strcpy(&equivalence[0], EQUIVALENCE_PCFF);
                strcpy(&bond_inc[0], BOND_INCREMENTS);
             }
          else if ( locate_string_in_string(CVFF_STRING, &line_string[0], num_of_chars) ) 
             {
                fprintf(output_fp,"The potential file contains terms for the cvff potential type.\n");
                pot_info.version= CVFF;
                version_found=TRUE;
                strcpy(&equivalence[0], EQUIVALENCE_CVFF);
                strcpy(&bond_inc[0], BOND_INCREMENTS);
             } 
          else if ( locate_string_in_string(CFF91_STRING, &line_string[0], num_of_chars) )
             {
                fprintf(output_fp,"The potential file contains terms for the cff91 potential type.\n");
                pot_info.version= CFF91;
                version_found=TRUE;
                strcpy(&equivalence[0], EQUIVALENCE_CFF91);
                strcpy(&bond_inc[0], BOND_INCREMENTS);
             }
          else if ( locate_string_in_string(AMBER_STRING, &line_string[0], num_of_chars) )
             {
                fprintf(output_fp,"The potential file contains terms for the AMBER potential type.\n");
                pot_info.version= AMBER;
                fprintf(output_fp,"This is potential version %d within zebedde.\n",pot_info.version);
                version_found=TRUE;
                strcpy(&h_bond[0], H_BOND_AMBER);
                strcpy(&equivalence[0], EQUIVALENCE_AMBER);
                strcpy(&bond_inc[0], BOND_INCREMENTS);
             }
          else if ( locate_string_in_string(OIE_STRING, &line_string[0], num_of_chars) )
             {
                fprintf(output_fp,"The potential file contains terms for the OIE potential type.\n");
                pot_info.version= OIE;
                fprintf(output_fp,"This is potential version %d within zebedde.\n",pot_info.version);
                version_found=TRUE;
                strcpy(&bond_inc[0], BOND_INCREMENTS);
             }
          else
             {
                fprintf(output_fp,"ERROR: Cannot interpret potential type from frc file.\n");
                fprintf(output_fp,"Offending line : %s\n", line_string);
                fflush(output_fp);
                fflush(stdout);
             }
       }
  }

/*************************************************************************/
/***** Get the potential information *************************************/
/*************************************************************************/
num_stretches=-1;
num_angles=-1;

while (num_of_chars != END_OF_FILE )
  {
    num_lines++;
    num_of_chars=read_line(file_fp, &line[0]);
    idummy= string_from_int(&line[0], &line_string[0]);

/**************************************************/
/******* Now try and read info from these *********/
/******* lines !!!                        *********/
/**************************************************/
/******* First see if we have one of the  *********/
/******* Myriad of possible stretch pots  *********/
/**************************************************/

    found_stretch= FALSE;
    which_stretch= FAILED_STRETCH;
     if ( locate_string_in_string(MORSE_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= MORSE_STRETCH;
        }
     else if (locate_string_in_string(QUARTIC_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= QUARTIC_STRETCH;
        }
     else if (locate_string_in_string(QUADRATIC_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= QUADRATIC_STRETCH;
        }

/**************************************************/
/** recognise angle bending potentials   **********/
/** and decide which one is present      **********/
/** Initially working with pcff frc file **********/
/**************************************************/
    found_angle  = FALSE;
    which_angle = FAILED_ANGLE;
    if (locate_string_in_string(QUARTIC_ANGLE_STRING, &line_string[0], num_of_chars))
      {
        found_angle = TRUE;
        which_angle = QUARTIC_ANGLE;
      }
    else if (locate_string_in_string(QUADRATIC_ANGLE_STRING, &line_string[0], num_of_chars))
      {
        found_angle = TRUE;
        which_angle = QUADRATIC_ANGLE;
      }

/**************************************************/
/** recognise torsion bending potentials **********/
/** and decide which one is present      **********/
/** Initially working with pcff frc file **********/
/**************************************************/
    found_torsion  = FALSE;
    which_torsion = FAILED_TORSION;
    if (locate_string_in_string(TORSION_1_STRING, &line_string[0], num_of_chars))
      {
        if (DEBUG) printf("GOT TORSION_1_STRING\n");
        found_torsion = TRUE;
        which_torsion = TORSION_1;
      }
    else if (locate_string_in_string(TORSION_3_STRING, &line_string[0], num_of_chars))
      {
        if (DEBUG) printf("GOT TORSION_3_STRING\n");
        found_torsion = TRUE;
        which_torsion = TORSION_3;
      }

/**************************************************/
/** recognise van der waals interaction potentials*/
/** and which one is present **********************/
/**************************************************/
    found_vdw=FALSE;
    this_type= FAILED_VDW;
    if (locate_string_in_string(NON_BOND_VDW_12_6, &line_string[0], num_of_chars))
      {
         found_vdw=TRUE;
         this_type=NB_12_6;
      }
    else if (locate_string_in_string(NON_BOND_VDW_9_6, &line_string[0], num_of_chars))
      {
         found_vdw=TRUE;
         this_type=NB_9_6;
      }
    else if (locate_string_in_string(NON_BOND_VDW_BUCK, &line_string[0], num_of_chars))
      {
         found_vdw=TRUE;
         this_type=NB_BUCK;
      }

    if (found_vdw)
      {
        if (this_type == FAILED_VDW)
          {
            printf("ERROR >> Unrecognised bond potential found in forcefield file\n");
            exit(0);
          }

/***************************************************/
/***** So we are at the grepped for title  *********/
/***** line !!!                            *********/
/***************************************************/
       
/***************************************************/
/***** Read info lines starting with @ symbols *****/
/***************************************************/

        while (!locate_string_in_string(INFO_LINE, &line_string[0], num_of_chars))
           {
              num_of_chars=read_line(file_fp, &line[0]);
              idummy= string_from_int(&line[0], &line_string[0]);
           }

        get_nonbond_info(file_fp, &line[0], 
                                 &line_string[0], &num_of_chars);  

/*****************************************************/
/**** Check the non-bond type is sensible ************/
/*****************************************************/

        if (strcmp(pot_info.type, R_EPS) == 0)
          {
            fprintf(output_fp,"type checks out: potentials in r, eps form\n");
          }
        else if (strcmp(pot_info.type, A_B) == 0)
          {
            fprintf(output_fp,"type checks out: potentials in A, B form\n");
          }
        else if (strcmp(pot_info.type, BUCK) == 0)
          {
            fprintf(output_fp,"type checks out: potentials in Buckingham form\n");
          }
        else
          {
             fprintf(output_fp,"ERROR: Do not understand the non-bond potential type >>%s<< found in the frc file\n",
                                                                                                      pot_info.type);
             fprintf(output_fp,"Currently available types are: r-eps and A-B\n");
             fflush(output_fp);
             fflush(stdout);
             exit(0);
          }

/*****************************************************/
/***** Now get the values !!! ************************/
/*****************************************************/

        *p_have_comb_rules=TRUE;
        if (strcmp(pot_info.combination, NONE) == 0) *p_have_comb_rules=FALSE;

        if (*p_have_comb_rules )
          {
            fprintf(output_fp,"Combinations rules read as: >>%s<<\n", pot_info.combination);
          }
        else
          {
            fprintf(output_fp,"No combining rules set\n");
          }

        get_nonbond_params(file_fp, &line[0], &line_string[0], 
                           &num_of_chars, *p_have_comb_rules,
                           this_type); 

/*****************************************************/
/**** Prepare values required by combining rules *****/
/**** Set bits that aren't required to zero      *****/
/*****************************************************/

       if (strcmp(pot_info.combination, SIXTH_POWER) == 0)
           {
             fprintf(output_fp,"Sixth Power combining rules in operation\n");
             for  (iloop=0; iloop <= num_potential_types; iloop++)
               {
                  a2 =  potent[iloop].a* potent[iloop].a;
                  potent[iloop].a3 = potent[iloop].a * a2;
                  potent[iloop].a6 = potent[iloop].a3 * potent[iloop].a3;
                  potent[iloop].sqrt_b = sqrt(potent[iloop].b);
                  potent[iloop].sqrt_a = 0;
               }  
           }
       else if (strcmp(pot_info.combination, ARITHMETIC) == 0)
           {
             for  (iloop=0; iloop <= num_potential_types; iloop++)
               {
                 potent[iloop].sqrt_a = sqrt(potent[iloop].a);
                 potent[iloop].sqrt_b = sqrt(potent[iloop].b);
                 potent[iloop].a3     = 0;
                 potent[iloop].a6     = 0;
               }
           }
       else if (strcmp(pot_info.combination, GEOMETRIC) == 0)
           {
             for  (iloop=0; iloop <= num_potential_types; iloop++)
               {
                 potent[iloop].sqrt_a = sqrt(potent[iloop].a);
                 potent[iloop].sqrt_b = sqrt(potent[iloop].b); 
                 potent[iloop].a3     = 0;
                 potent[iloop].a6     = 0;
               }
           }
       else if (strcmp(pot_info.combination, NONE) != 0)
           {
             fprintf(output_fp,
                     "ERROR: Do not understand the non-bond potential combining rules >>%s<< found in the frc file\n",
                                                                                                pot_info.combination);
             fprintf(output_fp,"Currently available types are: sixth-power and geometric\n");
             fflush(output_fp);
             fflush(stdout);
           }
        fprintf(output_fp,"Read %d van der Waals potential parameters\n", num_potential_types);
      }
/*****************************************************/
/******* Deal with bond stretch param. reading *******/
/*****************************************************/
    else if ( found_stretch )
      {
        if (which_stretch == FAILED_STRETCH)
          {
            printf("ERROR >> Unrecognised bond potential found in forcefield file\n");
            exit(0);
          }

/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

        while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars)) 
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_stretch_params( file_fp, &line[0], 
                            &line_string[0], &num_of_chars, which_stretch );

      }

    else if ( found_angle )
      {
         if (which_angle == FAILED_ANGLE)
           {
             printf("ERROR >> Unrecognised angle potential found in forcefield file\n");
             exit(0);
           }

/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

        while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars)) 
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_angle_params( file_fp, &line[0], 
                            &line_string[0], &num_of_chars, which_angle );

      }

    else if ( found_torsion )
      {
        if (which_torsion == FAILED_TORSION)
          {
            printf("ERROR >> Unrecognised torsion potential found in forcefield file\n");
            exit(0);
          }

/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

        while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars)) 
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_torsion_params( file_fp, &line[0], 
                            &line_string[0], &num_of_chars, which_torsion );

      }

/*****************************************************/
/******* Read in the equivalence table ***************/
/*****************************************************/

    else if ( locate_string_in_string(&equivalence[0], &line_string[0], num_of_chars) )
      {
/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

         while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars))
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_equivalences( file_fp, &line[0],
                                             &line_string[0], &num_of_chars );

      }
/*****************************************************/
/******* Read in the h_bond parameters (AMBER only) **/
/******* Requires sorting                           **/
/*****************************************************/

   else if ( locate_string_in_string(&h_bond[0], &line_string[0], num_of_chars) 
             && pot_info.version == AMBER )
      {
/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

        strcpy(pot_info.hbond, &h_bond[0]);

        while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars))
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_hbond_params( file_fp, &line[0],
                                             &line_string[0], &num_of_chars );
        
      }

/*****************************************************/
/******* Read in bond increment information     ******/
/******* Introduced for AMBER only March 99 DJW ******/
/******* Generalised for other potentials       ******/
/******* August 06                          DJW ******/
/*****************************************************/
   else if ( locate_string_in_string(&bond_inc[0], &line_string[0], num_of_chars))
      {
         if (DEBUG) printf("Found bond increment stuff in: >>%s<<\n", line_string);
         while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars))
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }
        get_bond_incs( file_fp, &line[0], &line_string[0], &num_of_chars );
      }


/*****************************************************/
/******* End of titles recognition if blocks *********/
/*****************************************************/

   }

/*****************************************************/
/** Report angle bonding potentials ******************/
/*****************************************************/

if (DEBUG) 
  {
    printf("Found %d angle potential types:\n", num_angles);
    for (iloop=0; iloop <= num_angles; iloop++)
      {
        if (intra_angle_potent[iloop].which == QUADRATIC_ANGLE)
          {
            printf("%d >>%s<< >>%s<< >>%s<< %10.6f %10.6f\n", iloop,
                                      intra_angle_potent[iloop].atom1,
                                      intra_angle_potent[iloop].atom2,
                                      intra_angle_potent[iloop].atom3,
                                      intra_angle_potent[iloop].A,
                                      intra_angle_potent[iloop].B);
          }
        else if (intra_angle_potent[iloop].which == QUARTIC_ANGLE)
          {
            printf("%d >>%s<< >>%s<< >>%s<< %10.6f %10.6f %10.6f %10.6f\n", iloop,
                                      intra_angle_potent[iloop].atom1,
                                      intra_angle_potent[iloop].atom2,
                                      intra_angle_potent[iloop].atom3,
                                      intra_angle_potent[iloop].A,
                                      intra_angle_potent[iloop].B,
                                      intra_angle_potent[iloop].C,
                                      intra_angle_potent[iloop].D);
          }
     }
  }


/*****************************************************/
/******* End of parameter reading while loop *********/
/******* Close the file behind you!!!        *********/
/*****************************************************/

fclose(file_fp);
return 0;
}


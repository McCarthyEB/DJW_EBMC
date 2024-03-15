/*********************************************************/
/***** Read the torsion potential parameters         *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  August 2006                     *****/
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "data.h"

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars );

int string_from_int(int *p_int, char *p_string);

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

void get_torsion_params(FILE *fp, int *p_line, char *p_line_string, 
                        int *p_num_of_chars, int which)
{
#include "header.h"
int ipoint, iend, itsanum;
int idummy, ref;
int num_digi, sign;

double version;

/**************************************************************/
/****** Ignore titles for parameters                    *******/
/**************************************************************/

   while (locate_string_in_string(TITLE_LINE,         p_line_string, *p_num_of_chars)
        || locate_string_in_string(ILLUSTRATION_LINE, p_line_string, *p_num_of_chars)
        || *p_num_of_chars == 0)
     {
        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }


/**************************************************************/
/****** Stay in the loop till we run out of parameter   *******/
/****** lines                                           *******/
/**************************************************************/

if (which == TORSION_1)
  {
    while(*p_num_of_chars > 0)
     {
        num_torsions++;
        intra_torsion_potent[num_torsions].which= TORSION_1;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

/***************************************************************/
/******* Read in potential type labels for this interaction ****/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom1, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom1[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom1[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom2, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom2[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom2[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom3, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom3[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom3[iend-ipoint-1]='\0';

        ipoint=iend;

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom4, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom4[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom4[iend-ipoint-1]='\0';

        ipoint=iend;

/***************************************************************/
/******* Read in data values for this interaction              */
/******* Assumed that the parameters are present in the order: */
/******* KPhi       n          Phi0              which become: */
/*******   A        B           C      D (0)   E (0)   F (0)   */ 
/***************************************************************/

        intra_torsion_potent[num_torsions].A
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].B
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].C
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].D = 0.0;
   
        intra_torsion_potent[num_torsions].E = 0.0;
   
        intra_torsion_potent[num_torsions].F = 0.0;
   
/***************************************************************/
/****** Read in next parameter line  ***************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
   }
else if (which == TORSION_3)
   {
    while(*p_num_of_chars > 0)
     {
        num_torsions++;
        intra_torsion_potent[num_torsions].which= TORSION_3;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

/***************************************************************/
/******* Read in potential type labels for this interaction ****/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom1, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom1[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom1[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom2, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom2[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom2[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom3, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom3[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom3[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_torsion_potent[num_torsions].atom4, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_torsion_potent[num_torsions].atom4[iend-ipoint-1]=='_')
                    intra_torsion_potent[num_torsions].atom4[iend-ipoint-1]='\0';

/***************************************************************/
/******* Read in data values for this interaction              */
/******* Assumed that the parameters are present in the order: */
/******* V(1)  Phi1(0) V(2) Phi2(0) V(3) Phi3(0) which become: */
/*******  A      B      C     D      E     F                   */
/***************************************************************/

        intra_torsion_potent[num_torsions].A 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].B 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].C 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].D 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].E 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_torsion_potent[num_torsions].F 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

/***************************************************************/
/****** Read in next parameter line  ***************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
   }
return;
}

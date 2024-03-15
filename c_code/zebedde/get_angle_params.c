/*********************************************************/
/***** Read the angle potential parameters           *****/
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

void get_angle_params(FILE *fp, int *p_line, char *p_line_string, 
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

if (which == QUARTIC_ANGLE)
  {
    while(*p_num_of_chars > 0)
     {
        num_angles++;
        intra_angle_potent[num_angles].which= QUARTIC_STRETCH;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

/***************************************************************/
/******* Read in potential type labels for this interaction ****/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom1, 
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom2, 
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom3, 
                                    p_line_string+ipoint, iend-ipoint);
        ipoint=iend;

/***************************************************************/
/******* Read in data values for this interaction              */
/******* Assumed that the parameters are present in the order: */
/******* THETA0     K2          K3          K4   which become: */
/*******   A        B           C            D                 */
/***************************************************************/

        intra_angle_potent[num_angles].A
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_angle_potent[num_angles].B
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_angle_potent[num_angles].C
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_angle_potent[num_angles].D
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
   
/***************************************************************/
/****** Read in next parameter line  ***************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
   }
else if (which == QUADRATIC_ANGLE)
   {
    while(*p_num_of_chars > 0)
     {
        num_angles++;
        intra_angle_potent[num_angles].which= QUADRATIC_ANGLE;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

/***************************************************************/
/******* Read in potential type labels for this interaction ****/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom1, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_angle_potent[num_angles].atom1[iend-ipoint-1]=='_')
                    intra_angle_potent[num_angles].atom1[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom2, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_angle_potent[num_angles].atom2[iend-ipoint-1]=='_')
                    intra_angle_potent[num_angles].atom2[iend-ipoint-1]='\0';

        ipoint=iend;
        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (intra_angle_potent[num_angles].atom3, 
                                    p_line_string+ipoint, iend-ipoint);

/*** strip trialing underscores ****/
        if (intra_angle_potent[num_angles].atom3[iend-ipoint-1]=='_')
                    intra_angle_potent[num_angles].atom3[iend-ipoint-1]='\0';

/***************************************************************/
/******* Read in data values for this interaction              */
/******* Assumed that the parameters are present in the order: */
/******* THETA0     K2                           which become: */
/*******    A        B       C (0)    D (0)                    */
/***************************************************************/

        intra_angle_potent[num_angles].A 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_angle_potent[num_angles].B 
                 = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

        intra_angle_potent[num_angles].C = 0;

        intra_angle_potent[num_angles].D = 0;
   
/***************************************************************/
/****** Read in next parameter line  ***************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
   }
return;
}

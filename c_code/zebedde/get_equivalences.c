/*********************************************************/
/***** Read the information about nonbond parameters *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  November 1995                   *****/
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

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

void get_equivalences(FILE *fp, int *p_line, char *p_line_string, 
                                int *p_num_of_chars)
{
#include "header.h"
int ipoint, iend, itsanum;
int idummy, ref;
int num_digi, sign;

double version;

/**************************************************************/
/****** Ignore titles for parameters                    *******/
/**************************************************************/


   while ( locate_string_in_string(TITLE_LINE,       p_line_string, *p_num_of_chars)
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

   num_equivalences=-1;
   while(*p_num_of_chars > 0)
     {

        num_equivalences++;
 
        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

/***************************************************************/
/******* Read in potential type that is to be checked **********/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].type, 
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/******* Read in Non-bond equivalence **************************/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].nonbond, 
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/******* Read in stretch equivalence ***************************/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].stretch,
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/******* Read in angle equivalence *****************************/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].angle,
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/******* Read in torsion equivalence ***************************/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].torsion,
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/******* Read in oop equivalence *******************************/
/***************************************************************/

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

        strncpy (equivalence_list[num_equivalences].oop,
                                    p_line_string+ipoint, iend-ipoint);

        ipoint=iend;

/***************************************************************/
/****** Read in next parameter line  ***************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
return;
}

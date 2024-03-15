/*********************************************************/
/***** Read the information about hbond parameters   *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  Febuaury 1999                   *****/
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

void get_hbond_params(FILE *fp, int *p_line, char *p_line_string, 
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
        || locate_string_in_string(JUST_RETURN      , p_line_string, *p_num_of_chars)
        || *p_num_of_chars == 0 || *p_num_of_chars == -1)
     {
        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }

/**************************************************************/
/****** Stay in the loop till we run out of parameter   *******/
/****** lines                                           *******/
/**************************************************************/

   num_hbonds=-1;
   while(*p_num_of_chars > 0)
     {
/**** Ignore commented parameter lines ********/
        if  (!locate_string_in_string(TITLE_LINE, p_line_string, *p_num_of_chars))
          {

             printf("Reading hbond parameters from: >>%s<<\n", p_line_string);

             num_hbonds++;
 
             ipoint= 0;

             version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
             ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

             ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
             iend= next_space(p_line, ipoint, *p_num_of_chars);
             strncpy (h_potent[num_hbonds].pot, 
                                    p_line_string+ipoint, iend-ipoint);

             ipoint= iend;
             ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
             iend= next_space(p_line, ipoint, *p_num_of_chars);
             strncpy (h_potent[num_hbonds].pot2, 
                                    p_line_string+ipoint, iend-ipoint);

             ipoint= iend;
             h_potent[num_hbonds].a
                      = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

             h_potent[num_hbonds].b= 
                        get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

             printf("Found: %s %s %10.6f %10.6f\n",h_potent[num_hbonds].pot,
                      h_potent[num_hbonds].pot2, h_potent[num_hbonds].a, 
                      h_potent[num_hbonds].b);
   
           }  
/***************************************************************/
/****** Read in next parameter line  **************************/
/***************************************************************/

         *p_num_of_chars=read_line(fp, p_line);
         idummy= string_from_int(p_line, p_line_string);
     }
return;
}

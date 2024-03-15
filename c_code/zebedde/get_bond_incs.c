/*********************************************************/
/***** Read the information about bond increment     *****/
/***** from lines begining with @ !                  *****/
/***** Dave Willock  March 1999                      *****/
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

void get_bond_incs(FILE *fp, int *p_line, char *p_line_string, 
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

   num_bond_incs=-1;
   while(*p_num_of_chars > 0)
     {
/**** Ignore commented parameter lines ********/
        if  (!locate_string_in_string(TITLE_LINE, p_line_string, *p_num_of_chars))
          {

             if (DEBUG) printf("Reading bond inc  parameters from: >>%s<<\n", p_line_string);

             num_bond_incs++;
 
             ipoint= 0;

             version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
             ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

             ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
             iend= next_space(p_line, ipoint, *p_num_of_chars);
             strncpy (bond_increments[num_bond_incs].pot, 
                                    p_line_string+ipoint, iend-ipoint);

             ipoint= iend;
             ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
             iend= next_space(p_line, ipoint, *p_num_of_chars);
             strncpy (bond_increments[num_bond_incs].pot2, 
                                    p_line_string+ipoint, iend-ipoint);

             ipoint= iend;
             bond_increments[num_bond_incs].a
                      = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

             bond_increments[num_bond_incs].b= 
                        get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

             if (DEBUG) printf("Found: Version %10.6f %s %s %10.6f %10.6f\n",version,
                        bond_increments[num_bond_incs].pot,
                        bond_increments[num_bond_incs].pot2, bond_increments[num_bond_incs].a, 
                        bond_increments[num_bond_incs].b);

             if (version < 1.0) exit(0);
   
           }  
/***************************************************************/
/****** Read in next parameter line  **************************/
/***************************************************************/

         *p_num_of_chars=read_line(fp, p_line);
         idummy= string_from_int(p_line, p_line_string);
     }
return;
}

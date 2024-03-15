/*********************************************************/
/***** Read the parameters for nonbond potentials    *****/
/***** Dave Willock  November 1995                   *****/
/***** Added Oie buckingham potentials Jan 06        *****/
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

void get_nonbond_params(FILE *fp, int *p_line, char *p_line_string, 
                        int *p_num_of_chars, int have_comb_rules, int type)
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

   while(*p_num_of_chars > 0)
     {

        if (DEBUG) printf("Reading vdw parameters from: >>%s<<\n", p_line_string);

        num_potential_types++;
 
        potent[num_potential_types].type = type;

        ipoint= 0;

        version = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);
        ref=  get_int(p_line, &ipoint, &itsanum, &num_digi,
                                               *p_num_of_chars, &sign);

        ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
        iend= next_space(p_line, ipoint, *p_num_of_chars);

/******************************************************/
/** If there are combining rules potentials are  ******/
/** made up of single element contributions so   ******/
/** will only read one label.                    ******/
/** If not there must be a unique potential for  ******/
/** every atom type pair and so need two labels  ******/
/** Each potential carries use_comb to refer to  ******/
/** later in the energy calculation. This will   ******/
/** allow mixing of potentials with and without  ******/
/** combination rules if someone wants it in the ******/
/** future! Dave Willock Jan 06.                 ******/
/******************************************************/

        if (have_comb_rules)
          {
            potent[num_potential_types].use_comb = TRUE;

            strncpy (potent[num_potential_types].pot, 
                                    p_line_string+ipoint, iend-ipoint);
            ipoint= iend;

            sprintf(potent[num_potential_types].pot2,"no");
   
            potent[num_potential_types].a
                     = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

            potent[num_potential_types].b= 
                   get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

            if (DEBUG) printf("Found: %s %10.6f %10.6f\n",potent[num_potential_types].pot,
                     potent[num_potential_types].a, potent[num_potential_types].b);
          }
        else
          {
            potent[num_potential_types].use_comb = FALSE;

            strncpy (potent[num_potential_types].pot, 
                                    p_line_string+ipoint, iend-ipoint);

            ipoint= iend;
            ipoint= next_none_space(p_line, ipoint, *p_num_of_chars);
            iend= next_space(p_line, ipoint, *p_num_of_chars);

            strncpy (potent[num_potential_types].pot2, 
                                    p_line_string+ipoint, iend-ipoint);

            potent[num_potential_types].a
                     = get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

            potent[num_potential_types].b= 
                       get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

/************************************************/
/** Read A,B,C for BUCK and just A and C for ****/
/** 12-6 or 9-6                              ****/
/************************************************/
            if (type == NB_BUCK)
             {
                potent[num_potential_types].c= 
                       get_doub(p_line, *p_num_of_chars, &ipoint, &itsanum);

                if (DEBUG) printf("Found: %s %s %10.6f %10.6f %10.6f\n",potent[num_potential_types].pot,
                                                         potent[num_potential_types].pot2,
                                                         potent[num_potential_types].a, 
                                                         potent[num_potential_types].b,
                                                         potent[num_potential_types].c);
             }
           else
             {
                if (DEBUG) printf("Found: %s %s %10.6f %10.6f\n",potent[num_potential_types].pot,
                                                         potent[num_potential_types].pot2,
                                                         potent[num_potential_types].a, 
                                                         potent[num_potential_types].b);
             }

          }
         
/***************************************************************/
/****** Read in next parameter line  **************************/
/***************************************************************/

        *p_num_of_chars=read_line(fp, p_line);
        idummy= string_from_int(p_line, p_line_string);
     }
return;
}

/***************************************************/
/*                                                 */
/*   update types list of template as new fragment */
/*                 is added                        */
/*                                                 */
/* expect num_from to be real number not upper     */
/* index. Dec. 2013                                */
/*                                                 */
/*                                                 */
/* Dave Willock August 1995                        */
/***************************************************/

#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int compare_strings( char *p_ichar1, char *p_ichar2 );

void copy_strings( char to[], char from[] );

void add_types( atom_number *p_types_to, int *p_num_to,
                atom_number *p_types_from, int num_types_from)
{
#include "header.h"
 int itypes_to, itypes_from;
 int index_to;

 int got_type;

 char *p_label_to, *p_label_from, *p_H;

 atom_number *p_this_to, *p_this_from;

 p_this_from= p_types_from-1;
 p_H= "H";
 p_label_to="none";

/**************************************************************************/
/************ Loop over the fragment types ********************************/
/**************************************************************************/
/*** Note that the num_types variables are true counters ******************/
/*** Setting num_to as equal to num_types_from puts it ready to index *****/
/*** the next new type, Dave Willock Feb 2014.                        *****/
/**************************************************************************/
/*** Note that if the number of types in the original guest (num_to)  *****/
/*** exceeds the number in the fragment we need to hold the number.   *****/
/*** Dave Willock Feb. 2014.                                          *****/
/**************************************************************************/

 index_to=num_types_from;
 for (itypes_from = 0; itypes_from < num_types_from; itypes_from++)
   {
     p_this_from++;
     p_this_to= p_types_to-1;
     p_label_from= &(p_this_from->atom_type[0]);
     got_type= FALSE;

//     printf("Adding type %s of which there are %d in fragment\n", p_label_from, p_this_from->num);

/**************************************************************************/
/************ Check to see if we already have any of these ****************/
/**************************************************************************/

     itypes_to = 0;
     while (itypes_to < *p_num_to && !got_type)
        {
           itypes_to++;
           p_this_to++;
           p_label_to= &(p_this_to->atom_type[0]);

           if ( compare_strings( p_label_to , p_label_from )) got_type= TRUE; 
        }

/**************************************************************************/
/******** Alter template list: if this type is already present ************/
/******** increase the number for it accordingly, if it is a   ************/
/******** new type then add it to the template types at the end************/
/******** adjust number of hydrogens assuming bond formed      ************/
/**************************************************************************/

     if (got_type)
        {
            p_this_to->num += p_this_from->num+1;

//            printf("Already got this as %s in guest, now total is %d\n", 
//                                                   p_this_to->atom_type,  p_this_to->num);

/*** Done later in code now, Feb 2014 Dave Willock   **************************/
/*            if (compare_strings( p_label_to , p_H )) (p_this_to->num) -= 2; */
/******************************************************************************/
        }
     else
        {
            *(p_types_to+*p_num_to) = *p_this_from;
        
//            printf("This is new adding as type %d label %s num %d \n",
//                                 *p_num_to, (p_types_to+*p_num_to)->atom_type,   
//                                 (p_types_to+*p_num_to)->num );      
            (*p_num_to)++;
        }
  }
return;
}


/***********************************************************************/
/* relabel.c :   Relabels atoms sequentially                           */
/*               Dave Willock 1995                                     */
/*               Last altered July 1998                                */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"


char int_to_char( int num);

void relabel(atom *p_molecule, int num_to_label, int *p_num_elems_used, int reset)
{

int iatom, ielem, single_char;
int *p_elem;
int ithou, ihuns, itens, iunits;

char element[3];

atom *p_atom;

/********************************************************************/
/**** Reset the number of elements used counter *********************/
/********************************************************************/

if (reset)
  {
    for (p_elem= p_num_elems_used; p_elem < p_num_elems_used+NUM_ELEMENTS; p_elem++) 
       {
          *p_elem=0;
       }
  }

for (iatom=0; iatom <= num_to_label; iatom++)
  {
     p_atom= p_molecule+iatom;

     element[0]= p_atom->elem[0];
     element[1]= p_atom->elem[1];
     element[2]= '\0';
     printf("DEBUG : For atom %d Making copy of element >>%s<< to >>%s<<\n", 
                       iatom, p_atom->elem, element);

/* check if there is only one character for this ones element type */

     single_char= TRUE;

     if (element[1] != '\0' ) 
        {
           single_char= FALSE;
        }


     for (ielem = 0; ielem < NUM_ELEMENTS; ielem++)
       {
          if ( strcmp( &element[0],  &(period_table[ielem].elem[0]))== 0 )
                                                                       break; 
       }

if (ielem == NUM_ELEMENTS) 
   {
      printf ("ERROR : An element which is not in the periodic table \n ");
      printf ("        was sent to relabel.c. The element label was : %s",
                                                                 element);
      printf ("        Zebedde is bailing out.\n\n");
      exit(1);
   }
  
    
     p_elem= p_num_elems_used+ielem;

     (*p_elem)++;

/*****************************************************************/
/**** Altered to cope with thousands of atoms July 98 DJW ********/
/*****************************************************************/

     ithou= (*p_elem)/1000;
     ihuns= (*p_elem)/100- 10*ithou;
     itens= (*p_elem)/10- 100*ithou-10*ihuns;
     iunits= (*p_elem)- 10*itens -100*ihuns- 1000*ithou;

     if (single_char)
       {
         p_atom->label[0] = toupper( element[0]);
         if (itens == 0 && ihuns == 0  && ithou == 0)
           {
              p_atom->label[1] = int_to_char(iunits);
              p_atom->label[2] = '\0';
           }
         else if (ihuns == 0 && ithou == 0)
           {
              p_atom->label[1] = int_to_char(itens);
              p_atom->label[2] = int_to_char(iunits);
              p_atom->label[3] = '\0'; 
           }
         else if (ithou == 0) 
           {
              p_atom->label[1] = int_to_char(ihuns);
              p_atom->label[2] = int_to_char(itens);
              p_atom->label[3] = int_to_char(iunits);
              p_atom->label[4] = '\0'; 
           }
         else if (ithou < 10) 
           {
              p_atom->label[1] = 'A'+10*(ithou-1)+ihuns;
              p_atom->label[2] = int_to_char(itens);
              p_atom->label[3] = int_to_char(iunits);
              p_atom->label[4] = '\0';
           }
         else 
           {
              printf("ERROR - Currently any molecule requiring relabelling cannot");
              printf("contain more than 2000 atoms of the same element, exiting\n");
           }
      }
    else
      {
         p_atom->label[0] = toupper( element[0]);
         p_atom->label[1] = toupper( element[1]);
         if (itens == 0 && ihuns == 0 && ithou == 0)
           {
              p_atom->label[2] = int_to_char(iunits);
              p_atom->label[3] = '\0';
           }
         else if (ihuns == 0 && ithou == 0)
           {
              p_atom->label[2] = int_to_char(itens);
              p_atom->label[3] = int_to_char(iunits);
              p_atom->label[4] = '\0';
           }
         else if (ithou == 0)   
           {
              p_atom->label[2] = int_to_char(ihuns);
              p_atom->label[3] = int_to_char(itens);
              p_atom->label[4] = int_to_char(iunits);
              p_atom->label[5] = '\0';
           }
         else if (ithou == 1)
           {
              p_atom->label[2] = 'A'-1+10*ithou+ihuns;
              p_atom->label[3] = int_to_char(itens);
              p_atom->label[4] = int_to_char(iunits);
              p_atom->label[5] = '\0';
           }
         else 
           {
              printf("ERROR - Currently any molecule requiring relabelling cannot");
              printf("contain more than 2000 atoms of the same element, exiting\n");
           }

      }
  }

return;
}










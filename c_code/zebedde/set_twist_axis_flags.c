/*******************************************************/
/***** routine to rotate section of molecule   *********/
/***** assumes first atom in the list is fixed *********/
/***** Dave Willock April 21st 1995            *********/
/***** Adapted for intra-molecular energy terms*********/
/***** Oct. 2006 Dave Willock                  *********/
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

/* required functions ------------------------------------*/

double real_random(int done);
 
void unit_vector(double *p_vector, double *p_size);

void join_atoms(atom *p_A, atom *p_B, double *p_A_to_B);
 
int find_chunk(atom *p_molecule, int num_template_atoms, int *p_flag_list,
                int atom1, int atom2);

void gather_molecule(atom *p_molecule, int num_atoms, int which_mol);

/*--------------------------------------------------------*/

void set_twist_axis_flags(atom *guest_ptrs[], int num_guests,
                          list_partition *p_guest_demarc,
                          int *p_flag_list,
                          links *p_link_atoms, 
                          double *p_origin, double *p_axis, int which_mol)
{
#include "header.h"

double size;
int atom1,atom2, flags_set;
int index;
atom *p_atom1, *p_atom2;

list_partition *p_demarc;

/****************************************************************/
/**** Start of Executable lines *********************************/
/****************************************************************/
 p_demarc=p_guest_demarc+which_mol;

if (DEBUG) printf("DEBUG>> Arrived in set_twist_axis_flags\n");

// printf("DEBUG>> which_test = %d\n", which_test);
// printf("DEBUG>> Gathering molecule consisting of %d atoms\n", p_demarc->num);

  if (symm_set) index=which_mol*(num_symm_ops+2);
                                         else index=which_mol;

// printf("DEBUG>> Calling gather_moelcule from set_twist_axis_flags\n");

  gather_molecule(guest_ptrs[index], p_demarc->end, -1);

/*** randomise which is atom1 and which atom2 ***/

if (real_random(1) < 0.5)
  {
   atom1= p_link_atoms->start;
   atom2= p_link_atoms->end;
  }
else
  {
   atom2= p_link_atoms->start;
   atom1= p_link_atoms->end;
  }

   p_atom1= guest_ptrs[index]+atom1;
   p_atom2= guest_ptrs[index]+atom2;

// printf("Atom1: %s %10.6f  %10.6f  %10.6f \n", p_atom1->label, p_atom1->x, p_atom1->y, p_atom1->z);
// printf("Atom2: %s %10.6f  %10.6f  %10.6f \n", p_atom2->label, p_atom2->x, p_atom2->y, p_atom2->z);

/*********** flag the chunk of molecule that atom2 belongs to ****/

   if ( !find_chunk(guest_ptrs[index], p_demarc->end, p_flag_list, atom1, atom2))
      {
         printf("ERROR: Thats a suprise these two are in a ring in set_twist_axis_flags\n");
         exit(0);
      }

// if (DEBUG) printf("twisting around bond %d (%s) - %d (%s)\n", atom1, p_atom1->label,
//                                                               atom2, p_atom2->label);

   flags_set= TRUE;

   *p_origin    = p_atom1->x;
   *(p_origin+1)= p_atom1->y;
   *(p_origin+2)= p_atom1->z;

// printf("Origin:  %10.6f  %10.6f  %10.6f \n", origin[0], origin[1], origin[2]);

   join_atoms( p_atom1, p_atom2, p_axis); 
   unit_vector(p_axis, &size);

//   printf("Set origin =  %10.6f  %10.6f  %10.6f \n", *p_origin, *(p_origin+1), *(p_origin+2));
//   printf("Set axis   =  %10.6f  %10.6f  %10.6f \n", *p_axis, *(p_axis+1), *(p_axis+2));
  
return;
}

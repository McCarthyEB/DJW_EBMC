/*************************************************************************/
/* assemble_vdw_list.c                                                   */
/* Assemble a list of all atom pairs in a molecule that are in range for */
/* van der Waals interactions.                                           */
/*                                                                       */
/* If the just_count variable is set the routine just returns the total  */
/* number of vdw interactions that would be listed.                      */
/*                                                                       */
/* nb_ctf_2 (from header.h) should be the square of the cut off required */
/*                                                                       */
/* if NO combining rules are being used the vdw_list array will hold     */ 
/* the index of the potential to use. If combining rules are in use      */ 
/* each atom will have its potential type in the nb_list part of the     */ 
/* atom structure and so this will be used as required and the vdw_list  */ 
/* will contain only the atom indices with the pot_index entry set to -1.*/ 
/*                                                                       */
/* Started Dave Willock 18th August 2006                                 */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void look_up_equiv_pots(atom *p_molecule, atom *p_atom, 
                        char *p_atom_pot, pot_types *p_neigh_pot_list,
                        int which_type, int need_neighbs);

void min_image( double *x, double *y, double *z);

int assemble_vdw_list(atom *p_molecule, int num_atoms, vdw_interact_list *p_vdw_list, int just_count)
{
#include "header.h"
int iatom, jatom, ipot, ineigh, neigh_index2;
int pot_index=-10;
int num_vdws_listed;
int neigh_index, n;
char atom2_pot[4], atom_pot[4];
pot_types neigh_pot_list[5];

double dx, dy, dz, r2;
	
BOOLEAN do_this;
atom *p_atom1, *p_atom2, *p_neigh;

 num_vdws_listed = -1;


 n=0;
 for (iatom=0; iatom<=num_atoms; iatom++) n+=iatom;

 if (DEBUG)
   {
      printf("Assembling intra-molecular vdw list for %d atoms with cut-off %10.6f\n",
                       num_atoms+1, sqrt(nb_ctf_2));
      printf("Exhaustive list would require %d elements\n",n);
   }

 if (just_count)
   {
/************************************************************************/
/***  if just counting still have to test all distances against cut off */
/************************************************************************/

     for (iatom=0; iatom<=num_atoms; iatom++)
       {
         p_atom1 = p_molecule+iatom;	
      
         for (jatom=iatom+1; jatom<=num_atoms; jatom++)
           {

/*******************************************************************/
/*** Check atom neighbours of atom1 and all second neighbours too **/
/*******************************************************************/

              do_this=TRUE;
              for (neigh_index=0; neigh_index <= p_atom1->num_neigh; neigh_index++)
                {
                  ineigh = p_atom1->neighb[neigh_index]; 
                  if (ineigh==jatom) do_this = FALSE;
                  
                  if (do_this)
                    {
                       p_neigh=p_molecule+ineigh;
                       for (neigh_index2=0; neigh_index2 <= p_neigh->num_neigh; neigh_index2++)
                         if (p_neigh->neighb[neigh_index2]==jatom) do_this = FALSE;
                    }
                }

              if (do_this)
                {
                   p_atom2 = p_molecule+jatom;	
                   
                   dx = p_atom1->x - p_atom2->x;
                   dy = p_atom1->y - p_atom2->y;
                   dz = p_atom1->z - p_atom2->z;

                   min_image( &dx, &dy, &dz);

                   r2 = dx*dx +dy*dy +dz*dz;

                   if (r2 <= nb_ctf_2) num_vdws_listed++;
                }
            }
       }
    return num_vdws_listed;
   }

/******************************************************************/
/******* Look for vdws in molecule and then assign potential   ****/	
/******* from potent list.                                     ****/	
/******************************************************************/

     for (iatom=0; iatom<=num_atoms; iatom++)
       {
         p_atom1 = p_molecule+iatom;	
         strcpy( atom_pot, p_atom1->pot);

/*** Get non-bond equivalence if we will be referencing potent here */
         if (strcmp(pot_info.combination, NONE) == 0 )
                 look_up_equiv_pots(p_molecule, p_atom1, &atom_pot[0], &neigh_pot_list[0],
                                                                            EQUIV_NONBOND, FALSE);

/*** Run over other atoms avoiding double accounting ***************/
      
         for (jatom=iatom+1; jatom<=num_atoms; jatom++)
           {

/*******************************************************************/
/*** Check atom neighbours of atom1 and all second neighbours too **/
/*** to eliminate jatom if it should be ignored.                  **/
/*******************************************************************/

              do_this=TRUE;
              for (neigh_index=0; neigh_index <= p_atom1->num_neigh; neigh_index++)
                {
                  ineigh = p_atom1->neighb[neigh_index]; 
                  if (ineigh==jatom) do_this = FALSE;
                  
                  if (do_this)
                    {
                       p_neigh=p_molecule+ineigh;
                       for (neigh_index2=0; neigh_index2 <= p_neigh->num_neigh; neigh_index2++)
                         if (p_neigh->neighb[neigh_index2]==jatom) do_this = FALSE;
                    }
                }

              if (do_this)
                {
                   p_atom2 = p_molecule+jatom;	

/************************************************/
/*** Test to see if there are combining rules ***/
/*** If not look up potential to use in potent **/
/*** array.                                    **/
/************************************************/
                   if (strcmp(pot_info.combination, NONE) == 0 )
                      {
/*** Get non-bond equivalence **************************************/
                         look_up_equiv_pots(p_molecule, p_atom2, &atom2_pot[0], &neigh_pot_list[0],
                                                                               EQUIV_NONBOND, FALSE);

/**************************************************************************/
/** Look up potential pair in pots list ***********************************/
/**************************************************************************/

                         ipot = -1;
                         pot_index=0;
                         while (ipot < 0 && pot_index <= num_potential_types)
                           {
                             if (strcmp(potent[pot_index].pot, p_atom1->pot)==0 &&
                                 strcmp(potent[pot_index].pot2, p_atom2->pot)==0)

                                ipot = 1;

                             else if (strcmp(potent[pot_index].pot2, p_atom1->pot)==0 &&
                                      strcmp(potent[pot_index].pot, p_atom2->pot)==0)

                                ipot = 1;

                             else
                                pot_index++;

                           }
                         if (ipot == -1)
                           {
                              printf("WARNING: No intra-molecular van der Waals for >>%s<< with >>%s<<\n",
                                                             p_atom1->pot, p_atom2->pot);
                              printf("         equivalents  >>%s<< and >>%s<<\n", atom_pot, atom2_pot);
                           }
                         else
                           {
                              printf("Intra-molecular van der Waals for >>%s<< with >>%s<<\n",
                                                             p_atom1->pot, p_atom2->pot);
                              printf("equivalents  >>%s<< and >>%s<<\n", atom_pot, atom2_pot);
                              printf("Set as       >>%s<< and >>%s<<\n", potent[pot_index].pot, 
                                                                           potent[pot_index].pot2);
                           }
                      }
                   
                   dx = p_atom1->x - p_atom2->x;
                   dy = p_atom1->y - p_atom2->y;
                   dz = p_atom1->z - p_atom2->z;

                   min_image( &dx, &dy, &dz);

                   r2 = dx*dx +dy*dy +dz*dz;

                   if (r2 <= nb_ctf_2)
                     {
                        num_vdws_listed++;

/*** Test to see if there are combining rules ***/
                        if (strcmp(pot_info.combination, NONE) == 0 )
                          {
                             (p_vdw_list+num_vdws_listed)->iatm1 = iatom;
                             (p_vdw_list+num_vdws_listed)->iatm2 = jatom;
                             (p_vdw_list+num_vdws_listed)->ivdwpot = pot_index;
                          }
                        else
                          {
                             (p_vdw_list+num_vdws_listed)->iatm1 = iatom;
                             (p_vdw_list+num_vdws_listed)->iatm2 = jatom;
                             (p_vdw_list+num_vdws_listed)->ivdwpot = -1;
                          }
                     }
                }
            } /* end of jatom loop */
       } /* end of iatom loop */

/****DEBUG DEBUG ****/
//    if (DEBUG)
//     {
//    printf("\n\nvdw list now: %d long\n", num_vdws_listed);
//    for (ipot=0; ipot <= num_vdws_listed; ipot++)
//      {
//        printf("%d %d %d\n", (p_vdw_list+ipot)->iatm1,  
//                             (p_vdw_list+ipot)->iatm2, 
//                             (p_vdw_list+ipot)->ivdwpot);
//      }
//    }
    return num_vdws_listed;
}


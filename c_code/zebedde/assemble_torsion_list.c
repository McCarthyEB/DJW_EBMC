/*************************************************************************/
/* assemble_torsion_list.c                                               */
/* Assemble a list of all torsions in a molecule and the potential to use*/
/* for the energy calculation from intra_torsion_potent.                 */
/*                                                                       */
/* Routine also makes up the links list for allowed torsional degrees    */
/* of freedom based on supplied rules.                                   */
/*                                                                       */
/* If the just_count variable is set the routine just returns the total  */
/* number of torsions in the molecule, i.e. the length of the possible   */
/* torsion potentials list.                                              */
/*                                                                       */
/* Started Dave Willock 18th August 2006                                 */
/* Altered Dave Willock  7th November 2006                               */
/* Last altered Dave Willock  March 2013                                 */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

int check_allowed_torsions(atom *p_molecule, int num_atoms, 
                           int iatom, int iatom2,
                           int ineigh1, int ineigh2, links *p_link_atoms,
                           int *p_num_links, int just_count);

int find_chunk(atom *p_molecule, int num_template_atoms, int *p_flag_list,
                int atom1, int atom2);

void look_up_equiv_pots(atom *p_molecule, atom *p_atom, char *p_atom_pot, 
                        pot_types *p_neigh_pot_list,
                        int which_type, int need_neighbs);

int assemble_torsion_list(atom *p_molecule, int num_atoms, 
                          torsion_interact_list *p_torsions_list, 
                          links *p_links_list, int *p_num_links, int just_count)
{
#include "header.h"
int iatom, iatom2, jpot, use_pot=-1, ineigh, neigh_index1, neigh_index2;
int ineigh1, ineigh2, num_torsions_listed, iallowed;
int neigh_index, ni, nj, i, in, in2, in3;

/**** flag_list should be malloced as needed ****/
int flag_list[MAXTEMPLATE];

char atom_pot[4];
char atom2_pot[4];
pot_types neigh_pot_list[5];
pot_types neigh2_pot_list[5];
	
BOOLEAN found_centreij, found_centreji;
BOOLEAN found_11, found_22;
BOOLEAN found_thisij, found_thisji, found_pot_flag;
BOOLEAN bother;
atom *p_atom, *p_atom2, *p_neigh, *p_neigh1, *p_neigh2;

/*** DEBUG ***/
/*  DEBUG=TRUE;  */
/*** DEBUG set false at end again! ***/

 num_torsions_listed = -1;
// if (just_count)
//   {
// printf("DEBUG>> In assemble_torsions_list, just counting\n");
// printf("DEBUG>> Have num_atoms = %d\n", num_atoms);
/************************************************************************/
/***  if just counting use torsions at a central pair i,j is given by ***/
/***  (ni - 1)(nj -1)                                                 ***/
/***  with ni and nj the number of neighbours at the centres.         ***/
/***  Remembering that num_neighs in the atom structure is zero       ***/
/***  referenced.                                                     ***/
/***  Only worth considering the pair if they each have more than     ***/
/***  one neighbour.                                                  ***/
/************************************************************************/

//     for (iatom=0; iatom<=num_atoms; iatom++)
//       {
//          p_atom = p_molecule+iatom;	
//
//          printf("DEBUG>> Looking at atom %s which has %d neighbours\n",
//                                          p_atom->label, p_atom->num_neigh);
//          if (p_atom->num_neigh > 0)
//            {
//              ni = p_atom->num_neigh;

/**** loop over neighbours of iatom looking for suitable pairs in a torsion ****/
//              for (in=0; in<=ni; in++)
//                 {
/**** Avoid double counting***/
//                    iatom2= p_atom->neighb[in];
//                    if (iatom2 < iatom)
//                       {
//                         p_atom2 = p_molecule + iatom2;
//
//                         if (p_atom2->num_neigh > 0)
//                           {
//                             nj = p_atom2->num_neigh;
//
//                           num_torsions_listed += ni * nj;

/*** loop over neighbours of atom2 *****/

//                           for (in2=0; in2<=nj; in2++)
//                              {
//                                ineigh2=p_atom2->neighb[in2];
//                        
//                                if (ineigh2 != iatom)
//                                  {
//
//                                    for (in3=0; in3<=ni; in3++)
//                                       {
//                                         ineigh=p_atom->neighb[in3];
//                                         if (ineigh != iatom2)
//                                           {
//                           check_allowed_torsions(p_molecule, num_atoms, 
//                                                  iatom, iatom2, 
//                                                  ineigh, ineigh2,
//                                                  p_links_list, p_num_links, 
//                                                  just_count);
//                                           }              
//                                      }
//                                  }
//                             }
//                         }
//                     }
//               }
//          }
//     }
//    printf("DEBUG>> returning %d torsions\n", num_torsions_listed);
//  return num_torsions_listed;
// }
//
//printf("DEBUG>> In assemble_torsions_list, ");
 
//if (just_count) printf("DEBUG>> just counting...\n");
//               else printf("DEBUG>> ... now filling\n");

/********************************************************************
/******* Look for torsions in molecule and then assign potential ****/	
/******* from list.                                              ****/	
/******* To do this will consider each atom in the list (iatom)  ****/	
/******* and treat it as the central atom of the torsion, so any ****/	
/******* two neighbours will give a valid torsion to set.        ****/	
/********************************************************************/

  for (iatom=0; iatom<=num_atoms; iatom++)
    {
      p_atom = p_molecule+iatom;	
/******************************************************************/
/**** Only bother at all if there are two or more neighbours   ****/
/**** on both the atom and the particular neighbour.           ****/
/**** atom and atom2 refer to the central atoms of the torsion ****/
/**** Remember num_neigh is zero reference so 1 means 2 neighs!****/
/******************************************************************/
      if (p_atom->num_neigh > 0)
        {
/******************************************************************/
/******* Check equivalence table for this atom and its neighbs ****/
/******************************************************************/

           look_up_equiv_pots(p_molecule, p_atom, &atom_pot[0], &neigh_pot_list[0], 
                              EQUIV_TORSION, TRUE ); 
     
/******************************************************************/
/*** Try out each neighbour as a possible partner at the centre ***/
/******************************************************************/

           for (neigh_index=0; neigh_index<=p_atom->num_neigh; neigh_index++)
             {
/**** Avoid double entries ***/
               if (p_atom->neighb[neigh_index] < iatom)
                 {
                    iatom2  = p_atom->neighb[neigh_index];
                    p_atom2 = p_molecule + iatom2;
                    if (p_atom2->num_neigh > 0)
                      {
/******************************************************************/
/******* Check equivalence table for this atom and its neighbs ****/
/******************************************************************/

                         look_up_equiv_pots(p_molecule, p_atom2, &atom2_pot[0], &neigh2_pot_list[0],
                                            EQUIV_TORSION, TRUE ); 

/******************************************************************/
/*** Loop over none repeating pairs of neighbours *****************/
/*** neigh1 is a neighbour from atom1             *****************/
/*** neigh2 is a neighbour from atom2             *****************/
/*** So we are looking for torsions of type       *****************/
/*** neigh1---atom1---atom2---neigh2              *****************/
/******************************************************************/
                         for (neigh_index1=0; neigh_index1 <= p_atom->num_neigh; 
                                                                    neigh_index1++)
                           {
                             ineigh1 = p_atom->neighb[neigh_index1];

                             if (ineigh1 != iatom2)
                               {
                                  p_neigh1 = p_molecule+ineigh1;
   
                                  for (neigh_index2=0; neigh_index2 <= p_atom2->num_neigh; neigh_index2++)
                                    {
                                       ineigh2 = p_atom2->neighb[neigh_index2];
                                       if (ineigh2 != iatom)
                                         {
                                            p_neigh2 = p_molecule+ineigh2;

/*******************************************************************/
/*** At this point can also check if the atom quartet has been  ****/
/*** indicated as flexible by the user, if so add to links list ****/
/*** Dave Willock October 2006                                  ****/
/*** This process has been moved to the subroutine below        ****/
/*** Dave Willock March 2013.                                   ****/
/*******************************************************************/

                                            check_allowed_torsions(p_molecule, num_atoms, 
                                                                   iatom, iatom2, ineigh1, ineigh2,
                                                                   p_links_list, p_num_links, 
                                                                   just_count);

//                                            if (DEBUG)
//                                              {
//                                                 printf("DEBUG>> Looking for a torsion potential for >>%s<< >>%s<< >>%s<< >>%s<<\n",
//                                                               p_neigh1->pot, p_atom->pot, p_atom2->pot, p_neigh2->pot);
//               
//                                                 printf("DEBUG>> assigned equivalences               >>%s<< >>%s<< >>%s<< >>%s<<\n",
//                                                                                   neigh_pot_list[neigh_index1].pot, atom_pot, atom2_pot, 
//                                                                                   neigh2_pot_list[neigh_index2].pot);
//                                              }
                                      
                                            found_pot_flag = FALSE;
/******************************************************************/
/***** test if potential matches database potential           *****/
/***** We may have i--j centre while database has j--i        *****/
/*****          found_centreij                 found_centreji *****/
/******************************************************************/
                                            for (jpot=0; jpot<=num_torsions; jpot++)
                                               {
                                                 found_centreij  =    strcmp(atom_pot, intra_torsion_potent[jpot].atom2) == 0
                                                                   && strcmp(atom2_pot, intra_torsion_potent[jpot].atom3) == 0;

                                                 found_centreji  =    strcmp(atom_pot, intra_torsion_potent[jpot].atom3) == 0
                                                                   && strcmp(atom2_pot, intra_torsion_potent[jpot].atom2) == 0;

                                                 found_thisij = FALSE;
                                                 found_thisji = FALSE;
 
                                                 if (found_centreij)
                                                   {
                                                      found_11 = strcmp(neigh_pot_list[neigh_index1].pot, intra_torsion_potent[jpot].atom1) == 0 
                                                               ||strcmp("*", intra_torsion_potent[jpot].atom1) == 0;

                                                      found_22 = strcmp(neigh2_pot_list[neigh_index2].pot, intra_torsion_potent[jpot].atom4) == 0
                                                               ||strcmp("*", intra_torsion_potent[jpot].atom4) == 0;

                                                      found_thisij = found_11 && found_22;
                                                   }
                                                 if (found_centreji)
                                                   {
                                                      found_11 = strcmp(neigh2_pot_list[neigh_index2].pot, intra_torsion_potent[jpot].atom1) == 0 
                                                               ||strcmp("*", intra_torsion_potent[jpot].atom1) == 0;

                                                      found_22 = strcmp(neigh_pot_list[neigh_index1].pot, intra_torsion_potent[jpot].atom4) == 0
                                                               ||strcmp("*", intra_torsion_potent[jpot].atom4) == 0;

                                                      found_thisji = found_11 && found_22;
                                                   }

                                                  if (found_thisij || found_thisji)
                                                    {
//                                                        if (DEBUG)
//                                                          {
//                                                             if (intra_torsion_potent[jpot].which==TORSION_1)
//                                                                     printf("Identified torsion_1 potential %s %s %s %s as a possible\n",
//                                                                                    intra_torsion_potent[jpot].atom1,
//                                                                                    intra_torsion_potent[jpot].atom2,
//                                                                                    intra_torsion_potent[jpot].atom3,
//                                                                                    intra_torsion_potent[jpot].atom4);
//
//                                                            if (intra_torsion_potent[jpot].which==TORSION_3)
//                                                                    printf("Identified torsion_3 torsion potential %s %s %s %s as a possible\n",
//                                                                                    intra_torsion_potent[jpot].atom1,
//                                                                                    intra_torsion_potent[jpot].atom2,
 //                                                                                   intra_torsion_potent[jpot].atom3,
 //                                                                                   intra_torsion_potent[jpot].atom4);
 //                                                          }
              
                                                       use_pot = jpot;
                                                       found_pot_flag = TRUE;
                                                    }
                                               } /* end jpot loop */

                                            if (found_pot_flag)
                                               {
                                                   if (use_pot == -1)
                                                     {
                                                        printf("ERROR >> assemble_torsions routine did not set use_pot when required.\n");
                                                        printf("         Please inform authors.\n");
                                                        exit(0);
                                                     }
                                                   if (DEBUG) printf("Using potential: >>%s<< >>%s<< >>%s<< >>%s<<\n",
                                                                                    intra_torsion_potent[use_pot].atom1,
                                                                                    intra_torsion_potent[use_pot].atom2,
                                                                                    intra_torsion_potent[use_pot].atom3,
                                                                                    intra_torsion_potent[use_pot].atom4);
/*********************************************************************************/
/*** add the lowest lying potential to the torsion pots list *********************/
/*********************************************************************************/
                                                   num_torsions_listed++;

                                                   if (!just_count)
                                                     {
                                                       (p_torsions_list+num_torsions_listed)->iatm1= ineigh1;
                                                       (p_torsions_list+num_torsions_listed)->iatm2= iatom;
                                                       (p_torsions_list+num_torsions_listed)->iatm3= iatom2;
                                                       (p_torsions_list+num_torsions_listed)->iatm4= ineigh2;
                                                       (p_torsions_list+num_torsions_listed)->itorpot= use_pot;
                                                     }
                                               }
                                             if (!found_pot_flag)
                                               {

/******* no  potential present for this element *******************/

                                                 if (DEBUG)
                                                   {
                                                   printf("WARNING: Cannot assign intra torsion potential for Atoms:");
                                                   printf("         %s (pot=%s), %s (pot=%s), %s (pot=%s), %s (pot=%s)\n", 
                                                                                            p_neigh1->label, p_neigh1->pot, 
                                                                                            p_atom->label,   p_atom->pot, 
                                                                                            p_atom2->label,  p_atom2->pot, 
                                                                                            p_neigh2->label, p_neigh2->pot); 

                                                   printf("atom %s neighs are :",p_atom->label);
                                                   for (i=0; i<=p_atom->num_neigh; i++)
                                                     {
                                                        p_neigh= p_molecule+p_atom->neighb[i];
                                                        printf(" %s", p_neigh->label);
                                                     }
                                                   printf("\n\n");

                                                   printf("atom %s neighs are :",p_atom2->label);
                                                   for (i=0; i<=p_atom2->num_neigh; i++)
                                                     {
                                                        p_neigh= p_molecule+p_atom2->neighb[i];
                                                        printf(" %s", p_neigh->label);
                                                     }
                                                   printf("\n\n");
                                                   
                                                   /*            exit(EXIT_FAILURE); */
                                                   }
                                               }
                                             if (DEBUG) printf("-----\n");
                                          }
                                      } /* end neigh_index2 loop */
                                  }
                              } /* end neigh_index1 loop */
                          }
                     }
                 } /* end neigh_index loop */
             }
         } /* end iatom loop */

/*** DEBUG ***/
/* DEBUG=FALSE;*/ 
/*** DEBUG set false at end again! ***/

return num_torsions_listed;
}


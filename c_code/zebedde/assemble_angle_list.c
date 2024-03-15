/*************************************************************************/
/* assemble_angle_list.c                                                 */
/* Assemble a list of all angles in a molecule and the potential to use  */
/* for the energy calculation from intra_angle_potent.                   */
/*                                                                       */
/* If the just_count variable is set the routine just returns the total  */
/* number of angles in the molecule, i.e. the length of the possible     */
/* angle potentials list.                                                */
/*                                                                       */
/* Started Dave Willock 18th August 2006                                 */
/*                                                                       */
/* Updated for multi-molecule version, March 2013, Dave Willock          */
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

int assemble_angle_list(atom *p_molecule, int num_atoms, angle_interact_list *p_angles_list,
                        int just_count)
{
#include "header.h"
int iatom, jpot, neigh_index1, neigh_index2;
int ineigh1, ineigh2, num_angles_listed;
int n, i, isum, iwarn;
char atom_pot[4];
pot_types neigh_pot_list[5];
	
BOOLEAN found_centre, have_warned;
BOOLEAN found_11, found_22, found_21, found_12;
BOOLEAN found_this, found_pot_flag;
atom *p_atom, *p_neigh1, *p_neigh2;

 num_angles_listed = -1;

/******************************************************************/
/******* Look for angles in molecule and then assign potential ****/	
/******* from list.                                            ****/	
/******* To do this will consider each atom in the list (iatom)****/	
/******* and treat it as the central atom of the angle, so any ****/	
/******* two neighbours will give a valid angle to set.        ****/	
/******************************************************************/

  for (iatom=0; iatom<=num_atoms; iatom++)
    {
      p_atom = p_molecule+iatom;	
/******************************************************************/
/**** Only bother at all if there are two or more neighbours   ****/
/**** Remember num_neigh is zero reference so 1 means 2 neighs!****/
/******************************************************************/
      if (p_atom->num_neigh > 0)
        {
           strcpy( atom_pot, p_atom->pot);

/******************************************************************/
/******* Check equivalence table **********************************/
/******************************************************************/
           look_up_equiv_pots(p_molecule, p_atom, &atom_pot[0], &neigh_pot_list[0],
                              EQUIV_ANGLE, TRUE);

/******************************************************************/
/*** Loop over none repeating pairs of neighbours *****************/
/******************************************************************/
      for (neigh_index1=0; neigh_index1 <= p_atom->num_neigh; neigh_index1++)
        {
          ineigh1 = p_atom->neighb[neigh_index1];
          p_neigh1 = p_molecule+ineigh1;
   
          for (neigh_index2=neigh_index1+1; neigh_index2 <= p_atom->num_neigh; neigh_index2++)
            {
              ineigh2 = p_atom->neighb[neigh_index2];
              p_neigh2 = p_molecule+ineigh2;
               
              if (DEBUG)
                {
                   printf("DEBUG HERE>> Looking for an angle potential for %s %s %s >>%s<< >>%s<< >>%s<<\n",
                                            p_neigh1->label, p_atom->label, p_neigh2->label,
                                            p_neigh1->pot, p_atom->pot, p_neigh2->pot);
               
                   printf("DEBUG>> assigned equivalences              >>%s<< >>%s<< >>%s<<\n",
                                                                neigh_pot_list[neigh_index1].pot, atom_pot, 
                                                                neigh_pot_list[neigh_index2].pot);
                }
                                      
              found_pot_flag = FALSE;

/******************************************************************/
/** First time look for exact matches only   **********************/
/** The outer 2 atoms may occur in any order **********************/
/******************************************************************/
              for (jpot=0; jpot<=num_angles; jpot++)
                {
                   found_centre = strcmp(atom_pot, intra_angle_potent[jpot].atom2) == 0;

                   if (found_centre)
                     {
                   
                       found_11 = strcmp(neigh_pot_list[neigh_index1].pot, intra_angle_potent[jpot].atom1) == 0;

                       found_12 = strcmp(neigh_pot_list[neigh_index1].pot, intra_angle_potent[jpot].atom3) == 0;

                       found_21 = strcmp(neigh_pot_list[neigh_index2].pot, intra_angle_potent[jpot].atom1) == 0;

                       found_22 = strcmp(neigh_pot_list[neigh_index2].pot, intra_angle_potent[jpot].atom3) == 0;

                       found_this = (found_11 && found_22) || (found_21 && found_12);

                       if (found_this)
                          {
                             if (DEBUG)
                               {
                                   if (intra_angle_potent[jpot].which==QUADRATIC_ANGLE)
                                     printf("Identified quadratic angle potential %s %s %s as a possible\n",
                                                             intra_angle_potent[jpot].atom1,
                                                             intra_angle_potent[jpot].atom2,
                                                             intra_angle_potent[jpot].atom3);

                                   if (intra_angle_potent[jpot].which==QUARTIC_ANGLE)
                                     printf("Identified quartic angle potential %s %s %s as a possible\n",
                                                             intra_angle_potent[jpot].atom1,
                                                             intra_angle_potent[jpot].atom2,
                                                             intra_angle_potent[jpot].atom3);
              
                                   printf("Because:\n");
                                   if (found_11) printf("found_11, >>%s<< >>%s<<\n", neigh_pot_list[neigh_index1].pot, 
                                                                                      intra_angle_potent[jpot].atom1);
                                   if (found_12) printf("found_12, >>%s<< >>%s<<\n",neigh_pot_list[neigh_index1].pot,
                                                                                      intra_angle_potent[jpot].atom3);
                                   if (found_21) printf("found_21, >>%s<< >>%s<<\n",neigh_pot_list[neigh_index2].pot,
                                                                                      intra_angle_potent[jpot].atom1);
                                   if (found_22) printf("found_22, >>%s<< >>%s<<\n",neigh_pot_list[neigh_index2].pot,
                                                                                intra_angle_potent[jpot].atom3);
                               }

/*********************************************************************************/
/*** add this atom set and potential to the list *********************************/
/*** Added "if" August 07 to prevent multiple potentials  ************************/
/***                                         for the same angle. *****************/
/*** This will always use the potential that is lowest in the frc list ***********/
/***                                                             Dave Willock ****/
/*********************************************************************************/
                             if (!found_pot_flag) num_angles_listed++;
                             found_pot_flag = TRUE;
 
                             if (!just_count)
                               {

                             (p_angles_list+num_angles_listed)->iatm1= ineigh1;
                             (p_angles_list+num_angles_listed)->iatm2= iatom;
                             (p_angles_list+num_angles_listed)->iatm3= ineigh2;
                             (p_angles_list+num_angles_listed)->iangpot= jpot;

                             if (DEBUG) printf("DEBUG >> list indices: %d %d %d with pot %d\n",
                                                        (p_angles_list+num_angles_listed)->iatm1,
                                                        (p_angles_list+num_angles_listed)->iatm2,
                                                        (p_angles_list+num_angles_listed)->iatm3,
                                                        (p_angles_list+num_angles_listed)->iangpot);
                                         
                              }
                          }
                     }

/******************************************************************/
/***** test if potential matches database potential ***************/
/******************************************************************/

                } /* end of jpot loop */
             if (DEBUG) printf("-----\n");

             if (!found_pot_flag)
               {
/*** See if we have already mentioned this one ***/

                 have_warned = FALSE;
                 for (iwarn=0; iwarn <= num_angle_warnings; iwarn++)
                   {
                      if (    (  strcmp(p_neigh1->pot, angle_warned[iwarn].atom1) == 0
                              && strcmp(p_neigh2->pot, angle_warned[iwarn].atom3) == 0 )
                          ||  (  strcmp(p_neigh1->pot, angle_warned[iwarn].atom3) == 0
                              && strcmp(p_neigh2->pot, angle_warned[iwarn].atom1) == 0 ) )
                          {
                               if ( strcmp(p_atom->pot, angle_warned[iwarn].atom2) == 0 )
                                  {
                                    have_warned = TRUE;
                                  }
                          }
                   }

/******* no  potential present for this element *******************/
                   if (!have_warned)
                     {
                 printf("Current warnings list:\n");
                 for (iwarn=0; iwarn <= num_angle_warnings; iwarn++)
                   {
                      printf(">>%s<< >>%s<< >>%s<<\n", angle_warned[iwarn].atom1, angle_warned[iwarn].atom2, angle_warned[iwarn].atom3);
                   }

                       num_angle_warnings++;
                       printf("WARNING: cannot assign an angle potential for %s %s %s >>%s<< >>%s<< >>%s<< %d unique angle warnings now issued.\n",
                                                p_neigh1->label, p_atom->label, p_neigh2->label,
                                                p_neigh1->pot, p_atom->pot, p_neigh2->pot, num_angle_warnings+1);

                       if (num_angle_warnings == MAX_ANGLE_WARN_LIST)
                         {
                            printf("ERROR: With this many potential warnings ZEBEDDE has to exit.....\n");
                            printf("       Check the atom potentials frc file matches your molecule files.\n");
                            exit(0);
                         }

                       strcpy(angle_warned[num_angle_warnings].atom1, p_neigh1->pot);
                       strcpy(angle_warned[num_angle_warnings].atom2, p_atom->pot);
                       strcpy(angle_warned[num_angle_warnings].atom3, p_neigh2->pot);
                     }
               
/*            exit(EXIT_FAILURE); */
              }
          }/* end loop neigh_index2 */
       }/* end loop neigh_index1 */
     }/* close of num_neigh if test */
  }/* end loop iatom */


return num_angles_listed;
DEBUG=0;
}


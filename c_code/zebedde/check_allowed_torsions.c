/*********************************************************/
/*** Check that a pair of atoms sent to the routine are **/
/*** suitable to be in the torsional degrees of freedom **/
/*** list.                                              **/
/*** Atoms are sent in the order:                       **/
/*** ineigh1----iatom---iatom2----ineigh2               **/
/***                                                    **/
/*** Added by Dave Willock, March 2013                  **/
/*********************************************************/
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int find_chunk(atom *p_molecule, int num_template_atoms, int *p_flag_list,
                int atom1, int atom2);

int check_allowed_torsions(atom *p_molecule, int num_atoms, int iatom, int iatom2,
                           int ineigh1, int ineigh2, links *p_link_atoms,
                           int *p_num_links, int just_count)
{
#include "header.h"

BOOLEAN bother, have_links;
BOOLEAN found_centreij, found_centreji; 
int *p_flags;
int i, iallowed;
char link_end1[4];
char link_end2[4];

links *p_this_link;

atom *p_atom, *p_atom2, *p_neigh1, *p_neigh2;

 have_links=FALSE;
 p_atom=p_molecule+iatom;
 p_atom2=p_molecule+iatom2;
 p_neigh1=p_molecule+ineigh1;
 p_neigh2=p_molecule+ineigh2;

   if (DEBUG)
    {

      printf("DEBUG>> Looking for %d >>%s<< %d >>%s<< %d >>%s<< %d >>%s<<",
                 ineigh1, p_neigh1->elem, iatom, p_atom->elem, iatom2,
                 p_atom2->elem, ineigh2, p_neigh2->elem);
      printf(" in allowed torsions list\n");
    }

/*******************************************************/ 
/*** gaurd against duplication *************************/ 
/*******************************************************/ 
   bother = TRUE;

   if (!just_count)
     {
       p_this_link= p_link_atoms;
       for (i=0; i<=*p_num_links; i++)
         {
           if (   p_this_link->start == iatom 
               && p_this_link->end   == iatom2 )
             {
                bother=FALSE;
                if (DEBUG) printf("No need link %d has start %d end %d\n", i, p_this_link->start,
                                                               p_this_link->end);
                break;
             }
           else if (   p_this_link->start == iatom2
                    && p_this_link->end   == iatom  )
             {
                bother=FALSE;
                if (DEBUG) printf("No need link %d has start %d end %d\n", i, p_this_link->start,
                                                               p_this_link->end);
                break;
             }
           p_this_link++;
         } 
     }

if (bother)
  {
/*******************************************************/
/*** malloc flags for this particular problem **********/
/*******************************************************/

    p_flags=(int*)malloc( (num_atoms+1) *sizeof(int));

    if (p_flags == NULL)
     { 
       printf("ERROR: Problem with memory assignment for %d flags in check_allowed_torsions.\n",
                                                num_atoms+1);
       exit(0);
     }

/****** Only bother with links that involve bonds not linked by a ring ******/
                
     if (!find_chunk(p_molecule, num_atoms, p_flags, iatom, iatom2))
       {
         bother=FALSE; 
         if (DEBUG) printf("No need these are in same ring\n");
       }

     free(p_flags);
    }

/***** So now see if this pair are allowed to be a moveable torsion ***/
 
 if (bother)
   {
/****** treat HA and HB as normal H atoms ***********************************/
/****** treat D as H for links, update added Jan. 09 Dave Willock ***********/

     if (   strcmp(p_neigh1->elem, "HA") == 0 || strcmp(p_neigh1->elem, "HB") == 0
                                              || strcmp(p_neigh1->elem, "D") == 0 )   
       {
          strcpy(link_end1, "H");
       }
     else
       {
          strcpy(link_end1, p_neigh1->elem);
       }

     if (   strcmp(p_neigh2->elem, "HA") == 0 || strcmp(p_neigh2->elem, "HB") == 0
                                              || strcmp(p_neigh2->elem, "D") == 0 )   
       {
          strcpy(link_end2, "H");
       }
     else
       {
          strcpy(link_end2, p_neigh2->elem);
       }

      if (DEBUG) 
       {
         printf("For %s %s as ends", p_neigh1->elem, p_neigh2->elem); 
         printf(" using %s %s\n", link_end1, link_end2); 
         printf("Testing for links against list of %d possibles\n", num_allowed_torsions);
         printf("\nFor torsion %s %s %s %s\n", link_end1, p_atom->elem, p_atom2->elem,
                                               link_end2);
       }

     for (iallowed = 0; iallowed<=num_allowed_torsions; iallowed++)
       {
          found_centreij  =    strcmp(p_atom->elem, allowed_torsions[iallowed].B) == 0
                            && strcmp(p_atom2->elem, allowed_torsions[iallowed].C) == 0;

          found_centreji  =    strcmp(p_atom2->elem, allowed_torsions[iallowed].B) == 0
                            && strcmp(p_atom->elem, allowed_torsions[iallowed].C) == 0;

          if (found_centreij)
            {
              if (  strcmp(link_end1, allowed_torsions[iallowed].A) == 0 
                  &&strcmp(link_end2, allowed_torsions[iallowed].D) == 0 )
                 {
                    if (DEBUG)
                       printf("Allowed case %d matches with: >>%s<<  >>%s<<  >>%s<<  >>%s<<\n",
                            iallowed, allowed_torsions[iallowed].A, 
                                      allowed_torsions[iallowed].B,
                                      allowed_torsions[iallowed].C,
                                      allowed_torsions[iallowed].D);
                                                               
                    ++*p_num_links;  

                    if (!just_count)
                      {
                         p_this_link= p_link_atoms+*p_num_links;
                         p_this_link->start=iatom;
                         p_this_link->end  =iatom2;
     
                         if (DEBUG) printf("DEBUG>> Just got ij link %d for atoms %d %d\n", 1+*p_num_links, p_this_link->start,
                                                                                 p_this_link->end);
                      }
                    break;
                 }
            }
          if (found_centreji)
            {
              if (  strcmp(link_end1, allowed_torsions[iallowed].D) == 0 
                  &&strcmp(link_end2, allowed_torsions[iallowed].A) == 0 )
                 {
                    if (DEBUG)
                         printf("Allowed case %d matches with: >>%s<<  >>%s<<  >>%s<<  >>%s<<\n",
                                iallowed, allowed_torsions[iallowed].A, 
                                          allowed_torsions[iallowed].B,
                                          allowed_torsions[iallowed].C,
                                          allowed_torsions[iallowed].D);

                    ++*p_num_links;  
                    if (!just_count)
                      {
                         p_this_link= p_link_atoms+*p_num_links;
                         p_this_link->start=iatom;
                         p_this_link->end  =iatom2;
     
                         if (DEBUG) printf("DEBUG>> Just got ji link %d for atoms %d %d\n", 1+*p_num_links, p_this_link->start,
                                                                                 p_this_link->end);
                      }
                    break;
                 }
            }
       }
    if (DEBUG) printf("DEBUG>>--------------------------------\n");
 }

return TRUE;
}


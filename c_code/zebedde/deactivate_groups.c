/************************************************************/
/* deactivate_groups                                        */
/* For a bond about to be made mark hydrogens in defined    */
/* groups with the bonding atoms for deletion from the      */
/* active hydrogens list.                                   */
/*                                                          */
/* Started Dave Willock Oct. 2006                           */
/************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "maxima.h"
#include "header.h"

int  deactivate_groups(atom *p_template, int temp_hyd_picked, int temp_bonding_atom,
                       atom *p_molecule, int mol_hyd_picked,  int mol_bonding_atom,
                       int  *p_template_hyd_list, int num_temp_hyds,
                       int  *p_molecule_hyd_list, int num_mol_hyds)
  {

atom *p_c_temp, *p_c_mol, *p_atom, *p_neigh, *p_neigh2, *p_neigh3;
atom *p_c_neigh;

int num_ends, num_done, ineigh, jneigh; 
int num_c_neigh, have_second_methyl;
int ends_list[MAX_ENDS], done_this, iii, iend;
int done_list[MAXTEMPLATE], *p_list, have_done;
int deactivate_rings;

/** define pointers to the atoms giving up their H's ***/

have_done= FALSE;

p_c_temp = p_template+temp_bonding_atom;
p_c_mol  = p_molecule+mol_bonding_atom;

/*****************************************************/
/*** deactivate other atoms on sp3 carbons ***********/
/*****************************************************/
   
if (p_c_temp->num_neigh == 3)
  {

//      printf("Dealing with an sp3 carbon in the current template, neigh of H to loose: %d\n", temp_hyd_picked);
      num_c_neigh=-1;
      for (ineigh=0; ineigh <= p_c_temp->num_neigh; ineigh++)
         {
           p_neigh= p_template+p_c_temp->neighb[ineigh];

           if (strcmp(&(p_neigh->elem[0]),"HA") ==0 
            || strcmp(&(p_neigh->elem[0]),"HB") ==0) 
                {
                  if (p_c_temp->neighb[ineigh] != temp_hyd_picked)
                    {
/****** Find atom in the template H list and flag with -1 for removal ****/
 
                       for (iii=0; iii <= num_temp_hyds; iii++)
                         {
                            p_list=p_template_hyd_list+iii;

                            if (p_c_temp->neighb[ineigh] == *p_list)
                              {
                                have_done=TRUE;
                                *p_list=-1;
                                break;
                              }
                         }
/*** Change atom element back to H ***/
//                     printf(".....reassigning elem %s to H for H atom %d\n", &(p_neigh->elem[0]), p_c_temp->neighb[ineigh]);
                     strcpy(&(p_neigh->elem[0]),"H");
                    }
                }
/****                                ****/
/**** check for second methyl groups ****/
/****                                ****/
            else if (strcmp(&(p_neigh->elem[0]),"C") ==0)
               {
                 num_c_neigh++;
                 p_c_neigh=p_neigh;
               }
         }
/****                                                               ****/
/**** Deal with secondary methyl for Glib's polymer builder example ****/
/**** i.e. for H3C-CH-CH3 deactivate both methyl groups             ****/
/****      at this point the p_c_neigh is the middle C atom.        ****/
/****                                                               ****/
//    printf("....it has %d C neighbours\n", num_c_neigh);
      if (num_c_neigh==0)
        {
//        printf("....checking for second methyl\n");
          for (ineigh=0; ineigh <= p_c_neigh->num_neigh; ineigh++)
             {
               p_neigh2= p_template+p_c_neigh->neighb[ineigh];

               if (strcmp(&(p_neigh2->elem[0]),"C") ==0)
                 {
                   if ( p_neigh2 != p_c_temp )
                     {
/*** see if this is the other methyl group ***/
                        for (jneigh=0; jneigh <= p_neigh2->num_neigh; jneigh++)
                           {
                             p_neigh3= p_template+p_neigh2->neighb[jneigh];

                             if (strcmp(&(p_neigh3->elem[0]),"HA") ==0 
                              || strcmp(&(p_neigh3->elem[0]),"HB") ==0) 
                                {
/****** Find atom in the template H list and flag with -1 for removal ****/
 
                                   for (iii=0; iii <= num_temp_hyds; iii++)
                                     {
                                       p_list=p_template_hyd_list+iii;

                                       if (p_neigh2->neighb[jneigh] == *p_list)
                                         {
                                           have_done=TRUE;
                                           *p_list=-1;
                                           break;
                                         }
                                    }
/*** Change atom element back to H ***/
//                            printf(".....reassigning elem %s to H\n", &(p_neigh3->elem[0])); 
                              strcpy(&(p_neigh3->elem[0]),"H"); 
                                }
                           }
                     }
                 }
             }
        }
  }

if (p_c_mol->num_neigh == 3)
  {
//  printf("Dealing with an sp3 carbon in the molecule being added, as neighbour of H to be deleted %d\n", mol_hyd_picked);
    num_c_neigh=-1;
    for (ineigh=0; ineigh <= p_c_mol->num_neigh; ineigh++)
      {
        p_neigh= p_molecule+p_c_mol->neighb[ineigh];

        if (strcmp(&(p_neigh->elem[0]),"HA")==0
         || strcmp(&(p_neigh->elem[0]),"HB")==0) 
           {
              if (p_c_mol->neighb[ineigh] != mol_hyd_picked)
                {

/****** Find atom in the molecule H list and flag with -1 for removal ****/

                   for (iii=0; iii <= num_mol_hyds; iii++)
                     {
                        p_list=p_molecule_hyd_list+iii;

                        if (p_c_mol->neighb[ineigh] == *p_list)
                          {
                             have_done=TRUE;
                             *p_list=-1;
                             break;
                          }
                     }
/*** Change atom element back to H ***/
//               printf(".....reassigning elem %s to H for atom %d\n", &(p_neigh->elem[0]), p_c_mol->neighb[ineigh]);
                 strcpy(&(p_neigh->elem[0]),"H"); 
                }
           }
/****                                ****/
/**** check for second methyl groups ****/
/****                                ****/
         else if (strcmp(&(p_neigh->elem[0]),"C") ==0)
            {
               num_c_neigh++;
               p_c_neigh=p_neigh;
            }
      }
/****                                                               ****/
/**** Deal with secondary methyl for Glib's polymer builder example ****/
/**** i.e. for H3C-CH-CH3 deactivate both methyl groups             ****/
/****      at this point the p_c_neigh is the middle C atom.        ****/
/****                                                               ****/
//    printf("....it has %d C neighbours\n", num_c_neigh);
      if (num_c_neigh==0)
        {
//        printf("....checking for second methyl\n");
          for (ineigh=0; ineigh <= p_c_neigh->num_neigh; ineigh++)
             {
               p_neigh2= p_molecule+p_c_neigh->neighb[ineigh];

               if (strcmp(&(p_neigh2->elem[0]),"C") ==0)
                 {
                   if ( p_neigh2 != p_c_mol )
                     {
/*** see if this is the other methyl group ***/
                        for (jneigh=0; jneigh <= p_neigh2->num_neigh; jneigh++)
                           {
                             p_neigh3= p_molecule+p_neigh2->neighb[jneigh];

                             if (strcmp(&(p_neigh3->elem[0]),"HA") ==0 
                              || strcmp(&(p_neigh3->elem[0]),"HB") ==0) 
                                {
/****** Find atom in the template H list and flag with -1 for removal ****/
 
                                   for (iii=0; iii <= num_temp_hyds; iii++)
                                     {
                                       p_list=p_template_hyd_list+iii;

                                       if (p_neigh2->neighb[jneigh] == *p_list)
                                         {
                                           have_done=TRUE;
                                           *p_list=-1;
                                           break;
                                         }
                                    }
/*** Change atom element back to H ***/
//                              printf(".....reassigning elem %s to H\n", &(p_neigh3->elem[0]));
                                strcpy(&(p_neigh3->elem[0]),"H");
                                }
                           }
                     }
                 }
             }
        }
//  printf("deactication done...\n");
 }

/*****************************************************/
/*** deactivate other atoms in same benzene ring   ***/
/*** template                                      ***/
/*** Added logical here to turn off ring           ***/
/*** deactivation for rings for Kent problems      ***/
/*** This should be made into an array of logicals ***/
/*** so that they can be defined in the input file.***/
/*** Added here and in next if block, Feb 2013     ***/
/*** by Dave Willock                               ***/
/*****************************************************/
deactivate_rings = FALSE;
if (deactivate_rings && p_c_temp->num_neigh == 2)
  {
     num_ends=0;
     num_done=-1;
     ends_list[0]= temp_bonding_atom;

     while (num_ends >= 0)
        {
           p_atom = p_template+ends_list[0];

           for (ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
             {
               p_neigh=p_template+p_atom->neighb[ineigh];
             
               if (p_neigh->num_neigh == 0)
                 {
                               
/***** If this neighbour is a hydrogen so mark for deactivation ****/
                    if (strcmp(&(p_neigh->elem[0]),"HA")==0
                     || strcmp(&(p_neigh->elem[0]),"HB")==0)
                        {
                           if (p_atom->neighb[ineigh] != temp_hyd_picked)
                             {
/****** Find atom in the template H list and flag with -1 for removal ****/
 
                                for (iii=0; iii <= num_temp_hyds; iii++)
                                  {
                                     p_list=p_template_hyd_list+iii;

                                     if (p_atom->neighb[ineigh] == *p_list)
                                       {
                                         have_done=TRUE;
                                         *p_list=-1;
                                         break;
                                       }
                                  }
/*** Change atom element back to H ***/
                                strcpy(&(p_neigh->elem[0]),"H");
                             }
                           else if (DEBUG)
                             {
                                printf("This is the hydrogen to go, left alone\n");
                             }
                        }
                    else if (DEBUG)
                        {
                           printf(" single neighboured but not HA or HB, %s\n",
                                      p_neigh->elem);
                        }
                  } 

/***** If this neighbour is another sp2 carbon add to ends list ****/
/***** make sure we do not go round the ring more than once!    ****/

                else if (p_neigh->num_neigh == 2 && p_atom->neighb[ineigh] != temp_bonding_atom)
                   {
                     done_this = FALSE;
                     for (iii=0; iii <= num_done; iii++) 
                        {
                          if (done_list[iii] == p_atom->neighb[ineigh])
                             {
                               done_this=TRUE;
                               break;
                             }
                        }

                     if (!done_this)
                        {
                          num_ends++;
                              
                          ends_list[num_ends]= p_atom->neighb[ineigh]; 

                          if (DEBUG) printf(" added to end of list\n");
                        }
                      else if (DEBUG)
                        {
                           printf("Already seen\n");
                        }
                    }
                  else if (DEBUG)
                    {
                      printf(" ignored\n");
                    }
                      
               }
             num_done++;

             if (num_ends >= MAX_ENDS || num_done > MAXTEMPLATE)
              {
                printf("ERROR: Too many ends in deactivate_groups\n");
                exit(0);
              }

              done_list[num_done]=ends_list[0];
/******* shuffle ends_list list to cover the one we have dealt with ****/

              for (iend= 0; iend < num_ends; iend++)
                {
                    ends_list[iend] =  ends_list[iend+1];
                }
                     num_ends--;
           }
       }
   
/*****************************************************/
/*** deactivate other atoms in same benzene ring *****/
/*** molecule                                    *****/
/*****************************************************/
if (deactivate_rings && p_c_mol->num_neigh == 2)
  {
     num_ends=0;
     num_done=-1;
     ends_list[0]= mol_bonding_atom;

     while (num_ends >= 0)
        {
           p_atom = p_molecule+ends_list[0];

           if (DEBUG) printf("Current end: %s\n",p_atom->label);

           for (ineigh=0; ineigh <= p_atom->num_neigh; ineigh++)
             {
               p_neigh=p_molecule+p_atom->neighb[ineigh];
             
               if (p_neigh->num_neigh == 0)
                 {
                               
/***** If this neighbour is a hydrogen so mark for deactivation ****/
                    if (strcmp(&(p_neigh->elem[0]),"HA")==0
                     || strcmp(&(p_neigh->elem[0]),"HB")==0)
                        {
                           if (p_atom->neighb[ineigh] != mol_hyd_picked)
                             {
/****** Find atom in the molecule H list and flag with -1 for removal ****/

                                for (iii=0; iii <= num_mol_hyds; iii++)
                                  {
                                     p_list=p_molecule_hyd_list+iii;

                                     if (p_atom->neighb[ineigh] == *p_list)
                                       { 
                                          have_done=TRUE;
                                          *p_list=-1;
                                          break;
                                       }
                                 }
                             }
                           else if (DEBUG)
                             {
                                printf("This is the hydrogen to go, left alone\n");
                             }
                        }
                    else if (DEBUG)
                        {
                           printf(" single neighboured but not HA or HB, %s\n",
                                      p_neigh->elem);
                        }
                  } 

/***** If this neighbour is another sp2 carbon add to ends list ****/
/***** make sure we do not go round the ring more than once!    ****/

                else if (p_neigh->num_neigh == 2 && p_atom->neighb[ineigh] != mol_bonding_atom )
                   {
                     done_this = FALSE;
                     for (iii=0; iii <= num_done; iii++) 
                        {
                          if (done_list[iii] == p_atom->neighb[ineigh])
                             {
                               done_this=TRUE;
                               break;
                             }
                        }

                     if (!done_this)
                        {
                          num_ends++;
                              
                          ends_list[num_ends]= p_atom->neighb[ineigh]; 

                          if (DEBUG) printf(" added to end of list\n");
                        }
                      else if (DEBUG)
                        {
                           printf("Already seen\n");
                        }
                    }
                  else if (DEBUG)
                    {
                      printf(" ignored\n");
                    }
                      
               }
             num_done++;

             if (num_ends >= MAX_ENDS || num_done > MAXTEMPLATE)
              {
                printf("ERROR: Too many ends in deactivate_groups\n");
                exit(0);
              }

              done_list[num_done]=ends_list[0];
/******* shuffle ends_list list to cover the one we have dealt with ****/

              for (iend= 0; iend < num_ends; iend++)
                {
                    ends_list[iend] =  ends_list[iend+1];
                }
                     num_ends--;
           }
       }
    return (have_done);
 }

/************************************************************/
/* join_frag_to_temp                                        */
/* add a molecule to the template                           */
/*                                                          */
/* Started DJW 11/7/95                                      */
/* Updated for malloc and realloc Dave Willock March 2013   */
/************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

/*****DEBUG******/
void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp);

/***END DEBUG****/

int pick_frm_wgt_list( int sum_weights, int *p_weights, int *p_chosen_weight,
                       int num_in_list );

int forbid_bond( char *p_type1,          char *p_type2,
                 bond *p_forbidden_bond);

void align_bonds(atom *p_template, atom *p_new_bit, atom *p_new_bonding_atom,
                 atom *p_new_hyd_picked, atom *p_temp_bonding_atom,
                 atom *p_temp_hyd_picked, int num_template_atoms,
                 int num_new_atoms);

int shuffle_atom_list(atom *p_molecule,
                      int num_atoms, int hyd_index1, int hyd_index2);

int fix_potentials( atom *p_molecule, int atom_number, int num_atoms);

void shuffle_neighbours(atom *p_template, int p_place_for_one,
                        int place_for_two, int num_template_atoms);

void shuffle_links(links *p_link_atoms, int place_for_one,
                        int place_for_two, int num_links);

int shuffle_hyd_list(int *p_list, int place_for_one,
                     int place_for_two, int num_in_list,
                     int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights,
                     int have_deactivated);

int  deactivate_groups(atom *p_template, int temp_hyd_picked, int temp_bonding_atom,
                       atom *p_molecule, int mol_hyd_picked,  int mol_bonding_atom,
                       int  *p_template_hyd_list, int num_temp_hyds,
                       int  *p_molecule_hyd_list, int num_mol_hyds);

int index_double_pot(int index_1, int index_2, int num_types, int *p_1_then_2);

int check_allowed_torsions(atom *p_molecule, int num_atoms, int iatom, int iatom2,
                           int ineigh1, int ineigh2, links *p_link_atoms,
                           int *p_num_links, dihedral *p_allowed_torsions, 
                           int num_allowed_torsions, int just_count);

/********************************************************************************/
/*** p_template is the current guest that has been selected for growth **********/
/*** p_molecule is the fragment to be added.                           **********/
/********************************************************************************/

int join_frag_to_temp(atom *p_template, int *p_num_template_atoms,
                      atom *p_molecule, int num_molecule_atoms,
                      int *p_place_for_one, int *p_place_for_two,
                      int *p_template_hyd_list, int *p_num_template_hyds,
                      int *p_molecule_hyd_list, int num_molecule_hyds,
                      bond *p_forbidden_bond, int *p_temp_hyd_weights,
                      int *p_sum_temp_hyd_weights,
                      int increment_hydrogen_weight, int have_AB)
{
#include "header.h"
int mol_hyd_picked;

int allowed_mol_hyds[MAXTEMPLATE], have_mol_hyds, sum_mol_hyd_wgts;

int iflag, temp_hyd_picked, ineigh, num_old_template_atoms,iatom;
int num_old_template_hyds, ihyd;
int mol_bonding_atom, temp_bonding_atom;
int *p_new_list, *p_mol_hyd_list;
int can_bond, good_AB,dummy;
int temp_parent_hyd_weight, mol_parent_hyd_weight;
int *p_new_temp_hyd_weights;
int pot_index_t, pot_index_f, pot_index_th, pot_index_fh;
int index, f_then_fh, t_then_th, t_then_f;
int have_deactivated, idebug, idebug2;
int just_count, mol_max_index, debug_bail;

atom *p_mol_hyd_picked, *p_mol_bonding_atom, *p_temp_hyd_picked;
atom *p_temp_bonding_atom, *p_new_atom, *p_old_atom, *p_new_bit;
atom *p_debug_atom;

have_deactivated= FALSE;

/**** define the maximum index for the new molecule to add ***/
mol_max_index= num_molecule_atoms-1;

DEBUG=FALSE;
/***************************** pick a template valance ******************************/
if (DEBUG)
  {
     debug_bail=0;
     printf("In join_frag_to_temp picking template hydrogen\n");
     printf("number of hydrogens in template = %d\n",*p_num_template_hyds);

     for (ihyd=0; ihyd <= *p_num_template_hyds; ihyd++)
       {
          p_new_atom = p_template+ *(p_template_hyd_list+ihyd);
          p_new_temp_hyd_weights= p_temp_hyd_weights+ihyd;
          printf("\n%d => %s (elem %s) weight: %d ", ihyd, p_new_atom->label, 
                                            p_new_atom->elem, *p_new_temp_hyd_weights );
       }
     printf("\nSum H-weights: %d\n",  *p_sum_temp_hyd_weights);
     printf("\n\n");
     printf("number of hydrogens in fragment = %d\n", num_molecule_hyds);

     printf("\nFragment being added looks like:\n");
     p_debug_atom=p_molecule;
     for (index=0; index <= mol_max_index; index++)
       {
         printf("%d  ->  %s (%s)\n", index, p_debug_atom->label, p_debug_atom->elem);
         p_debug_atom++;
       }
  }

/********** sort out which hydrogens are allowed ************************************/

have_mol_hyds= FALSE;
sum_mol_hyd_wgts=0;

/******************************************************************************************/
/****** Need to calculate the template hydrogen weights according to user supplied ********/
/****** rules now!!                                                                ********/
/******************************************************************************************/
/****** Initially flag all H atoms as possible connection points **************************/
/******************************************************************************************/

for (iflag=0; iflag <= num_molecule_hyds; iflag++)
  {
    allowed_mol_hyds[iflag] =1;
    have_mol_hyds=TRUE;
    sum_mol_hyd_wgts += allowed_mol_hyds[iflag];
  }
	
if (*p_sum_temp_hyd_weights < 0)
   {
     printf("No template hydrogens avaiable for continued building template! Oh Err\n");
     exit (0);
   }
if (!have_mol_hyds)
   {
     printf("This fragment has no available hydrogens!!!!! Oh Err\n");
     exit (0);
   }

/*******************************************************************************************/
/************** Only accept bonds that are not forbidden ***********************************/
/*******************************************************************************************/

can_bond = FALSE;
while (!can_bond)
  {

/*******************************************************************************************/
/************************************* pick a template valence *****************************/
/*** The hyd_list and hyd_weights list are indexed in the same way and so picking from *****/
/*** the weights list gives a valid index for both lists ***********************************/
/*******************************************************************************************/

   if (DEBUG) 
     {
       printf("SI DEBUG>> Choosing temp_hyd_picked using pick_frm_wgt_list\n");
       printf("SI DEBUG>> sending sum_wgts=%d num=%d \n", *p_sum_temp_hyd_weights, *p_num_template_hyds);
       printf("SI DEBUG>> Weights:\n");
       for (idebug=0; idebug<=*p_num_template_hyds; idebug++)
          {
             idebug2=*(p_template_hyd_list + idebug);
             printf("%d: %d %s (elem:%s) wgt: %d\n", idebug, idebug2,  
                                                     (p_template+idebug2)->label, 
                                                     (p_template+idebug2)->elem, 
                                                     *(p_temp_hyd_weights+idebug)); 
          }
     }

   temp_hyd_picked   = pick_frm_wgt_list(*p_sum_temp_hyd_weights, p_temp_hyd_weights,
                                          &temp_parent_hyd_weight, *p_num_template_hyds );

//   printf("DEBUG Chosen %d\n", temp_hyd_picked);

   temp_hyd_picked   = *(p_template_hyd_list +temp_hyd_picked);
   p_temp_hyd_picked = p_template+ temp_hyd_picked;

//   printf("DEBUG Corresponds to atom %d\n", temp_hyd_picked);
//   printf("DEBUG which is %s (elem: %s)\n", p_temp_hyd_picked->label, p_temp_hyd_picked->elem);

//   printf("DEBUG 1. first hyd weight: %d\n", *p_temp_hyd_weights);

/*******************************************************************************************/
/************** find the template neighbour of the picked hydrogen *************************/
/*******************************************************************************************/

   temp_bonding_atom   = p_temp_hyd_picked->neighb[0];
   p_temp_bonding_atom = p_template+ temp_bonding_atom; 

/*******************************************************************************************/
/************************************* pick a fragment valence *****************************/
/** allowed_mol_hyds has been set evenly to 1 above and so acts as a weights list here *****/
/** in which all H atoms have the same probability of being chosen.                    *****/
/*******************************************************************************************/

   if (DEBUG) 
     {
       printf("SI DEBUG>> Choosing mol_hyd_picked using pick_frm_wgt_list\n");
       printf("SI DEBUG>> sending sum_wgts=%d num=%d \n", sum_mol_hyd_wgts, num_molecule_hyds);
       printf("SI DEBUG>> List of allowed molecule hydrogens:\n");
  
       for (idebug=0; idebug <= num_molecule_hyds; idebug++)
         {
           printf("%d\n", allowed_mol_hyds[idebug]); 
         }
     }

   mol_hyd_picked = pick_frm_wgt_list(sum_mol_hyd_wgts, &allowed_mol_hyds[0], 
                                               &mol_parent_hyd_weight, num_molecule_hyds);
   mol_hyd_picked = *(p_molecule_hyd_list +mol_hyd_picked);
   p_mol_hyd_picked = p_molecule + mol_hyd_picked; 

   mol_bonding_atom   = p_mol_hyd_picked->neighb[0];
   p_mol_bonding_atom = p_molecule+ mol_bonding_atom; 

   can_bond= !forbid_bond(p_temp_bonding_atom->elem,
                          p_mol_bonding_atom->elem,
                          p_forbidden_bond);

//   printf("DEBUG 2. first hyd weight: %d\n", *p_temp_hyd_weights);
/*********************************************************************************************/
/*** If using HA..HB hydrogen element types only allow HA from 1 and HB from other ***********/
/*** added October 2006 Dave Willock                                               ***********/
/*********************************************************************************************/
   if (have_AB)
     {
        good_AB =    strcmp(p_temp_hyd_picked->elem, "HA") == 0 
                  && strcmp(p_mol_hyd_picked->elem, "HB")  == 0; 

        good_AB = good_AB || (strcmp(p_temp_hyd_picked->elem, "HB") == 0 
                          && strcmp(p_mol_hyd_picked->elem, "HA")   == 0 );

        if (DEBUG) printf("DEBUG>> have_AB so testing temp_hyd %s and mol_hyd %s\n",
                              p_temp_hyd_picked->elem, p_mol_hyd_picked->elem);

        can_bond = can_bond && good_AB;

        if (can_bond) 
          {
            if (DEBUG)
              {
                 printf("That will do for a new bond\n");

/****** Flag other atoms for de-activation as Cl *****/
/****** firstly for atoms bonded to same C       *****/

                  printf("temp_hyd_picked: %d %s (elem %s)\n",temp_hyd_picked,
                              p_temp_hyd_picked->label, p_temp_hyd_picked->elem);

                  printf("mol_hyd_picked: %d %s (elem %s)\n",mol_hyd_picked,
                              p_mol_hyd_picked->label, p_mol_hyd_picked->elem);
              }

            have_deactivated= deactivate_groups(p_template, temp_hyd_picked, temp_bonding_atom,
                                                p_molecule, mol_hyd_picked,  mol_bonding_atom,
                                                p_template_hyd_list, *p_num_template_hyds,
                                                p_molecule_hyd_list, num_molecule_hyds);
         }
     }

   if (DEBUG)
     {
        printf("DEBUG 3. first hyd weight: %d\n", *p_temp_hyd_weights);
        printf("For Fragment sent:\n sum_weights= %d\n weights are:\n",sum_mol_hyd_wgts);
        for (ihyd=0; ihyd <= num_molecule_hyds; ihyd++)
                                 printf("%d => %d\n",ihyd, allowed_mol_hyds[ihyd]); 
        printf("DEBUG >> Picked hydrogen from template label=%s\n", 
                                                        p_temp_hyd_picked->label);
        printf("Neighbour of chosen template hydrogen = %s\n\n",
                                                      p_temp_bonding_atom->label); 
        printf("picked hydrogen from fragment = %d label=%s\n ",
                                         mol_hyd_picked, p_mol_hyd_picked->label);
        printf("Neighbour of chosen molecule hydrogen = %s\n\n",
                                                       p_mol_bonding_atom->label); 

        if ( strncmp(p_temp_hyd_picked->label, "H",1) != 0 )
          {
             printf("ERROR: guest H picked looks strange?????\n");
             printf("ERROR: exiting exiting..................\n");
             exit(0);
          }
        if ( strncmp(p_mol_hyd_picked->label, "H",1) != 0 )
          {
             printf("ERROR: fragment H picked looks strange?????\n");
             printf("ERROR: exiting exiting..................\n");
             exit(0);
          }

        if (can_bond)
          { 
             printf("This bond will be acceptable\n");
             debug_bail=0;
          }
        else
          {
             printf("This bond will NOT be acceptable\n");
             debug_bail++;
//             if (debug_bail > 10) exit(0);
          }
     }
  }

num_old_template_atoms = *p_num_template_atoms;
p_new_atom= p_template+ *p_num_template_atoms;
p_old_atom= p_molecule-1;

if (DEBUG) printf("DEBUG>> num_molecule_atoms = %d max index = %d\n", num_molecule_atoms, mol_max_index);

/***** p_new_atom points to the end of the current guest list *****/
/***** p_old_atom points to the list of atoms to be added     *****/

for (iatom= 0; iatom <= mol_max_index; iatom++)
  {
     p_new_atom++;
     p_old_atom++;

    *p_new_atom= *p_old_atom;

/*** make good the molecule index of the new atoms ***/
/*** added March 07 Dave Willock                   ***/

    p_new_atom->mol = p_temp_bonding_atom->mol;

     if (p_old_atom == p_mol_hyd_picked ) 
       {
          p_mol_hyd_picked   = p_new_atom;
          mol_hyd_picked     = *p_num_template_atoms+ iatom +1;
       }
     if (p_old_atom == p_mol_bonding_atom) 
       {
          p_mol_bonding_atom = p_new_atom;
          mol_bonding_atom   = *p_num_template_atoms+ iatom +1;
       }
  }

if (DEBUG)
         printf("Leaving loop num_molecule_atoms = %d, current number of template atoms = %d\n", 
                                                               num_molecule_atoms, *p_num_template_atoms);

*p_num_template_atoms+= num_molecule_atoms;

if (DEBUG)
  {
    printf("\nBefore joining molecule to template num atoms = %d, after = %d\n",
                   num_old_template_atoms, *p_num_template_atoms);

    printf("\nReindexed _mol_hyd_picked= %s, mol_bonding_atom = %s\n\n",
                     p_mol_hyd_picked->label, p_mol_bonding_atom->label);
  }

/*********************************************************************************/
/*************************Redefine the neighbours for the new bit ****************/
/*********************************************************************************/

p_new_atom= p_template+num_old_template_atoms;

for (iatom= 0; iatom <= mol_max_index; iatom++)
  {
     p_new_atom++;

     for (ineigh = 0; ineigh <= p_new_atom->num_neigh; ineigh++)
        {
           p_new_atom->neighb[ineigh]+= num_old_template_atoms+1;
        }
  }
 

//if (DEBUG)
//  {
//    printf("\n\nNeighbours immediately after joining new bit to template list:\n\n");
//    print_neighbours( p_template, *p_num_template_atoms, stdout);
//  }

/****************Add all the molecules hydrogens to the templates*****************/

num_old_template_hyds= *p_num_template_hyds;
p_mol_hyd_list= p_molecule_hyd_list-1;

//printf("In join have %d old template hydrogens and %d from the new part\n", num_old_template_hyds, num_molecule_hyds);

p_new_temp_hyd_weights= p_temp_hyd_weights+ num_old_template_hyds;
p_new_list= p_template_hyd_list+num_old_template_hyds;

if (DEBUG) printf("DEBUG>> Parent H weight = %d\n", temp_parent_hyd_weight);

for (ihyd= 0; ihyd <= num_molecule_hyds; ihyd++)
  {
    p_mol_hyd_list++;
    p_new_list++;

/***** Carry over deactivation flags *****/ 
    if (*p_mol_hyd_list == -1)
      {
        *p_new_list= *p_mol_hyd_list;
      }
    else
      {
        *p_new_list= (*p_mol_hyd_list) + num_old_template_atoms+1;
//      printf("New H index = %d, %s (elem %s)\n", *p_new_list, 
//                                                (p_template+*p_new_list)->label,
//                                                (p_template+*p_new_list)->elem);
      }

/****** Alter hydrogen weights **************************************************/
//   printf("DEBUG 4. first hyd weight: %d\n", *p_temp_hyd_weights);

    p_new_temp_hyd_weights++;
    (*p_new_temp_hyd_weights) = temp_parent_hyd_weight + 
                                         increment_hydrogen_weight;

    if (*p_new_temp_hyd_weights < 0) *p_new_temp_hyd_weights=0;

    *p_sum_temp_hyd_weights += *p_new_temp_hyd_weights;
  }

*p_num_template_hyds += num_molecule_hyds+1;

if (DEBUG)
  {
    printf("Currently have %d hydrogens with indices and weights:\n",
                                                         *p_num_template_hyds);
    p_new_list= p_template_hyd_list;
    p_new_temp_hyd_weights= p_temp_hyd_weights;
    for (ihyd= 0; ihyd <= *p_num_template_hyds; ihyd++)
      {
        if (*p_new_list > 0)
          {
        printf("%d hyd index= %d, label=%s (elem=%s), weight = %d  ",
                            ihyd, *p_new_list, (p_template+*p_new_list)->label, 
                            (p_template+*p_new_list)->elem, 
                            *p_new_temp_hyd_weights);
          }
        else
          {
        printf("%d hyd index= %d, so will be deactivated ",
                            ihyd, *p_new_list);
          }

        if   (*p_new_list == temp_hyd_picked
           || *p_new_list == mol_hyd_picked)
                printf("One that will be deleted on bonding\n");
        else
                printf("\n");

        p_new_list++;
        p_new_temp_hyd_weights++;
      }
  }

/****************************** align the bonds **********************************/

p_new_bit= p_template+num_old_template_atoms+1;

align_bonds(p_template, p_new_bit, p_mol_bonding_atom, p_mol_hyd_picked,
            p_temp_bonding_atom, p_temp_hyd_picked, *p_num_template_atoms,
            mol_max_index);

/************** Adjust neighbours information accordingly ***********************/


 for (ineigh=0; ineigh <= p_temp_bonding_atom->num_neigh; ineigh++)
    {
      if (p_temp_bonding_atom->neighb[ineigh] == temp_hyd_picked)
        {
           p_temp_bonding_atom->neighb[ineigh]= mol_bonding_atom;
        }
    }

 for (ineigh=0; ineigh <= p_mol_bonding_atom->num_neigh; ineigh++)
    {
      if (p_mol_bonding_atom->neighb[ineigh] == mol_hyd_picked)
        {
           p_mol_bonding_atom->neighb[ineigh]= temp_bonding_atom;
        }
    }

/**************** Adjust potentials if AMBER forcefield in use *******************/

 if (DEBUG) printf ("version= %d, AMBER= %d, newly bonded atoms: %s (%s) and %s (%s)\n",
               pot_info.version, AMBER, p_mol_bonding_atom->label, p_mol_bonding_atom->pot,
                                        p_temp_bonding_atom->label, p_temp_bonding_atom->pot);

 if (pot_info.version == AMBER)
    {

/**************** Sort out atom charges using bond increments  *******************/
/**************** take off hydrogen delta first                *******************/

      pot_index_t  = p_temp_bonding_atom->bi_list;
      pot_index_f  = p_mol_bonding_atom->bi_list;
      pot_index_th = p_temp_hyd_picked->bi_list;
      pot_index_fh = p_mol_hyd_picked->bi_list;

if (DEBUG)
  {
printf ("In join_frag_to_temp about to make bond between temp:%s (pot %d) and frag:%s (pot %d)\n", 
                 p_temp_bonding_atom->elem, pot_index_t, p_mol_bonding_atom->elem, pot_index_f );
printf ("Loosing element %s (%d) from temp and %s (%d) from frag\n", 
                 p_temp_hyd_picked->elem, pot_index_th, p_mol_hyd_picked->elem, pot_index_fh);
  }

      index= index_double_pot(pot_index_f, pot_index_fh, num_bond_inc_types, &f_then_fh);

      if (f_then_fh)
        {
          p_mol_bonding_atom->part_chge -= bond_increments[index].a;

          if (DEBUG)
                  printf("to get rid of hydrogen charge for mol removing %10.6f\n", 
                                                                   bond_increments[index].a);
        }
      else
        {
          p_mol_bonding_atom->part_chge -= bond_increments[index].b;

          if (DEBUG)
                  printf("to get rid of hydrogen charge for mol removing %10.6f\n", 
                                                                   bond_increments[index].b);
        }

      
      index= index_double_pot(pot_index_t, pot_index_th, num_bond_inc_types, &t_then_th);
      if (DEBUG) printf("For template and its hydrogen using index %d\n",index);

      if (t_then_th)
        {
          p_temp_bonding_atom->part_chge -= bond_increments[index].a;

          if (DEBUG)
                  printf("to get rid of hydrogen charge for temp removing %10.6f\n", 
                                                                  bond_increments[index].a);
        }
      else
        {
          p_temp_bonding_atom->part_chge -= bond_increments[index].b;

          if (DEBUG)
                  printf("to get rid of hydrogen charge for temp removing %10.6f\n", 
                                                                  bond_increments[index].b);
        }
   
/**************************************************************************************************/
/****** Now change potential type of the nitrogen in the backbone *********************************/
/**************************************************************************************************/

      if (strcmp(p_mol_bonding_atom->pot, "NT") == 0 ) 
        {
           strcpy(p_mol_bonding_atom->pot, "N");
           for (index=0; index<=num_h_pot_types; index++) 
             {
               if (strcmp( h_pot_types[index].name, p_mol_bonding_atom->pot ) == 0 )
                 {
                    p_mol_bonding_atom->hb_list = index;
                 }
             }
           for (index=0; index<=num_bond_inc_types; index++) 
             {
               if (strcmp( bond_inc_types[index].name, p_mol_bonding_atom->pot ) == 0 )
                 {
                    p_mol_bonding_atom->bi_list = index;
                 }
             }
         }

      if (strcmp(p_temp_bonding_atom->pot, "NT") == 0 ) 
        {
           strcpy(p_temp_bonding_atom->pot, "N");
           for (index=0; index<=num_h_pot_types; index++) 
             {
               if (strcmp( h_pot_types[index].name, p_temp_bonding_atom->pot ) == 0 )
                 {
                    p_temp_bonding_atom->hb_list = index;
                 }
             }
           for (index=0; index<=num_bond_inc_types; index++)     
             {
               if (strcmp( bond_inc_types[index].name, p_temp_bonding_atom->pot ) == 0 )
                 {
                    p_temp_bonding_atom->bi_list = index;
                 }
             }
        }

/**************************************************************************************************/
/****** Finally do delta of new bond **************************************************************/
/**************************************************************************************************/
      pot_index_t  = p_temp_bonding_atom->bi_list;
      pot_index_f  = p_mol_bonding_atom->bi_list;

      index= index_double_pot(pot_index_t, pot_index_f, num_bond_inc_types, &t_then_f);
      if (DEBUG) printf("For template and fragment using index %d\n",index);

      if (t_then_f)
        {
          p_temp_bonding_atom->part_chge += bond_increments[index].a;
          p_mol_bonding_atom->part_chge += bond_increments[index].b;
        
          if (DEBUG)
                  printf("to account for new bond temp gets %10.6f frag %10.6f\n", 
                                                 bond_increments[index].a, bond_increments[index].b);
        }
      else
        {
          p_temp_bonding_atom->part_chge += bond_increments[index].b;
          p_mol_bonding_atom->part_chge += bond_increments[index].a;
         
          if (DEBUG)
                  printf("to account for new bond temp gets %10.6f frag %10.6f\n", 
                                                 bond_increments[index].b, bond_increments[index].a);
        }
    }

/*************** delete template and fragment valence hydrogens *****************/

if (DEBUG)
  {
     printf("\nentering atom list shuffle loop\n");
     printf("num_template_atoms before shuffle = %d\n",*p_num_template_atoms);
  }

/******************************************************************************/
/******************update lists ***********************************************/
/******************************************************************************/

 *p_num_template_atoms= shuffle_atom_list(p_template,
                                          *p_num_template_atoms, 
                                          temp_hyd_picked, mol_hyd_picked);

if (DEBUG)
  {
     printf("num_template_atoms after hyd deletion = %d\n",
                                               *p_num_template_atoms);
  }


/*** Deal with link atoms in assemble_intra after build, not here ******/
/* link_atoms[num_links].start = temp_bonding_atom;*/
/* link_atoms[num_links].end   = mol_bonding_atom; */
/*******************shuffle links to cover up for lst hydrogens ****************/
/* shuffle_links(&link_atoms[0], temp_hyd_picked, mol_hyd_picked, num_links); */

/*******************shuffle neighbours to cover up for lst hydrogens ***********/

 shuffle_neighbours(p_template, temp_hyd_picked, mol_hyd_picked,
                                                *p_num_template_atoms);

/*******************shuffle hydrogen list to cover up for lst hydrogens ***********/
/*******************and to remove deactivated H atoms if appropriate    ***********/

 *p_num_template_hyds= shuffle_hyd_list(p_template_hyd_list, temp_hyd_picked, 
                                        mol_hyd_picked, *p_num_template_hyds, 
                                        p_temp_hyd_weights, p_sum_temp_hyd_weights,
                                        have_deactivated);

  if (DEBUG)
   {
      printf("After shuffle %d hydrogens with indices and weights:\n",
                                                         *p_num_template_hyds);
      p_new_list= p_template_hyd_list;
      p_new_temp_hyd_weights= p_temp_hyd_weights;
      for (ihyd= 0; ihyd <= *p_num_template_hyds; ihyd++)
        {
          printf("%d hyd index= %d, label=%s (elem=%s), weight = %d\n",
                              ihyd, *p_new_list, (p_template+*p_new_list)->label,
                             (p_template+*p_new_list)->elem, 
                             *p_new_temp_hyd_weights);
          p_new_list++;
          p_new_temp_hyd_weights++;
        }

      printf("Total weight: %d\n", *p_sum_temp_hyd_weights);
   }

/*******************Patch potentials***************************************/
/*******************moved to after neighbours to get connectivity right ***/
/*******************before patching potentials.                         ***/

 dummy = fix_potentials( p_template, temp_bonding_atom, *p_num_template_atoms);
 dummy = fix_potentials( p_template, mol_bonding_atom, *p_num_template_atoms);
   
/*******************Alter tether atom index if required ***************************/

if (have_tethers && tether.index_b > temp_hyd_picked ) (tether.index_b)--;

if (DEBUG)
  {
     printf("\nHydrogens in template at start = %d, after addition = %d\n\n",
                             num_old_template_hyds,*p_num_template_hyds);

     for (ihyd=0; ihyd <= *p_num_template_hyds; ihyd++)
       {
          p_new_atom = p_template+ *(p_template_hyd_list+ihyd);
          p_new_temp_hyd_weights= p_temp_hyd_weights+ihyd;
          printf("\n%d => %s weight: %d ", ihyd, p_new_atom->label, *p_new_temp_hyd_weights );
       }
     printf("\nNew sum of weights = %d",*p_sum_temp_hyd_weights);
     printf("\n\n");
     print_neighbours( p_template, *p_num_template_atoms, stdout);
  }

DEBUG=FALSE;
return(0); 
}


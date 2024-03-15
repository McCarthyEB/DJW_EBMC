/************************************************************/
/* ring_maker.c                                             */
/*  New ring maker. Attempts to include all sizes of ring   */
/*  Uses some of Daves original code	                    */
/*                                                          */
/* Started AJWL 4/09                                        */
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

double real_random(int done);

void print_neighbours( atom *p_molecule, int num_atoms, FILE *fp);

void shuffle_neighbours(atom *p_template, int place_for_one,
                        int place_for_two, int num_template_atoms);

void shuffle_links(links *p_link_atoms, int place_for_one,
                        int place_for_two, int num_links);

int shuffle_hyd_list(int *p_list, int place_for_one,
                     int place_for_two, int num_in_list,
                     int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights,
                     int have_deactivated);

int shuffle_atom_list(atom *p_molecule, int num_atoms, 
                                       int hyd_index1, int hyd_index2);
 
int atom_has_hyds(atom *p_molecule, int atom);

int forbid_bond( char *p_type1,          char *p_type2,
                 bond *p_forbidden_bond);

int neighbour_order(atom *p_molecule, int num_atoms, int atom1, int atom2);

int fix_potentials( atom *p_molecule, int atom_number, int num_atoms);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc);

int get_random_int(int min, int max);

int ring_maker(atom *p_molecule, int *p_num_atoms, bond *p_forbidden_bond,
                int *p_template_hyd_list, int *p_num_template_hyds,
                int *p_temp_hyd_weights, int *p_sum_temp_hyd_weights)
{
#include "header.h"
int i=0;
int j=0;
int num_pairs=0;
int can_bond;
int near_neighs[128][2];
int atom1,atom2;
int allowed_bond1,allowed_bond2;
int order;
int rand_no;
int ineigh1, ineigh2, index1, index2;
int hyd_index1=-1, hyd_index2=-1, neigh_of_1=-1, neigh_of_2=-1;
int this_pair;
int dummy, itemp;
int low_hyd, hi_hyd;
int have_made_ring;
int start=0;

char dummy2[40];
char dummy3[40];

double r,r2;
double least_hyd_sep_sqrd;

atom *p_atom1, *p_atom2, *p_neigh1, *p_neigh2;


/*************************************************************/
/**** Loop over all atoms to find atoms which are within *****/
/**** cutoff and are able to form a ring		 *****/
/*************************************************************/
have_made_ring = FALSE;

for (atom1=0; atom1 < *p_num_atoms; atom1++)
	{
	allowed_bond1 = FALSE;
	p_atom1= p_molecule+atom1;
	
	if (atom_has_hyds(p_molecule, atom1) && p_atom1->num_neigh > 0) allowed_bond1 = TRUE;
		
	for (atom2=atom1+1; atom2 <= *p_num_atoms; atom2++)
		{
                p_atom2= p_molecule+atom2;
		can_bond = FALSE;
		allowed_bond2 = FALSE;
		if (atom_has_hyds(p_molecule, atom2) && p_atom2->num_neigh > 0) allowed_bond2 = TRUE;
		
		r2=atom_separation_squared(p_atom1, p_atom2, pbc);
		r=sqrt(r2);
		
		can_bond = !forbid_bond(p_atom1->elem,
                                p_atom2->elem,
                                p_forbidden_bond);

		order= neighbour_order(p_molecule, *p_num_atoms, atom1, atom2);
		
		if (r2 <= ring_ctf_2 && can_bond && allowed_bond1 && allowed_bond2 && order >=4) //store that pair of atoms
			{
			printf("DB>> Atom 1: %i, atom 2: %i\n",atom1,atom2);
			near_neighs[i][0]=atom1;
			near_neighs[j][1]=atom2;
			i++;
			j++;
			num_pairs++;
			}
		}
	}

if (num_pairs == 0)
	{
	printf("No pairs of atoms available to make rings. Returning\n");
	return FALSE;
	}
/****************************************************************/
/***** Now go through those pairs and see what size rings	*/
/***** we can form rings				        */
/****************************************************************/
/*i=0;
j=0;

for (this_pair=0;this_pair <= num_pairs; this_pair++)
	{
	atom1 = near_neighs[i][0];
	atom2 = near_neighs[j][1];
	order= neighbour_order(p_molecule, *p_num_atoms, atom1, atom2);
	printf("Order for pair %i %i is %i\n", atom1,atom2,order);
	i++;
	j++;
	}
*/
rand_no = get_random_int(0, num_pairs);

printf("The Random number is %i\n",rand_no);

i=0;
j=0;

for (this_pair=0;this_pair < num_pairs; this_pair++)
        {
        atom1 = near_neighs[i][0];
        atom2 = near_neighs[j][1];
        
	if (this_pair == rand_no)
		{
        	printf("Pair to bind will be %i %i \n", atom1,atom2);
        	break;
		}
	
	else
		{
		printf("Not this pair %i %i \n", atom1, atom2);
		}
	i++;
        j++;
        }

/*********************************************************************/
/**** Rest of code is same as written by DJW *************************/
/*********************************************************************/


least_hyd_sep_sqrd = -1.0;

p_atom1 = p_molecule+atom1;
p_atom2 = p_molecule+atom2;

for (ineigh1 = 0; ineigh1 <= p_atom1->num_neigh; ineigh1++)
	{
	index1= p_atom1->neighb[ineigh1];
	p_neigh1= p_molecule + index1;

	if ( p_neigh1->elem[0] == 'H' && p_neigh1->elem[1] == '\0')
		{
		for (ineigh2 = 0; ineigh2 <= p_atom2->num_neigh; ineigh2++)
			{
			index2= p_atom2->neighb[ineigh2];
			p_neigh2= p_molecule + index2;

			if ( p_neigh2->elem[0] == 'H' && p_neigh2->elem[1] == '\0')
				{
				r2 = atom_separation_squared(p_neigh1, p_neigh2, pbc);
				if (r2 < least_hyd_sep_sqrd || least_hyd_sep_sqrd < 0)
					{
					least_hyd_sep_sqrd = r2;
					hyd_index1= index1;
					hyd_index2= index2;
					neigh_of_1= ineigh1;
					neigh_of_2= ineigh2;
					}
				}
			}
		}
	}

/***************************************************************************/
/****** Now have hyd_index1 and hyd_index2 as the atom positions     *******/
/****** to be removed from the list on bonding. Their postions in    *******/
/****** the atom1 and atom2 neighbours lists are given by neigh_of_1 *******/
/****** and neigh_of_2 respectivily.                                 *******/
/****** ERROR checking added December 2006, Dave Willock             *******/
/***************************************************************************/

if (   hyd_index1 == -1 || hyd_index2 == -1 || neigh_of_1 == -1 || neigh_of_2 == -1)
	{
        printf("ERROR >> ring_maker function failed to select a valid atom set\n");
        printf("         Please report this error to the authors.\n");
        exit(0);
        }
else 
	{
	have_made_ring = TRUE;
	}

p_atom1->neighb[neigh_of_1]= atom2;
p_atom2->neighb[neigh_of_2]= atom1;

/***************************************************************************/
/****** Shuffle indicies to accomadate the loss of the hydrogens     *******/
/***************************************************************************/
low_hyd= hyd_index1;
hi_hyd = hyd_index2;

if (low_hyd > hi_hyd)
	{
	itemp= low_hyd;
	low_hyd= hi_hyd;
	hi_hyd= itemp;
	}

/******************* Patch up potentials *******************************************/

dummy = fix_potentials( p_molecule, atom1, *p_num_atoms);
dummy = fix_potentials( p_molecule, atom2, *p_num_atoms);

*p_num_atoms= shuffle_atom_list(p_molecule, *p_num_atoms, hyd_index1, hyd_index2);

/*******************shuffle neighbours to cover up for lost hydrogens ***********/

shuffle_neighbours(p_molecule, low_hyd, hi_hyd, *p_num_atoms);

/*******************shuffle links to cover up for lost hydrogens ****************/

// shuffle_links(&link_atoms[0], low_hyd, hi_hyd, num_links);

/*******************shuffle hydrogen list to cover up for lost hydrogens ***********/

*p_num_template_hyds= shuffle_hyd_list(p_template_hyd_list, low_hyd, hi_hyd,
                                                           *p_num_template_hyds,
                                                            p_temp_hyd_weights,
                                                            p_sum_temp_hyd_weights,
                                                            FALSE);

print_neighbours(p_molecule, *p_num_atoms, stdout);

if (have_made_ring == TRUE) return TRUE;
else return FALSE;

}		

/**** get_random_int is designed to return a random number ****/
/**** between min and max.                                 ****/

int get_random_int(int min, int max)
{
#include "header.h"
int random_no;

random_no = (real_random(0) * (max - min + 1) + min);

return (random_no);
}























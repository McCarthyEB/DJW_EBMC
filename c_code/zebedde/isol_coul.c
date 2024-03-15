/****************************************************************************************/
/***** Simple summation electrostatics for clusters and removing self interaction *******/
/***** from interaction energy calcs ****************************************************/
/***** Since we only care about the basic molecule only molecule A is assigned an *******/
/***** energy                                                                     *******/
/***** Started 25 July 1996 Dave Willock ************************************************/
/****************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "maxima.h"
#include "structures.h"
#include "ewald.h"
#include "own_maths.h"
#include "global_values.h"

void isol_coul(atom *p_mol_A, int num_mol_A_atoms, 
               atom *p_mol_B, int num_mol_B_atoms, int zero_first, int self)
{
#include "header.h"
int mol_index_A, mol_index_B;
int next_neigh, kneigh;

atom *p_atom_A, *p_atom_B;
double dx, dy, dz, separation;
double scale;

/***DEBUG DEBUG ****/
int jneigh, ineigh, is_neigh;
int neigh_index;
atom *p_neigh;
atom *p_next_neigh;
/***DEBUG DEBUG ****/

if (zero_first) 
  {
    p_atom_A= p_mol_A;
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
         p_atom_A->electrostatic_pot = 0;
         p_atom_A++;
      }
  }

/***************************************************************************************/
/**** Case when mol_A not mol_B ********************************************************/
/***************************************************************************************/

if (!self && !pbc)
  {

    p_atom_A= p_mol_A;
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
        p_atom_B= p_mol_B;
        for (mol_index_B=0; mol_index_B <= num_mol_B_atoms; mol_index_B++)
          {

            dx= (p_atom_A->x)-(p_atom_B->x);
            dy= (p_atom_A->y)-(p_atom_B->y);
            dz= (p_atom_A->z)-(p_atom_B->z);
 
            separation = sqrt (dx*dx + dy*dy + dz*dz);

            p_atom_A->electrostatic_pot += (p_atom_B->part_chge)/separation;

            p_atom_B++;
          }
        p_atom_A++;
      }
  }

/*****************************************************/
/*** proper internal electrostatic ignore atoms ******/
/*** in the same angle potential                ******/
/*** with 0.5 scaling for atoms in the same     ******/
/*** torsion potential DJW Feb 99               ******/
/*****************************************************/
else if (self && !pbc)
  {
    p_atom_A= p_mol_A;
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
        p_atom_B= p_mol_B+mol_index_A+1;
        for (mol_index_B=mol_index_A+1; mol_index_B <= num_mol_B_atoms; mol_index_B++)
          {
/**** ingnore nearest neighs ****/
             is_neigh=FALSE;
             scale=1.0;

             for (ineigh=0; ineigh <= p_atom_B->num_neigh; ineigh++) 
               {
                 neigh_index= p_atom_B->neighb[ineigh]; 
                 if (neigh_index == mol_index_A) is_neigh=TRUE;

/**** and next nearest neighs ***/
                  p_neigh= p_mol_B+neigh_index; 
                  for (jneigh=0; jneigh <= p_neigh->num_neigh; jneigh++) 
                   {
                     next_neigh= p_neigh->neighb[jneigh];
                     if (next_neigh == mol_index_A) 
                       {
                          is_neigh=TRUE;
                       }
                     else
                       {
/**** scale interaction over torsion ends ******/
                          p_next_neigh= p_mol_B+next_neigh;
                          for (kneigh=0; kneigh <= p_neigh->num_neigh; kneigh++) 
                            {
                               if (p_next_neigh->neighb[kneigh] == mol_index_A ) scale = 0.5;
                            }
                       }
                   }  
               }              

             if (!is_neigh)
               {

                  printf("In isol_coul doing %s (%10.6f)  with %s (%10.6f)",
                           p_atom_A->label, p_atom_A->part_chge, 
                           p_atom_B->label, p_atom_B->part_chge);

                  dx= (p_atom_A->x)-(p_atom_B->x);
                  dy= (p_atom_A->y)-(p_atom_B->y);
                  dz= (p_atom_A->z)-(p_atom_B->z);

                  separation = sqrt (dx*dx + dy*dy + dz*dz);

                  p_atom_A->electrostatic_pot += scale*(p_atom_B->part_chge)/separation;
               }

             p_atom_B++;
          }
        p_atom_A++;
      }
  }

/***************************************************************************************/
/**** Case when we wish to remove mol_A/mol_A interactions from the Ewald sum **********/
/***************************************************************************************/
else
  {
    p_atom_A= p_mol_A;
    for (mol_index_A=0; mol_index_A <= num_mol_A_atoms; mol_index_A++)
      {
        p_atom_B= p_mol_B+mol_index_A+1;
        for (mol_index_B=mol_index_A+1; mol_index_B <= num_mol_B_atoms; mol_index_B++)
          {

             dx= (p_atom_A->x)-(p_atom_B->x);
             dy= (p_atom_A->y)-(p_atom_B->y);
             dz= (p_atom_A->z)-(p_atom_B->z);
  
             separation = sqrt (dx*dx + dy*dy + dz*dz);

             p_atom_A->electrostatic_pot -= (p_atom_B->part_chge)/separation;

             p_atom_B++;
          }
        p_atom_A++;
      }
  }

return;
}

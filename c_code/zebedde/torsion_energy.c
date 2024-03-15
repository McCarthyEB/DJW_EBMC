/**************************************************************************/
/******* Intra-molecular angle torsion calculation based on       *********/
/******* pre-assembled list.                                      *********/
/******* Dave Willock Aug 2006                                    *********/
/*******                                                          *********/
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"
#include "constants.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double torsion_energy(atom *p_molecule, int num_atoms,
                      torsion_interact_list *p_torsions_list, int num_torsions_listed)
{
#include "header.h"
atom *p_atom1, *p_atom2, *p_atom3, *p_atom4;

double vec1[3], vec2[3], norm1[3], norm2[3], phi, size1, size2;
double dot;

int itorsion, j, pot_index;

double energy;

/*******DEBUG*********/
/* DEBUG=TRUE; */
/*******DEBUG set false at end again*********/

energy = 0.0;

//if (DEBUG) printf("DEBUG >> doing torsions\n");

for (itorsion=0; itorsion<=num_torsions_listed ;itorsion++)
   {
      pot_index = p_torsions_list->itorpot;

/******************************************************/
/** set atom pointers so torsion 1-2-3 is what we want **/
/******************************************************/
//      if (DEBUG) printf("Torsion listed as: %d %d %d %d\n",
//                              p_torsions_list->iatm1,
//                              p_torsions_list->iatm2,
//                              p_torsions_list->iatm3,
//                              p_torsions_list->iatm4);


      p_atom1 = p_molecule +  p_torsions_list->iatm1;
      p_atom2 = p_molecule +  p_torsions_list->iatm2;
      p_atom3 = p_molecule +  p_torsions_list->iatm3;
      p_atom4 = p_molecule +  p_torsions_list->iatm4;

      p_torsions_list++;

//       if (DEBUG) printf("Thats atoms %s >>%s<< %s >>%s<< %s >>%s<< %s >>%s<<\n",
//                                  p_atom1->label, p_atom1->pot,
//                                  p_atom2->label, p_atom2->pot,
//                                  p_atom3->label, p_atom3->pot,
//                                  p_atom4->label, p_atom4->pot);

/******************************************************/
/*** Work out torsion angle  1-2-3-4 ******************/
/******************************************************/

      vec1[0] = p_atom1->x - p_atom2->x;
      vec1[1] = p_atom1->y - p_atom2->y;
      vec1[2] = p_atom1->z - p_atom2->z;

      vec2[0] = p_atom3->x - p_atom2->x;
      vec2[1] = p_atom3->y - p_atom2->y;
      vec2[2] = p_atom3->z - p_atom2->z;

/***** vector cross product of 2->1 and 2->3 *********/ 

      vec_cross(&vec1[0], &vec2[0], &norm1[0]);

      vec1[0] = p_atom4->x - p_atom3->x;
      vec1[1] = p_atom4->y - p_atom3->y;
      vec1[2] = p_atom4->z - p_atom3->z;

/*** note order of vectors                   *********/
/***** vector cross product of 3->4 and 2->3 *********/ 
/***** equivalent to cross  of 3->2 and 3->4 *********/ 
      
      vec_cross(&vec1[0], &vec2[0], &norm2[0]);

      dot = norm1[0]*norm2[0]+norm1[1]*norm2[1]+norm1[2]*norm2[2];

      size1 = sqrt(norm1[0]*norm1[0]+norm1[1]*norm1[1]+norm1[2]*norm1[2]);
      size2 = sqrt(norm2[0]*norm2[0]+norm2[1]*norm2[1]+norm2[2]*norm2[2]);

//      if (DEBUG)
//        {
//          printf("vectors for torsion: norm1 %10.6f %10.6f %10.6f, size %10.6f\n",
//                                    norm1[0], norm1[1], norm1[2], size1);
//          printf("                     norm2 %10.6f %10.6f %10.6f, size %10.6f\n",
//                                    norm2[0], norm2[1], norm2[2], size2);
//        }


      dot = dot/(size1*size2);

      if ( dot > 1.0 && dot < 1.0001) dot = 1.0;
      if ( dot < -1.0 && dot > -1.0001) dot = -1.0;

      if (dot < -1.0 || dot > 1.0)
        {
          printf("ERROR: Unreasonable torsion sent to torsion_energy routine.");
          printf("       dot product of normals is %10.6f not in range for cosine",
                                                                             dot);
          printf("       function.");
          exit(0);
        }

      phi = acos(dot);

//      if (DEBUG) printf("measured torsion %s %s %s %s is %10.6f degrees\n",
//                     p_atom1->label, p_atom2->label,
//                     p_atom3->label, p_atom4->label,
//                     phi*RAD_TO_DEG);

/*******************************************************/
/*** Use appropriate energy function *******************/
/*******************************************************/

      if (intra_torsion_potent[pot_index].which==TORSION_1)
        {
//           if (DEBUG) printf("Using torsion_1\n");
/***************************************************************/
/******* In Torsion_1 parameters stored as                     */
/******* KPhi       n          Phi0                            */
/*******   A        B           C      D (0)   E (0)   F (0)   */
/***************************************************************/

           energy += intra_torsion_potent[pot_index].A * 
                     ( 1.0 + cos(   intra_torsion_potent[pot_index].B * phi
                                  - intra_torsion_potent[pot_index].C / RAD_TO_DEG));

        }
      else if (intra_torsion_potent[pot_index].which==TORSION_3)
        {
//            if (DEBUG) 
//              {
//                printf("Using torsion_3\n");
//             
//                printf("with V(1)= %10.6f Phi(1)= %10.6f,  V(2)= %10.6f Phi(2)= %10.6f, V(3)= %10.6f Phi(3)= %10.6f\n",
//                                       intra_torsion_potent[pot_index].A,
//                                       intra_torsion_potent[pot_index].B,
//                                       intra_torsion_potent[pot_index].C,
//                                       intra_torsion_potent[pot_index].D,
//                                       intra_torsion_potent[pot_index].E,
//                                       intra_torsion_potent[pot_index].F);
//              }
/***************************************************************/
/******* In Torsion_3 parameters stored as                     */
/******* V(1)  Phi0(1) V(2) Phi0(2) V(3) Phi0(3)               */
/*******  A      B      C     D      E     F                   */
/***************************************************************/

/************************************************************/
/*** Note use of minus sign in torsion potential        *****/
/*** this gives correct potential shape for c-c-c-c        **/
/*** staggered vs eclipsed conformations.                  **/
/***                                                       **/
/*** original definition in new_pcff.frc file is           **/
/*** in error having the text book standard:               **/
/***  E = SUM(n=1,3) { V(n) * [ 1 + cos(n*Phi - Phi0(n)) ] **/
/***  rather than the implemented:                         **/
/***  E = SUM(n=1,3) { V(n) * [ 1 - cos(n*Phi - Phi0(n)) ] **/
/************************************************************/

           energy += intra_torsion_potent[pot_index].A * 
                     ( 1.0 - cos(  phi
                                  - intra_torsion_potent[pot_index].B / RAD_TO_DEG))
                    +intra_torsion_potent[pot_index].C * 
                     ( 1.0 - cos(  2.0*phi
                                - intra_torsion_potent[pot_index].D / RAD_TO_DEG))
                  +intra_torsion_potent[pot_index].E * 
                     ( 1.0 - cos(  3.0*phi
                                - intra_torsion_potent[pot_index].F / RAD_TO_DEG));

        }
      
   }
/*******DEBUG*********/
/* DEBUG=FALSE; */
/*******DEBUG set false at end again*********/
  return energy;
}

/**********************************************************************/
/* fix_potentials_cff91.c: for Discover run, reset potentials at      */
/*                       : atom                                       */
/* started Dewi 16/6/95 DWL/DJW 28/2/96                               */
/* Modified Alan May 09 to make it more generic. Counts number of     */
/* heavies and assigns a potential according to list in frc file      */
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int fix_potentials_cff91(atom *p_molecule, int this_atom)
{
#include "header.h"
int fixed=FALSE ;
int dummy;
int iatom;
int heavy;

atom *p_atom;


p_atom = p_molecule + this_atom;
if (DEBUG) printf("DB>> Num_neigh %i\n",p_atom->num_neigh);
if (DEBUG) printf("DB>> Label %s\n",p_atom->label);

/*** These types only valid for sp3 Carbons so only do this         ***/
/*** if there are 4 neighbours in total. added July 09 Dave Willock ***/

if (strcmp(p_atom->elem,"C") == 0 && p_atom->num_neigh == 3)    /**sp3 Carbon**/
  {

/*** Count the number of heavy neighbours for this carbon           ***/
/*** Many options for a hydrogen element type, Dave Willock July 09 ***/

     heavy=0;
     for (dummy=0; dummy <= p_atom->num_neigh; dummy++)
       {
          iatom = (p_atom->neighb[dummy]);
          if (DEBUG) printf("DB>> Elem %s\n",(p_molecule + iatom)->elem);

          if (  strcmp((p_molecule + iatom)->elem,"H") != 0 
             && strcmp((p_molecule + iatom)->elem,"HA") != 0 
             && strcmp((p_molecule + iatom)->elem,"HB") != 0 
             && strcmp((p_molecule + iatom)->elem,"D") != 0  ) heavy++;
       }
      if (DEBUG) printf("DB>> number heavies %i\n",heavy);

		 if (heavy == 0)
		 		 {
		 		 printf("WARNING: Atom %d labelled %s elem %s potential %s appears to be isolated carbon atom with only %d H neighbours\n",
                                                              this_atom, p_atom->label, p_atom->elem, p_atom->pot, (p_atom->num_neigh)+1);
                                 printf("Neighbours:\n");
                                 for (dummy=0; dummy <= p_atom->num_neigh; dummy++)
                                   {
                                      iatom = (p_atom->neighb[dummy]);
                                      printf("%d index %d label %s Elem %s\n",dummy, iatom, 
                                               (p_molecule + iatom)->label, (p_molecule + iatom)->elem);
                                   }
		 		 strcpy(p_atom->pot, "c");
		 		 }

		 if (heavy == 1)
		 		 {
		 		 if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);

		 		 strcpy(p_atom->pot, "c3");
		 		 if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
		 		 fixed = TRUE;
		 		 }		 

		 if (heavy == 2)
		 		 {
		         if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
		         strcpy(p_atom->pot, "c2");
		 		 if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
		         fixed = TRUE;
		 		 }

		 if (heavy == 3)
        		 {
		         if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
        		 strcpy(p_atom->pot, "c1");
		 		 if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
        		 fixed = TRUE;
		         }
		 }

else if (strcmp(p_atom->elem,"N") == 0)    /**Nitrogen**/
        {
        heavy=0;
        for (dummy=0; dummy <= p_atom->num_neigh; dummy++)
                {
                iatom = (p_atom->neighb[dummy]);
                if (DEBUG) printf("DB>> Elem %s\n",(p_molecule + iatom)->elem);
                        if (strcmp((p_molecule + iatom)->elem,"H") != 0)
                                {
                                heavy++;
                                }
                if (DEBUG) printf("DB>> number heavies %i\n",heavy); /*Dont really need this for N*/
                }
		 
		 if (p_atom->num_neigh == 2) /*Is an amine*/
		 		 {
		        if (DEBUG)  printf("DB>> Atom type is %s\n",p_atom->pot);
		 		 strcpy(p_atom->pot, "npc");
                if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
		 		 fixed = TRUE;		 
		 		 }
		 if (p_atom->num_neigh == 3) /*Is an ammonium*/
                {
                if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
                strcpy(p_atom->pot, "n4");
                if (DEBUG) printf("DB>> Atom type is %s\n",p_atom->pot);
                fixed = TRUE;   
                }
		 }

else if (strcmp(p_atom->elem, "O") == 0) /**Oxygen**/
		 {
		 if (strcmp(p_atom->pot, "oh") == 0)
		 		 {
                /**was bonded to H **/
		 		 strcpy(p_atom->pot, "o");
		 		 fixed = TRUE;
		 		 }
		 }

else if (strcmp(p_atom->elem, "S") == 0) /**sulphur*/
		 {
		 if (strcmp(p_atom->pot, "sh") == 0)
		 		 {
                /**was bonded to H **/
		 		 strcpy(p_atom->pot, "s");
		 		 fixed = TRUE;
		 		 }
		 }
return(fixed);
}


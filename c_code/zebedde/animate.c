/*****************************************************************************/
/* animate.c : writes .arc file of template during growth                    */
/*             and writes a insight .log to read and play                    */
/*             (Cos Insight is cack and can't read a growing molecule)       */
/*                                                                           */
/* STarted 28/3/95 Dewi                                                      */
/* Need to pass separate frame counters according to which animation file!   */
/* 6/6/96 DWL Mods to handle both Biosym and Xmol animations                 */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void print_molecule(atom *p_molecule, int num_atoms, FILE *output_fp, int D_to_H);

void print_molecule_xmol(atom *p_molecule, int num_atoms, FILE *output_fp);

void print_pbc_header(FILE *file_fp, double *p_abc);

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec, int which_mol);

void min_image( double *x, double *y, double *z);

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z, 
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void animate(atom *guest_ptrs[], int num_guests,
             list_partition *p_guest_demarc, FILE *an_fp, FILE *lr_fp, 
	     FILE *ls_fp, char *p_info, char *p_animation_file,
             int *p_frame_counter, atom *pore)
{
#include "header.h"

double cofm[3], vec[3], cell_centre[3], total_mass;
int imol, start_mol, index, tot_atoms;
int isymm, num_unique_atoms;

list_partition *p_demarc;

atom *p_image;

/* increment counter */
(*p_frame_counter)++;

if (DEBUG) printf("DB>>> in animate\n");

if (animate_flag == BIOSYM_ANIMATION)
	{

     if (DEBUG) printf("DB>> BIOSYM animation num_guests=%d, num_symm_ops=%d\n", num_guests, num_symm_ops);

/**** Move molecule by lattice translations so that its centre of mass is in the unit cell ***/
/**** Added July 09 Dave Willock.                                                          ***/
/**** Note that if no symmetry operations are set num_symm_ops = -1                        ***/

//        p_demarc= p_guest_demarc;
//        for (imol=0; imol< num_guests; imol++)
//          {
//            for (isymm=0; isymm <= num_symm_ops+1; isymm++)
//              {
//                index=imol*(num_symm_ops+2) + isymm;
//
//                printf("Index=%d\n", index);
//
//                centre_of_mass(&cofm[0], &total_mass, guest_ptrs[index], 
//                               p_demarc->end, -1 );
//
//                printf("cofm %10.6f %10.6f %10.6f\n", cofm[0], cofm[1], cofm[2] );
/*** shift molecule to origin ***/
//                vec[0] = -cofm[0];
//                vec[1] = -cofm[1];
//                vec[2] = -cofm[2];
//                move_molecule(guest_ptrs[index], p_demarc->end, &vec[0], -1);
//
//                printf("molecule moved\n", index);
//
//                min_image( &cofm[0], &cofm[1], &cofm[2]);
//
//                printf("min_image\n", index);
/*** shift molecule to centre of mass in unit cell ****/
//                move_molecule(guest_ptrs[index], p_demarc->end, &cofm[0], -1);
//                printf("molecule moved again\n", index);
//             }
//            p_demarc++;
//          }

      if (DEBUG) printf("DB>>> in animate BIOSYM_ANIMATION\n");
	/* write title and date */
	print_frame_header(p_info ,*p_frame_counter, an_fp);

/* Also print PBC data (update October 2006, Dave Willock) */
        print_pbc_header( an_fp, &abc[0]);

	/* append the molecules to the archive */
        p_demarc= p_guest_demarc;
        for (imol=0; imol<num_guests; imol++)
          {
            for (isymm=0; isymm <= num_symm_ops+1; isymm++)
              {
                index=imol*(num_symm_ops+2) + isymm;

	        print_molecule( guest_ptrs[index], p_demarc->end, an_fp, FALSE);
        	fprintf(an_fp,"end\n");
              }
            p_demarc++;
          }
	fprintf(an_fp,"end\n");

	/* now write out the log to read and display this frame */

        if (logfile_needed)
          {
              fprintf(lr_fp,"Get Molecule Archive Frame %d",*p_frame_counter);
              fprintf(lr_fp," %s anim%i Reference_Object pore\n", 
				p_animation_file, *p_frame_counter);

              fprintf(ls_fp,"display on anim%i\n",*p_frame_counter);
	      fprintf(ls_fp,"display off anim%i\n",*p_frame_counter);
          }
	}

else if (animate_flag == XMOL_ANIMATION_WITHPORE )
	{
           tot_atoms=0;
           p_demarc= p_guest_demarc;
           for (imol=0; imol< num_guests; imol++)
             {
               for (isymm=0; isymm <= num_symm_ops+1; isymm++) tot_atoms+=p_demarc->num;
               p_demarc++;
             }

	fprintf(an_fp,"%i\n", num_pore_atoms+1 + tot_atoms);
	fprintf(an_fp,"%s\n", p_info);

	print_molecule_xmol(&pore[0], num_pore_atoms, an_fp);

        p_demarc= p_guest_demarc;
        for (imol=0; imol< num_guests; imol++)
          {
            for (isymm=0; isymm <= num_symm_ops+1; isymm++)
              {
                index=imol*(num_symm_ops+2) + isymm;
	        print_molecule_xmol(guest_ptrs[index], p_demarc->end, an_fp);
              }
            p_demarc++;
          }
	}

else if (animate_flag == XMOL_ANIMATION_NOPORE )
	{
           tot_atoms=0;
           p_demarc= p_guest_demarc;
           for (imol=0; imol< num_guests; imol++)
             {
               for (isymm=0; isymm <= num_symm_ops+1; isymm++) tot_atoms+=p_demarc->num;
               p_demarc++;
             }
	   fprintf(an_fp,"%i\n", tot_atoms); 
   	   fprintf(an_fp,"%s\n", p_info);
	
           p_demarc= p_guest_demarc;
           for (imol=0; imol< num_guests; imol++)
             {
               for (isymm=0; isymm <= num_symm_ops+1; isymm++)
                 {
                   index=imol*(num_symm_ops+2) + isymm;
                   print_molecule_xmol(guest_ptrs[index], p_demarc->end, an_fp);
                 }
              p_demarc++;
            }
	}
 else if (animate_flag == FALSE)
       {
       printf("No animation requested\nReturning\n");
       }
else
    {
	printf("SERIOUS ERROR: Unknown Animation type %d\n", animate_flag);
        fprintf(output_fp,"SERIOUS ZEBEDDE ERROR: Exiting\n");
        fflush(output_fp);
        exit(EXIT_FAILURE);
    }

return;
}


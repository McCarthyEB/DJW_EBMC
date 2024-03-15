/*****************************************************************************/
/*   initialise_animation.c                                                  */
/*   16.10.95 started : initialises animation file N and writes headers      */
/*   29/2/96          : modify to handle a seed number as well               */
/*   14/11/06         : modify to stop logfiles being openned if not required*/
/*   NB!!!! -1 as a use_number is the flag for initialise_gooduns       */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"


void print_biosym_header(FILE *file_fp, int pbc_flag);
void print_log_header(char *p_title, FILE *file_fp);
int open_file(FILE **p_file, char *p_filename, char *p_status);
char int_to_char (int num);

void initialise_animation(char *p_this_animation_file, char *p_animation_file,
                          int seed_type, int seed_number, int use_number, 
                          int build_number)
{

#include "header.h"
char f_title[80];
char animation_readlog[FILELEN_MAX];
char animation_showlog[FILELEN_MAX];
int initialise_gooduns;

if (DEBUG)
  {
    printf("DB>>>> In initialise_animation: seed_type %d, seed_number %d, use_number %d\n", 
                                               seed_type,seed_number,use_number);
    printf("DB>>>> In initialise_animation: animation file stem >>%s<<\n", p_animation_file);
  }
/*****************************************************************************/
/****** work out if this is a call to do gooduns or just the latest arc file */
/*****************************************************************************/
    initialise_gooduns = FALSE;


    if (use_number < -1 )
      {
        fprintf(stderr, 
			"Negative use_number in initialise animation: aborting!!!\n");
        exit(1);
      }
    else if (use_number < 0)
      {
      initialise_gooduns = TRUE;
      }
    else if (use_number >= 1000)
      {
        fprintf(stderr, "Error: Too many templates requested cannot cope in initialise_animation\n");
        fprintf(stderr, "Please restrict yourself to 999 templates or less! EXITING\n");
        exit(1);
      }

    if (initialise_gooduns)
      {
        strcpy(p_this_animation_file, p_animation_file);
      }
    else
      {
        strcpy(p_this_animation_file, p_animation_file);
        if (seed_type != MOLE)
          {
             sprintf(p_this_animation_file,"%s%d_%d_%d",p_animation_file,
                                       use_number, seed_number, build_number);
          }
        else
          {
             sprintf(p_this_animation_file,"%s%d_%d",p_animation_file,
                                                seed_number, build_number);
          }
      }


/*****************************************************************************/
/****** Make up a new .arc file by attaching the seed and template number ****/
/*****************************************************************************/

    strcpy(animation_readlog, p_this_animation_file);
    strcpy(animation_showlog, p_this_animation_file);
    strcat(animation_readlog, "_read.log");
    strcat(animation_showlog, "_show.log");
    strcat(p_this_animation_file, ".arc");

/*****************************************************************************/
/****** Assign the file pointers *********************************************/
/*****************************************************************************/
    if (!initialise_gooduns)
		{
        if (open_file(&anim_file_fp, p_this_animation_file,"w") == EXIT_FAILURE)
          {
             printf("ERROR>> Cannot open file %s for writing animations.", p_this_animation_file);
             exit(0);
          }

        if (logfile_needed)
          {
            if (open_file(&anim_read_fp, animation_readlog, "w") == EXIT_FAILURE)
              {
                 printf("ERROR>> Cannot open file %s for read log", animation_readlog);
                 exit(0);
              }

            if (open_file(&anim_show_fp, animation_showlog,"w") == EXIT_FAILURE)
              {
                 printf("ERROR>> Cannot open file %s for show log", animation_showlog);
                 exit(0);
              }
          }

		/********print header*********************************/
	    print_biosym_header(anim_file_fp, pbc);

		/********print header in log files if required! ********************/
                if (logfile_needed) 
                  { 
		    sprintf(f_title,"#To read in animation %s",p_this_animation_file);
                    print_log_header(f_title, anim_read_fp);
                    fprintf(anim_read_fp,"Get Molecule Archive Frame 1 %s pore -Reference_Object\n",
					pore_file);
                    sprintf(f_title,"#To show the  animation %s",p_animation_file);
	   	    print_log_header(f_title, anim_show_fp);
		    fprintf(anim_show_fp,"display off anim*\n");
                 }
	
	    }	
    else
		{

		if (open_file(&gooduns_fp, p_this_animation_file,"w") == EXIT_FAILURE)
                  {
                     printf("ERROR>> Cannot open file %s for gooduns", p_this_animation_file);
                     exit(0);
                  }

		/********print header ********************************/
                print_biosym_header(gooduns_fp, pbc);

                if (logfile_needed)
                  {
                    if (open_file(&gooduns_read_fp, animation_readlog, "w")== EXIT_FAILURE)
                      {
                         printf("ERROR>> Cannot open file %s for gooduns read log", animation_readlog);
                         exit(0);
                      }

                    open_file(&gooduns_show_fp, animation_showlog,"w");
                      {
                         printf("ERROR>> Cannot open file %s for gooduns show log", animation_showlog);
                         exit(0);
                      }

		/********print header in log files********************/
	    	    sprintf(f_title,"#To read in animation %s",p_this_animation_file);
		    print_log_header(f_title, gooduns_read_fp);
                    fprintf(gooduns_read_fp,"Get Molecule Archive Frame 1 %s pore -Reference_Object\n",
					pore_file);
	    	    sprintf(f_title,"#To show the  animation %s",p_animation_file);
		    print_log_header(f_title, gooduns_show_fp);
		    fprintf(gooduns_show_fp,"display off anim*\n");
                  }


		}	
	return;
}


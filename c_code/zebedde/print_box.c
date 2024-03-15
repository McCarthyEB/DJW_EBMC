/*******************************************************/
/***** print_box.c :     Writes out a cat/mdf    *******/
/*****                   for the box limits used *******/
/***** Dewi 17/5/95                              *******/
/*******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void box_n_print(FILE *b_fp ,int a, int b,int c);

void print_biosym_header(FILE *file_fp, int pbc_flag);

void car_end(FILE *fp);

void print_mdf_header(char *p_assembly, FILE *fp);

void print_frame_header(char *p_title, int num_frame, FILE *file_fp);

void print_mdf_end(FILE *fp);

/*** NOTE p_bl short for p_box_limits for convenience*****/

void print_box(double *p_bl, char *box_file)
{
#include "header.h"
int i,x,y,z;
FILE *box_fp;
char b_title[80], b_file[FILELEN_MAX];
char dummy[10];

printf("print_box\n");
for (i=0;i<=5;i++) printf("%f \n",p_bl[i]);

/******************write out the .car file***************************/
printf("print_box.car\n");

strcpy(b_file, box_file);
strcat(b_file,".car");

if (!(box_fp = fopen(b_file,"w")))
	{
	printf("ERROR : Could not open car file %s in print_box\n", box_file);
	exit(1);
	}

print_biosym_header(box_fp, 0);

sprintf(b_title,"BOX for run %s",outputfile);
print_frame_header(b_title, 1, box_fp);

i=1;               /* label counter */
for (x=0;x<=1;x++)
	{
	for (y=2;y<=3;y++)
   	 	{
		for (z=4;z<=5;z++)
   	 		{
			sprintf(dummy,"BOX%i",i++);
			fprintf(box_fp,"%-5s %14.9f %14.9f %14.9f %-4s 1      %-3s     %-2s %6.3f\n", 
                                                    dummy, p_bl[x],p_bl[y],p_bl[z], "GRP","h", "D",0.0);
			}
		}
	}

car_end(box_fp);

/******************write out the .mdf file***************************/
printf("print_box.mdf\n");

strcpy(b_file, box_file);
strcat(b_file,".mdf");

if (!(box_fp = fopen(b_file,"w")))
        {
        printf("ERROR : Could not open mdf file %s in print_box\n", box_file);
        exit(1);
        }


print_mdf_header(NULL, box_fp);
fprintf(box_fp, "@molecule box\n\n");

i=1; /*label counter */
for (x=0;x<=1;x++)
        {
        for (y=0;y<=1;y++)
                {
                for (z=0;z<=1;z++)
                        {
						sprintf(dummy,"GRP_1:BOX%i",i);
						fprintf(box_fp,"%-20s",dummy);
						fprintf(box_fp, "%-2s %-3s","D","h");
						fprintf(box_fp, "     %-4s","GRP");
						fprintf(box_fp, " 0  0");
  						fprintf(box_fp, "  %7.4f",0.0);
						fprintf(box_fp, " 0 0 8");
						fprintf(box_fp, "  1.0000  0.0000");

						/* print neighbours */
						switch (i)
							{
							case 1: box_n_print(box_fp,2,3,5);break;
							case 2: box_n_print(box_fp,1,4,6);break;
							case 3: box_n_print(box_fp,1,4,7);break;
							case 4: box_n_print(box_fp,2,3,8);break;
							case 5: box_n_print(box_fp,1,6,7);break;
							case 6: box_n_print(box_fp,2,5,8);break;
							case 7: box_n_print(box_fp,3,5,8);break;
							case 8: box_n_print(box_fp,4,6,7);break;
							}
						i++;
						}
				}
		}
fprintf(box_fp,"\n\n");

print_mdf_end(box_fp);

return;
}

void box_n_print(FILE *b_fp ,int a, int b,int c)
{
fprintf(b_fp, " BOX%i BOX%i BOX%i\n",a,b,c);

return;
}

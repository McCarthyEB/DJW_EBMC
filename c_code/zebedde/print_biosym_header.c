/************************************************************/
/*      print_biosym_header.c : utilities for               */
/*      printing the tops of diff                           */
/*      Biosym files                                        */
/* dewi 28/3/95	                                            */
/* Adapted for Materials Studio to always print Date line   */
/* dave Oct. 2006                                           */
/************************************************************/
#include <stdio.h>
#include <time.h>

/****** prints the header bit in an .arc file *****/
void print_biosym_header(FILE *file_fp, int pbc_flag)
{
fprintf(file_fp,"!BIOSYM archive 3\n");
if (!pbc_flag) fprintf(file_fp,"PBC=OFF\n");
else  fprintf(file_fp,"PBC=ON\n");

return;
}

/****** prints the frame header in an .arc file *****/
void print_frame_header(char *p_title, int num_frame, FILE *file_fp)
{
time_t idate;
idate = time(NULL);
if (num_frame<0) /***ie. not an animation ***/
	{
	fprintf(file_fp,"%s\n",p_title);
	}
else
	{
	 fprintf(file_fp,"%s frame number %d\n",p_title,  num_frame);
	}
/* fprintf(file_fp,"!DATE %24.24s\n",ctime(&idate)); */
 fprintf(file_fp,"!DATE Mon Oct 12 12:17:26 1998\n"); 

return;
}

void print_pbc_header(FILE *file_fp, double *p_abc)
{

fprintf(file_fp,"PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f (P1)\n",
		*p_abc, *(p_abc+1), *(p_abc+2), *(p_abc+3), *(p_abc+4), *(p_abc+5));

return;
}

/****** prints the header for a .log file       *****/
void print_log_header(char *p_title, FILE *file_fp)
{

time_t idate;
idate = time(NULL);
fprintf(file_fp,"#Running: insightII\n");
fprintf(file_fp,"#Log Created: %24.24s\n",ctime(&idate));
fprintf(file_fp,"%s\n",p_title);
return;
}

void car_end(FILE *fp)
{
fprintf(fp,"end\n");
fprintf(fp,"end\n");
fclose(fp);
}

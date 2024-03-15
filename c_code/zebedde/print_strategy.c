/************************************************************/
/*      print_strategy.c : routines for writing strategy    */
/*      files for DISCOVER                                  */
/*              		                          	*/
/* dewi 3/5/95												*/
/************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "maxima.h"
#include "structures.h"
#include "header.h"

#define SWTDIS 1.5 /* switching distance (a sensible number according to manual) */
#define GOOD_CTF 10.0 /* Use this cutoff is user is daft and gives nb_ctf < SWTDIS */

void print_strategy_header(FILE *file_fp);
void print_strategy_minimize(FILE *file_fp);
void print_strategy_fixhost(FILE *file_fp);

void print_template_strategy(char *file)
{

FILE *file_fp;

printf("Writing template strategy file %s (None supplied)\n", file);

if (!(file_fp = fopen(file,"w")))
    {
    printf("error opening Discover template strategy file %s.inp\n", file);
    printf("Errno = %d\n",errno);
    printf("Serious ERROR\n");
    exit(-1);
    }


print_strategy_header(file_fp);
print_strategy_minimize(file_fp);
fclose(file_fp);
return;
}

void print_inpore_strategy(char *file)
{

FILE *file_fp;

printf("Writing inpore strategy file %s (None supplied)\n", file);


if (!(file_fp = fopen(file,"w")))
    {
    printf("error opening Discover inpore strategy file %s.inp\n", file);
    printf("Errno = %d\n",errno);
    printf("Serious ERROR\n");
    exit(-1);
    }

print_strategy_header(file_fp);
print_strategy_fixhost(file_fp);
print_strategy_minimize(file_fp);
 
fclose(file_fp);
return;
}

/****** prints the top bit of a strategy file *******/
/****** including the cutoffs and other flags *******/

void print_strategy_header(FILE *file_fp)
{
fprintf(file_fp,"! Discover input file written from ZEBEDDE\n");
fprintf(file_fp,"\toverlap =      0.0\n");
fprintf(file_fp,"\tigrpck  =      0\n");
fprintf(file_fp,"\titrap   =      1\n");
if (nb_ctf > SWTDIS)
  {
    fprintf(file_fp,"\tcutoff              =      %10.6f\n", nb_ctf);
    fprintf(file_fp,"\tcutdis              =      %10.6f\n", nb_ctf - SWTDIS);
  } 
else
  {
    fprintf(output_fp, "\n\n WARNING: nb_ctf is less than a sensible value for DISCOVER");
    fprintf(output_fp, "\n WARNING: Setting to %10.4f for Minimisations\n", GOOD_CTF);
    fprintf(file_fp,"\tcutoff              =      %10.6f\n", GOOD_CTF);
    fprintf(file_fp,"\tcutdis              =      %10.6f\n", GOOD_CTF - SWTDIS);
  }
fprintf(file_fp,"\tswtdis              =      %10.6f\n", SWTDIS);
fprintf(file_fp,"\treduce\n");
fprintf(file_fp,"\tmipbc               =      0\n");
fprintf(file_fp,"      begin simulation\n");
fprintf(file_fp,"     *   add-automatic bond torsion valence out-of-plane\n");
fprintf(file_fp,"     reduce\n");
fprintf(file_fp,"     set dielectric  = 1.000000\n");


return;
}

void print_strategy_fixhost(FILE *file_fp)
{
/***** fixes MOLECULE ONE (the framework)****/
fprintf(file_fp,"      fixed atom list generation\n");
fprintf(file_fp,"     *    add all atoms\n");
fprintf(file_fp,"     *    molecule 1\n");

return;
}

void print_strategy_minimize(FILE *file_fp)
{
fprintf(file_fp,"      minimize\n");
fprintf(file_fp,"     *    no cross terms\n");
fprintf(file_fp,"     *    no morse\n");
fprintf(file_fp,"     *    and no charges\n");
fprintf(file_fp,"     *    for 250 iterations\n");
fprintf(file_fp,"     *    using steepest\n");
fprintf(file_fp,"     *    until the maximum derivative is less than 10.0000 kcal/A\n");
fprintf(file_fp,"      end\n");

return;
}

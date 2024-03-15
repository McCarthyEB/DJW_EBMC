/**********************************************************/
/* Banner.c                                               */
/**********************************************************/
#include <stdio.h>


void banner(FILE *fp)
{


fprintf(fp," ------------------------------------------------------------------------ \n");
fprintf(fp,"|                                                                        |\n");
fprintf(fp,"|                                                             TTTTT      |\n");
fprintf(fp,"|                                                               T        |\n");
fprintf(fp,"|      ZZZZ    EEEE    BBB     EEEE    DDD     DDD     EEEE     T        |\n");
fprintf(fp,"|        Z     E       B  B    E       D  D    D  D    E        T        |\n");
fprintf(fp,"|       Z      EEEE    BBBB    EEEE    D  D    D  D    EEEE              |\n");
fprintf(fp,"|      Z       E       B  B    E       D  D    D  D    E                 |\n");
fprintf(fp,"|      ZZZZ    EEEE    BBB     EEEE    DDD     DDD     EEEE              |\n");
fprintf(fp,"|                                                                        |\n");
fprintf(fp,"|           ZEolites By Evolutionary De-novo DEsign of Templates         |\n");
fprintf(fp,"|           ==       =  =            =       ==        =                 |\n");
fprintf(fp,"|                                                                        |\n");
fprintf(fp,"|      Version 3.1 alpha realease : last altered 30th July 2007          |\n");
fprintf(fp,"|      Modifications from Version 2.0 include:                           |\n");
fprintf(fp,"|      1) AMBER forcefield implemented in costfunction                   |\n");
fprintf(fp,"|      2) c2discover calling implemented                                 |\n");
fprintf(fp,"|      3) line constraints                                               |\n");
fprintf(fp,"|      4) Preset dihedrals for new bonds                                 |\n");
fprintf(fp,"|      5) vdw forcefield with Monte Carlo checking of actions            |\n");
fprintf(fp,"|                                                                        |\n");
fprintf(fp,"|      Authors: Dewi W. Lewis and Dave J. Willock                        |\n");
fprintf(fp,"|                                                                        |\n");
fprintf(fp," ------------------------------------------------------------------------ \n");

return;
}

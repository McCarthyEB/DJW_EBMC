/***************************************************/
/***** Random number generator to give a      ******/
/***** random double 0.0 -> 1.0 using the     ******/
/***** standard library routines random and   ******/
/***** srandom with optional seeding on the   ******/
/***** current calander time                  ******/
/*****                                        ******/
/***** Dave Willock May 95                    ******/
/*****                                        ******/
/***** Modified 26.10.96 DWL                  ******/
/***** To allow user specification of no      ******/
/***** random seeding via input file          ******/
/***** Previously if done = 0, then seed on   ******/
/***** calender time and if time_now=1 then   ******/
/***** seed at 1. Now allow -1 to set to seed ******/
/***** at 1.                                  ******/
/***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double real_random(int done)
{
long irand;
double real_rand;
time_t time_now;

/**** seed if done = 0 **************************************/
if (!done)
   {
      time_now= time(&time_now); 
      /* for debugging */
/*       time_now= -1;  */

      if (time_now == -1)
        {
          printf("DEBUGGING::::::: Warning: Calander time not available seeding with 1\n");
          time_now = 1;
        }

     printf("Random seed is: %d\n", (unsigned int) time_now);
         
     srand((unsigned int) time_now ); 
   }
else if (done == -1)
   {
    time_now = 1;
    srand((unsigned int) time_now );
   }

    irand = rand();
    real_rand = irand;

return real_rand/RAND_MAX;
}

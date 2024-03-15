/***************************************************************/
/* timer.c timer routine nicked from Marvin                    */
/***************************************************************/
#include <time.h>

double timer(time_t reference_t)
{

  time_t     tick;
  double     elapsed;

  tick = time(NULL);
  elapsed = difftime(tick, reference_t);
  if (elapsed == 0) elapsed = 1E-6;
  return(elapsed);
}


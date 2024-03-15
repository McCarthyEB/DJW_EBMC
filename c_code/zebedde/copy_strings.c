#include <stdio.h>

/* subroutine for copying one string to another */

void copy_strings( char to[], char from[] )

{

int iloop;

   iloop=0;
   while ( (to[iloop] = from[iloop]) != '\0')
      {
      iloop++;
      }

  return;
}





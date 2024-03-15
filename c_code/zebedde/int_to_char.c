/******************************************************/
/***** convert an integer to its string equivalent ****/
/***** restricted to single characters (0->9) for now */
/***** Dave Willock                               *****/
/******************************************************/

#include <stdio.h>
#include <stdlib.h>

char int_to_char( int num)
{
   
   if (num < 0 || num > 9)
     {
        printf("Out of range number sent to int_to_char: num= %d\n",num);
        exit(1);
     } 
   return  '0'+num;
}


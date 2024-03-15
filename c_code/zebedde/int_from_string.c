#include <stdio.h>

int int_from_string(char *p_string, int *p_int)
  {
   
int iloop;

    iloop = 0;
    while (*p_string != '\0')
       {
          *p_int = *p_string;
          p_int++;
          p_string++;
          iloop++;
       }
   return iloop;
  }

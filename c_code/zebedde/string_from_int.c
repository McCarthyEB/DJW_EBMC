#include <stdio.h>

int string_from_int(int *p_int, char *p_string)
  {
   
int iloop;

iloop=0;
    while (*p_int != '\0')
       {
          *p_string= *p_int;
          p_int++;
          p_string++;
          iloop++;
       }
   *p_string= *p_int;
   return iloop;
  }

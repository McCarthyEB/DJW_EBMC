#include <stdio.h>
#include <limits.h>
#include "global_values.h"

/******************************************************/
/***** p_key is searched for in p_char. key can be ****/
/***** of any length, p_char should be of          ****/
/***** num_of_chars length                         ****/
/***** Dave Willock November 1995                  ****/
/******************************************************/

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars )

{

int iloop;
char *p_start;

  p_start= p_key;

  for (iloop=0; iloop <= num_of_chars; ++iloop)
    {
      if ( *(p_char+iloop) == *p_key ) p_key++; 
                                     else p_key = p_start;
      if (*(p_key) == '\0') return TRUE;

/* Bail out if end of test string found **********************/

      if (*(p_char+iloop) == '\0') return FALSE;
    }

  return FALSE;
}


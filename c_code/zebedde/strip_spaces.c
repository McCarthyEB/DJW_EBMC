/***************************************************/
/****** Routine to remove trailing spaces from *****/
/****** a string                              ******/
/****** DWL 3/7/96                            ******/
/****** Simon pinched 4/7/96                  ******/
/***************************************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>

void strip_spaces(char *string)
{
  register int i = 0;

  while (isspace(string[i]))
    i++;

  if (i)
    strcpy(string, string + i);

  i = strlen(string) - 1;

  while ((i > 0) && isspace(string[i]))
    i--;

  string[++i] = '\0';

return;
}

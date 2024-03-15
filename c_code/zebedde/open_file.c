/******************************************************************************/
/* open_file.c opens file from file pointer                                   */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

int open_file(FILE **p_file, char *p_filename, char *p_status)
{
// printf("p_filename %s, p_status %s\n",p_filename,p_status);
if ((*p_file = fopen(p_filename,p_status)) == NULL)
          {
            fprintf(stderr, "Unable to open file %s.\tAborting\n",
                        p_filename);
            return EXIT_FAILURE;
          }
/* printf("p_filename %s, p_status %s\n",p_filename,p_status); */
return 0;
}

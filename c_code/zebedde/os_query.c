/******************************************************/
/******** Code to find out what OS we are using *******/
/******** AJWL 07/08                            *******/
/******** Achieved using the uname command      *******/
/******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void os_query()

{
#include "header.h"
FILE *uname_fp;
char command[256];
char output[256];
int size;
sprintf(command, "uname -o > uname_temp");
system(command);

uname_fp = fopen("uname_temp", "r");
fgets(buffer,BUFFER,uname_fp);
strcpy(output,buffer);
printf("%s\n",output);
size = strlen(output);

if (size ==0)	
	{
	printf("OS = Windows\n");
	os_windows = TRUE;
	os_linux = FALSE;
	}
else if ((strstr(output, "Linux")) !=NULL)
	{
	printf("OS = LINUX\n");
	os_windows = FALSE;
        os_linux = TRUE;
	}
else 
	{
	printf("OS unknown - assuming unix based OS\n");
	os_windows = FALSE;
	os_linux = TRUE;
	}

fclose(uname_fp);
if (os_linux)		
	{
	sprintf(command, "rm uname_temp");
	system(command);
	}
if (os_windows)
	{
	sprintf(command, "del uname_temp");
        system(command);
        }
return;
}

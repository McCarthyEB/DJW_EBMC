#include <stdio.h>
#include <errno.h>
void delete_file(char *fileroot, char *extension)
{

char dummy[256];
int error;

sprintf(dummy,"%s%s",fileroot,extension);
error = remove(dummy);

if (error == -1)
	{
	printf("Warning: Error %i whilst deleting file %s\n",
				errno, dummy);
	}
return;
}

/*******************************************************************************/
/* find_zebedde_file (UNIX version)                                            */
/*  Searches various directories trying to find the file file_name.  If the    */
/*  file is found its file pointer is returned.  The directories are searched  */
/*  in the order (i) the current directory (ii) the directory in shell variable*/
/*  ZEBEDDE_PATH (iii) the directory prepended to the command name in argv[0]  */
/*  (iv) the directories contained in the $PATH shell variable and (v) a hard  */
/*  coded directory in FINAL_PATH.                                             */
/*                                                                             */
/* Pinched from find_marvin_file (D Gay / A Rohl) 26.10.96 DWL                 */
/* char *executing_file is the full filename of the ZEBEDDE executable         */
/*                                                                             */
/*                                                                             */
/*******************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
#include "data.h"
#include "global_values.h"



char *find_zebedde_file(char *file_name, char *executing_file)
{
  FILE *test_file(char *path_name, char *file_name);

  char  *env;
  FILE  *fp;
  char   path_name[ZEB_PATH_MAX+1];
  char   path[ZEB_PATH_MAX+1];
  char  *path_list;
  char  *last_slash, *ind_path;
  
  if ((fp=test_file(NULL, file_name)) != NULL)
    {
    /*printf("Returning with %s\n", file_name);*/
    return(file_name);
    }
  
  env = getenv("ZEBEDDE_PATH");
  if (env != NULL)
    if ((fp = test_file(env, file_name)) != NULL)
      {
      strncpy(path, env, ZEB_PATH_MAX);
      strncat(path, "/", ZEB_PATH_MAX);
      strncat(path, file_name, ZEB_PATH_MAX);
      /*printf("Returning with %s\n", path);*/
      strcpy(file_name, path);
      return(file_name);
      }
  
  strncpy(path_name, executing_file, ZEB_PATH_MAX);
  last_slash = strrchr(path_name, '/');
  if (last_slash != NULL)
    *last_slash = '\0';
  else
    strncpy(path_name, ".", ZEB_PATH_MAX);
  if ((fp=test_file(path_name, file_name)) != NULL)
      {
      strncpy(path, path_name, ZEB_PATH_MAX);
      strncat(path, "/", ZEB_PATH_MAX);
      strncat(path, file_name, ZEB_PATH_MAX);
      /*printf("Returning with %s\n", path);*/
      strcpy(file_name, path);
      return(file_name);
      }
  
  env = getenv("PATH");
  /* path_list = (char *)our_alloc(NULL, (strlen(env)+1)*sizeof(char), "copying PATH environment variable"); */
  /** malloc space for path_list **/
  path_list = malloc(strlen(env)+1);
  strcpy(path_list, env);
  /*printf("PATH_LIST: %s\n",path_list);*/
  ind_path = strtok(path_list, ":");
  if ((fp=test_file(ind_path, file_name)) != NULL)
      {
      strncpy(path, ind_path, ZEB_PATH_MAX);
      strncat(path, "/", ZEB_PATH_MAX);
      strncat(path, file_name, ZEB_PATH_MAX);
      /*printf("Returning with %s\n", path);*/
      strcpy(file_name, path);
      return(file_name);
      }

  while ((ind_path = strtok(NULL, ":")) != NULL)
    if ((fp=test_file(ind_path, file_name)) != NULL)
      {
      strncpy(path, ind_path, ZEB_PATH_MAX);
      strncat(path, "/", ZEB_PATH_MAX);
      strncat(path, file_name, ZEB_PATH_MAX);
      /*printf("Returning with %s\n", path);*/
      strcpy(file_name, path);
      return(file_name);
      }

  
  free(path_list);
  
  if ((fp=test_file(FINAL_PATH, file_name)) != NULL)
      {
      strncpy(path, FINAL_PATH, ZEB_PATH_MAX);
      strncat(path, "/", ZEB_PATH_MAX);
      strncat(path, file_name, ZEB_PATH_MAX);
      /*printf("Returning with %s\n", path);*/
      strcpy(file_name, path);
      return(file_name);
      }

  
  return(NOT_SET); /* return for fatal condition */
}


/******************************************************************************
 * test_file
 *  returns the file pointer to the file file_name in directory path_name if
 *  it exists.  If it doesn't exist, NULL is returned.
 ******************************************************************************/

FILE *test_file(char *path_name, char *file_name)
{
  char   path[ZEB_PATH_MAX+1];
  FILE  *fp;
  
  if (path_name != NULL) {
    strncpy(path, path_name, ZEB_PATH_MAX);
    strncat(path, "/", ZEB_PATH_MAX);
    strncat(path, file_name, ZEB_PATH_MAX);
  }
  else
    strncpy(path, file_name, ZEB_PATH_MAX);
  /*printf("trying %s... ", path);*/
  if ((fp = fopen(path, "r")) != NULL) { 
    /*printf("found!\n");*/
    return(fp);
  }
  else {
    /*printf("not found\n");*/
    return(NULL);
  }
}


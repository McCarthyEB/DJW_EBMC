/******************************************************/
/** Subroutine to submit partially docked structures **/
/** to the clinton cluster ****************************/
/** Started AJWL 15/09/08  ****************************/
/******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void exxon_cluster()

{
#include "header.h"
char filename[FILELEN_MAX];

sprintf(filename,"run_zebedde%i",car_file_count);

if (!(exxon_submit_fp = fopen(filename,"w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening new Exxon cluster script file : %s\n",
                      filename);
        return;
    }

fprintf(exxon_submit_fp, "echo cd $PBS_O_WORKDIR\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_HOST =  $PBS_O_HOST \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_QUEUE =  $PBS_O_QUEUE \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_WORKDIR =  $PBS_O_WORKDIR\"\n",car_file_count);
fprintf(exxon_submit_fp, "echo \"PBS_ENVIRONMENT =  $PBS_ENVIRONMENT \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBID =  $PBS_JOBID \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBNAME =  $PBS_JOBNAME \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_NODEFILE =  $PBS_NODEFILE \"\n");
fprintf(exxon_submit_fp, "time echo \"PBS_QUEUE =  $PBS_QUEUE \"\n");
fprintf(exxon_submit_fp, "cat $PBS_NODEFILE\n");
fprintf(exxon_submit_fp, "echo \"running on node `hostname` \"\n");
fprintf(exxon_submit_fp, "cd $PBS_O_WORKDIR\n",car_file_count);
fprintf(exxon_submit_fp, "/home/alobo/Zebedde_EXE/zebedde %s\n" ,zebedde_input_filename);
fprintf(exxon_submit_fp, "exit\n");
fclose(exxon_submit_fp);
sprintf(exxon_submit_filename, filename);

return;
}

void disco_cluster(char *discover_root)
{
#include "header.h"

char filename[FILELEN_MAX];

sprintf(filename,"discover");

if (!(exxon_submit_fp = fopen(filename,"w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening new Exxon cluster script file \n",
                      filename);
        return;
    }

fprintf(exxon_submit_fp, "echo cd $PBS_O_WORKDIR\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_HOST =  $PBS_O_HOST \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_QUEUE =  $PBS_O_QUEUE \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_WORKDIR =  $PBS_O_WORKDIR\"\n");
fprintf(exxon_submit_fp, "echo \"PBS_ENVIRONMENT =  $PBS_ENVIRONMENT \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBID =  $PBS_JOBID \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBNAME =  $PBS_JOBNAME \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_NODEFILE =  $PBS_NODEFILE \"\n");
fprintf(exxon_submit_fp, "time echo \"PBS_QUEUE =  $PBS_QUEUE \"\n");
fprintf(exxon_submit_fp, "cat $PBS_NODEFILE\n");
fprintf(exxon_submit_fp, "echo \"running on node `hostname` \"\n");
fprintf(exxon_submit_fp, "cd $PBS_O_WORKDIR\n");
fprintf(exxon_submit_fp, "%s %s\n" ,discover_path,discover_root);
fprintf(exxon_submit_fp, "exit\n");
fclose(exxon_submit_fp);
sprintf(exxon_submit_filename, filename);

return;
}

void ucl_cluster()
{
#include "header.h"
char filename[FILELEN_MAX];

sprintf(filename,"run_zebedde%i",car_file_count);

if (!(exxon_submit_fp = fopen(filename,"w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening new UCL cluster script file \n",
                      filename);
        return;
    }

fprintf(exxon_submit_fp, "echo cd $PBS_O_WORKDIR\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_HOST =  $PBS_O_HOST \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_QUEUE =  $PBS_O_QUEUE \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_O_WORKDIR =  $PBS_O_WORKDIR\"\n");
fprintf(exxon_submit_fp, "echo \"PBS_ENVIRONMENT =  $PBS_ENVIRONMENT \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBID =  $PBS_JOBID \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_JOBNAME =  $PBS_JOBNAME \"\n");
fprintf(exxon_submit_fp, "echo \"PBS_NODEFILE =  $PBS_NODEFILE \"\n");
fprintf(exxon_submit_fp, "time echo \"PBS_QUEUE =  $PBS_QUEUE \"\n");
fprintf(exxon_submit_fp, "cat $PBS_NODEFILE\n");
fprintf(exxon_submit_fp, "echo \"running on node `hostname` \"\n");
fprintf(exxon_submit_fp, "cd $PBS_O_WORKDIR\n");
fprintf(exxon_submit_fp, "/home/alan/Zebedde_EXEs/zebedde64 $PBS_O_WORKDIR/%s\n" ,zebedde_input_filename);
fprintf(exxon_submit_fp, "exit\n");
fclose(exxon_submit_fp);
sprintf(exxon_submit_filename, filename);

return;
}


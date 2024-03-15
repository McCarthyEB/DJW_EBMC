#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "data.h"

void condor_submit()

{
#include "header.h"

char filename[FILELEN_MAX];

sprintf(filename,"zebedde%i.sub",car_file_count);

if (!(condor_submit_fp = fopen(filename,"w")))
    {
    fprintf(output_fp,"ZEBEDDE WARNING: Error opening new CONDOR input file: %s\n",
                      filename);
        return;
    }

fprintf(condor_submit_fp,"universe = vanilla\n");

fprintf(condor_submit_fp," Executable     = /home/alan/Zebedde_EXEs/zeb.$$(OpSys)\n\
Requirements = (Arch == \"INTEL\" && OpSys == \"LINUX\") || \\\n\
                (Arch == \"INTEL\" && OpSys == \"WINNT51\") || \\\n\
                (Arch == \"INTEL\" && OpSys == \"WINNT52\")\n");

fprintf(condor_submit_fp,"arguments = \"%s\"\n",zebedde_input_filename);

fprintf(condor_submit_fp,"input   = %s \n",zebedde_input_filename);

fprintf(condor_submit_fp,"  output  = zeb_screen.out.$(process)\n\
  error   = error.log.$(process)\n\
  Log     = log.log.$(process)\n");

fprintf(condor_submit_fp,"transfer_files = always\n");

fprintf(condor_submit_fp,"transfer_input_files = %s,%s,%s,%s,%s,%s\n", zebedde_input_filename, fragment_file, forcefield_library,
									temp_car_output, gulp_pots_file,pore_file);

fprintf(condor_submit_fp,"should_transfer_files = yes\n\n\
when_to_transfer_output = ON_EXIT_OR_EVICT\n\
getenv = true\n\n\
initial_dir = Run%i\n\
queue\n",car_file_count);
sprintf(condor_submit_filename, filename);
fclose(condor_submit_fp);
return;
}





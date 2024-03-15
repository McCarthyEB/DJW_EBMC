/******************************************************************/
/*** Print_statistics.c : Print rudimentary stats on a         ****/
/***                      make_a_template run                  ****/
/*** Started 7/5/96 DWL                                        ****/
/******************************************************************/
#include <stdio.h>
#include "maxima.h"
#include "structures.h"
#include "data.h"

void print_dashes(int ndashes, FILE *fp);

void print_statistics(stats *p_statistics, int number_of_actions, 
					    int *p_action_weights, FILE *output_fp)
{
double rate;
int i;
int tot_tries, tot_success;
double percent_good;
stats *p_action;

fprintf(output_fp,"Statistical Analysis of Run:\n");
print_dashes(40, output_fp);

fprintf(output_fp,"                              before   after      before after\n");
fprintf(output_fp,"                              Build    Build      Build  Build\n");
fprintf(output_fp,"Action              Weight    Tries    Tries      Accept Accept       %%Success\n");
print_dashes(80, output_fp);

tot_tries=0;
tot_success=0;

for (i=0;i<number_of_actions; i++)
	{
	p_action = p_statistics+i;

	tot_tries += (p_action->tries) + (p_action->tries_after_build);
	tot_success += (p_action->accepted) + (p_action->accepted_after_build);

    if ((p_action->tries) != 0)
		{
    	rate = ((p_action->accepted) *100) / (p_action->tries);
		}
	else
		{
		rate = 0;
		}

	fprintf(output_fp,"%-20s %5d    %5d    %5d    %5d    %5d      %6.1f\n",
                                                    action_list[i].name, 
						*(p_action_weights+i),
						p_action->tries, 
						p_action->tries_after_build, 
						p_action->accepted,	
						p_action->accepted_after_build,	
						rate);
	}
print_dashes(80, output_fp);
percent_good= tot_success/tot_tries*100;
fprintf(output_fp,"Totals                                %6d            %6d      %6.1f\n", 
								tot_tries, tot_success, percent_good);
print_dashes(80, output_fp);
fprintf(output_fp,"\n\n");

return;
}

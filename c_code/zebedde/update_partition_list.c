#include <stdio.h>
#include "maxima.h"
#include "structures.h"

 int update_partition_list(list_partition *p_list_partition, 
                           int index )
{

    p_list_partition->start= index;
    index += p_list_partition->num;
    p_list_partition->end= index;
    index++;

return index;
}

#include <stdio.h>

void print_product(double *p_matrix, double *p_vector, int dimension)
{
int iloop, jloop;

for (iloop=0; iloop < dimension; iloop++)
  {
     for (jloop=0    ; jloop < iloop    ; jloop++)
       {
         printf("           ");
       }
     for (jloop=iloop; jloop < dimension; jloop++)
       {
         printf("%10.6f ",*p_matrix); 
         p_matrix++;
       } 
     printf(" | %10.6f\n",*p_vector);
     p_vector++;
  }

return;
}

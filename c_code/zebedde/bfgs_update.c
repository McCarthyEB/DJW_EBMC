#include <stdio.h>
#include <time.h>

void matrix_x_ket( double *p_matrix, double *p_a, double *p_b,
                   int dimension );

void ket_bra_to_matrix( double *p_matrix, double *p_ai, double *p_bi,
                        double factor, int dimension );

double bra_x_ket( double *p_a, double *p_b, int dimension );

double bra_mat_ket( double *p_bra, double *p_matrix, double *p_ket, int dimension );

void print_product(double *p_matrix, double *p_vector, int dimension);

void bfgs_update(double *p_G, double *p_gamma, double *p_delta, int dimension)
{

/***********************************************************************/
/**** Test triangular matrix programs with arrays upto 10x10 ***********/
/***********************************************************************/

double temp_bra[10], temp_mat[55];

/***********************************************************************/
/**** Test times for triagular matrix and square matrix ****************/
/***********************************************************************/

double denominator, denominator_2, factor;

int iloop, dimension_2, triang_length; 

/*************************************************************************/
/***** BFGS update formula:                                  *************/
/*****                                                       *************/
/***** G => G+ (<d|g> + <g|G|g>) |d><d|                      *************/
/*****             |<d|g>| 2                                 *************/
/*****                                                       *************/
/*****       - ( |d><g|G + G|g><d| )                         *************/
/*****              <d|g>                                    *************/
/*************************************************************************/

/*************************************************************************/
/***** Work out the length of array for the triangular matrix ************/
/*************************************************************************/

/**************************************************************************/
/****** this should be calculated at the start of the minimisation !!******/

 dimension_2= dimension*dimension;
 triang_length= dimension/2 + dimension_2/2; 

 for (iloop=0; iloop < dimension; iloop++) temp_bra[iloop]=0;
 for (iloop=0; iloop < triang_length; iloop++) temp_mat[iloop]=0;

/*****odd case*****/
 if ((2*(dimension/2) -dimension) != 0) triang_length++;

/****** this should be calculated at the start of the minimisation !!******/
/**************************************************************************/

 denominator   = bra_x_ket(p_delta, p_gamma, dimension);
 denominator_2 =  denominator * denominator;

 printf("denominator= %10.6f denominator_2=%10.6f\n", denominator, denominator_2);
 matrix_x_ket(p_G, p_gamma, &temp_bra[0], dimension);

 printf("result of first matrix_x_ket operation: %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",
                     temp_bra[0],temp_bra[1],temp_bra[2],temp_bra[3],temp_bra[4],temp_bra[5]); 

 factor = denominator 
        + bra_mat_ket(p_gamma, p_G, p_gamma, dimension);

 factor = 0.5*factor/ denominator_2;

printf("Before first ket_bra_to_matrix:\n\n");
 print_product(p_G, p_delta, 6);

 ket_bra_to_matrix(&temp_mat[0], p_delta, p_delta, factor, dimension);

printf("After first ket_bra_to_matrix:\n\n");
 print_product(p_G, p_delta, 6);
 print_product(&temp_mat[0], p_delta, 6);

 factor = -1/ denominator;

 ket_bra_to_matrix(p_G, p_delta, &temp_bra[0], factor, dimension); 
printf("After second ket_bra_to_matrix:\n\n");
 print_product(p_G, p_delta, 6);


/***** Add in first bit **********/

 for (iloop=0; iloop < triang_length; iloop++)
   {
      *p_G += temp_mat[iloop];
      p_G++;
   }

return;
}

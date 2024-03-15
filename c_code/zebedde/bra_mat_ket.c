#include <stdio.h>
/*****************************************************************/
/****** To perform c = <bra|A|ket> where : ***********************/
/****** <bra| and |ket> are vectors A is a symmetric matrix ******/
/****** sent as its triangular segment, of arb. length dimension */
/*****************************************************************/

double bra_mat_ket( double *p_brai, double *p_matrix, double *p_keti,
                    int dimension )
{

int i,j;

double *p_braj, *p_ketj, ans;

ans= 0;
for (i=0; i < dimension; i++)
  {
    ans += *p_brai * *p_matrix * *p_keti;

    p_matrix++;

    p_braj= p_brai;
    p_ketj= p_keti;
    
    for (j=i+1; j < dimension; j++)
      { 
        p_braj++;
        p_ketj++;

        ans += *p_brai * *p_matrix * *p_ketj;      
        ans += *p_braj * *p_matrix * *p_keti;

        p_matrix++;

      }
    p_brai++;
    p_keti++;
  }
return ans;
}

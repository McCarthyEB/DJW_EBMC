/******************************************************************************/
/** index_sort.c : Sorts array indx on the contents of arrin                ***/
/**              : Lifted from indexx.c in Num Recip in C                   ***/
/**              : DWL 4/7/96                                               ***/
/******************************************************************************/
#include <stdio.h>
void index_sort(n,arrin,indx)
int n,indx[];
float arrin[];
{
	int l,j,ir,indxt,i;
	float q;

	for (j=1;j<=n;j++) indx[j]=j;
for (j=1;j<=n;j++) printf("index %d value %10.4f\n",indx[j], arrin[j]);
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
for (j=1;j<=n;j++) 
printf("index %d value %10.4f\n",indx[j], arrin[indx[j]]);
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}

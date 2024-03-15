#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "global_values.h"

/* subroutine for getting next integer and its sign SEPERATLY */

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign)
{
      int idave,k;
      int signif,sign_set;

      idave=0;

/* look for number */

      *point_itsa=FALSE;

/* don't lose first character */

      if (! *point_j) --(*point_j);

      while ( *point_j<=max_chars && !*point_itsa )

      {    ++*point_j;

/* look for sign */
           *sign=1;
           *point_itsa=TRUE;
           sign_set=0;

           if (*(p_ichar+*point_j) == '-' 
               && ( (*(p_ichar+*point_j+1) >= '0' && *(p_ichar+*point_j+1) <= '9')
                    || *(p_ichar+*point_j+1) == '.') )
           {
               *sign= -1;
               ++*point_j;
               sign_set=1;
           }

           if (*(p_ichar+*point_j) == '+'
                && ( (*(p_ichar+*point_j+1) >= '0' && *(p_ichar+*point_j+1) <= '9')
                    || *(p_ichar+*point_j+1) == '.') )
           {
               *sign= 1;
               ++*point_j;
               sign_set=1;
           }

/* look for digits */
           *point_ndigi=0;
           if (*(p_ichar+*point_j) < '0' || *(p_ichar+*point_j) > '9' ) 
                                                        *point_itsa= FALSE;

           while (*(p_ichar+*point_j) >= '0' && *(p_ichar+*point_j) <= '9'
                                             && *point_j < MAX_LINE_LEN )
           {
               ++*point_ndigi;
               ++*point_j;
           }

/* work out what number this string is and put in idave */
           idave=0;
           if ( *point_ndigi != 0 && *point_itsa )
           {
               signif=1;
               for (k=0; k<= *point_ndigi-1; ++k)
               {
                    idave= idave+ (*(p_ichar+*point_j-k-1)-'0')*signif;
                    signif=signif*10;
               }
           }

/* catch case where a minus sign is used without digits as in -.45
   for when this routine is used as part of reading doubles            */

           else if ( *point_ndigi == 0 && sign_set ) 
           {
                 idave=0;
                 *point_itsa= TRUE;
	   }
     }
      return idave;
}



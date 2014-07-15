/************************************************************************/
/**

   \file       pearson.c
   
   \version    V1.3
   \date       07.07.14
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1996-2014
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  29.01.96 Original   By: ACRM
-  V1.3  07.07.14 Use bl prefix for functions By: CTP


*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"


/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>REAL blpearson(REAL *x, REAL *y, int NItem)
   -------------------------------------------
*//**

   \param[in]     *x     Array of x items
   \param[in]     *y     Array of y items
   \param[in]     NItem  Number of items
   \return                 Pearson correlation coefficient

   This version makes 2 passes through the data

-  15.07.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP

*/
REAL blpearson(REAL *x, REAL *y, int NItem)
{
   REAL MeanX,
        MeanY,
        SumXSq,
        SumYSq,
        SumXY,
        xt, yt,
        r;
   int  i;

   /* Calculate means of x and y                                        */
   MeanX = MeanY = (REAL)0.0;
   for(i=0; i<NItem; i++)
   {
      MeanX += x[i];
      MeanY += y[i];
   }
   MeanX /= NItem;
   MeanY /= NItem;
   
   /* Calculate correlation coefficient                                 */
   SumXSq = SumYSq = SumXY = (REAL)0.0;
   for(i=0; i<NItem; i++)
   {
      xt = x[i] - MeanX;
      yt = y[i] - MeanY;
      SumXSq += (xt * xt);
      SumYSq += (yt * yt);
      SumXY  += (xt * yt);
   }

   r = SumXY / (REAL)sqrt((double)(SumXSq * SumYSq));
   
   return(r);
}



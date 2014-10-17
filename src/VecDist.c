/************************************************************************/
/**

   \file       VecDist.c
   
   \version    V1.3
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996-2014
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
-  V1.1  18.06.96 Added vector routines
-  V1.2  06.10.98 Added VecAdd3()
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Vectors
   #FUNCTION blVecDist()
   Finds the distance between two vectors of arbitrary length
*/
/************************************************************************/
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
/*>REAL blVecDist(REAL *a, REAL *b, int len)
   -----------------------------------------
*//**

   \param[in]     *a     An arbitrary length vector (as an array)
   \param[in]     *b     An arbitrary length vector (as an array)
   \param[in]     len    The dimensionality of the vectors (array
                           length)
   \return                   The distance between the points described by
                           the two vectors

   Finds the distance between two vectors of arbitrary length

-  28.07.95 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blVecDist(REAL *a, REAL *b, int len)
{
   REAL sumsq = 0.0;
   int  i;
   
   for(i=0; i<len; i++)
      sumsq += (a[i] - b[i]) * (a[i] - b[i]);

   return((REAL)sqrt((double)sumsq));
}



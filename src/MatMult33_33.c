/************************************************************************/
/**

   \file       MatMult33_33.c
   
   \version    V1.7
   \date       07.07.14
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, University of Reading, 2002-14
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
-  V1.0  06.09.91 Original
-  V1.0a 01.06.92 Documented
-  V1.1  30.09.92 Matrix multiplication added
-  V1.2  10.06.93 void return from matrix multiplication
-  V1.3  22.07.93 Added CreateRotMat()
-  V1.4  03.08.93 Changed matrix multiplication to standard direction
-  V1.5  28.07.95 Added VecDist()
-  V1.6  27.09.95 Added MatMult33_33()
-  V1.7  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
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
/*>void blMatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
   ---------------------------------------------------------------
*//**

   \param[in]     a            Matrix to be multiplied
   \param[in]     b            Matrix to be multiplied
   \param[out]    out          Output matrix

   Multiply two 3x3 matrices

-  27.09.95 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blMatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
{
   int  i, j, k;
   REAL ab;
   
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         ab = (REAL)0.0;
         for(k=0; k<3; k++)
         {
            ab += a[i][k]*b[k][j];
         }
         out[i][j]=ab;
      }
   }
}

         
   
      
   

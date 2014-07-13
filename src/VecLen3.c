/************************************************************************/
/**

   \file       VecLen3.c
   
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
/* Includes
*/
#include <math.h>
#include "MathUtil.h"

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
/*>REAL blVecLen3(VEC3F Vec)
   -------------------------
*//**

   \param[in]     Vec       Vector
   \return                      Length of vector

-  18.06.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blVecLen3(VEC3F Vec)
{
   return((REAL)sqrt((double)((Vec.x * Vec.x) + 
                              (Vec.y * Vec.y) + 
                              (Vec.z * Vec.z))));
}



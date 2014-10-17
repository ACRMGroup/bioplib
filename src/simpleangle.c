/************************************************************************/
/**

   \file       simpleangle.c
   
   \version    V1.6
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin, 1993-2014
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

-  V1.6  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Geometry
   #FUNCTION  blSimpleangle()
   Simplifies a signed angle to an unsigned angle <=2*PI
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"
#include "macros.h"

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
/*>REAL blSimpleangle(REAL ang)
   ----------------------------
*//**

   \param[in]     ang         An angle
   \return                    Simplified angle
   
   Simplifies a signed angle to an unsigned angle <=2*PI

-  07.02.89 Original    By: ACRM
-  04.03.91 Fixed return value
-  16.06.93 Changed float to REAL
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blSimpleangle(REAL ang)
{
   /* Reduce to less than 360 degrees                                   */
   while(ang > 2*PI) ang -= 2*PI;
   
   if(ang >= 0.0 && ang <= PI)         /* 1st & 2nd quadrant +ve        */
      return(ang);
   else if(ang > PI)                   /* 3rd & 4th quadrant +ve        */
      ang = 2*PI - ang;
   else if(ang < 0.0 && ang > -PI)     /* 1st & 2nd quadrant -ve        */
      ang = ABS(ang);
   else                                /* 3rd & 4th quadrant -ve        */
      ang = 2*PI + ang;
   
   return(ang);
}


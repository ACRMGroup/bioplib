/************************************************************************/
/**

   \file       TrueAngle.c
   
   \version    V1.6
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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

   REAL TrueAngle(REAL opp, REAL adj)

   \param[in]     opp         Length of opposite side
   \param[in]     adj         Length of adjacent side
   \return                    The angle from 0 to 2PI

   Returns the true positive angle between 0 and 2PI given the opp and
   adj lengths

**************************************************************************

   Revision History:
   =================
-  V1.0  07.02.91 Original
-  V1.1  17.02.91 Corrected comments to new standard and added phi()
-  V1.2  04.03.91 angle() and phi() now return _correct_ values!
-  V1.3  01.06.92 ANSIed
-  V1.4  08.12.92 Changed abs() to ABS() from macros.h
-  V1.5  27.03.95 Added TrueAngle()
-  V1.6  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Geometry
   #ROUTINE  blTrueAngle()
   Return the +ve angle between 0 and 2PI given the opp and adj values.
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
/*>REAL blTrueAngle(REAL opp, REAL adj)
   ------------------------------------
*//**

   \param[in]     opp     Opposite length
   \param[in]     adj     Adjacent length
   \return                     Angle between 0 and 2PI

   Return the +ve angle between 0 and 2PI given the opp and adj values.

-  25.07.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blTrueAngle(REAL opp, REAL adj)
{
   REAL ang;
   
   if(adj != 0.0)
   {
      ang = (REAL)atan((double)(opp/adj));

      /* 4th quadrant; ang -ve so add 2PI                             */
      if(opp < 0.0 && adj > 0.0) ang += 2*PI;

      /* 2nd & 3rd quadrant; add PI                                     */
      if(adj < 0.0) ang += PI;
   }
   else
   {
      if(opp>0.0)                /* 1st->2nd quadrant boundary          */
         ang = PI/2.0;
      else                       /* 3rd->4th quadrant boundary          */
         ang = 3.0*PI/2.0;
   }
   
   if(opp == 0.0)
   {
      if(adj > 0.0)              /* 4th->1st quadrant boundary          */
         ang = 0.0;
      else                       /* 2nd->3rd quadrant boundary          */
         ang = PI;
   }

   return(ang);
}


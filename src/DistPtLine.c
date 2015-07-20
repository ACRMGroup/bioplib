/************************************************************************/
/**

   \file       DistPtLine.c
   
   \version    V1.5
   \date       20.07.15
   \brief      General maths/stats/vector functions
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996-2015
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
-  V1.4  20.07.15 Renamed from blDistPtVect()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Geometry
   #FUNCTION  blDistPtLine()
   Calculate the distance from a point to a line vector described by two
   vector end points. See also blPointLineDistance()
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"
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
/*>REAL blDistPtLine(VEC3F Point, VEC3F End1, VEC3F End2)
   ------------------------------------------------------
*//**

   \param[in]     Point     The coordinates of a point
   \param[in]     End1      Coordinates of one end of vector
   \param[in]     End2      Coordinates of other end of vector
   \return                  The distance from pt to line

   Calculate the distance from a point to a line described by two
   vector end points

-  18.06.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
-  20.07.15 Renamed from blDistPtVect()  By: ACRM
*/
REAL blDistPtLine(VEC3F Point, VEC3F End1, VEC3F End2)
{
   VEC3F Vec,
         UVec,
         PQVec,
         PRVec;
   REAL  len;
   
   /* Find the vector from End1 to End2                                 */
   blVecSub3(&Vec, End2, End1);

   /* Find the length of this vector                                    */
   len = blVecLen3(Vec);

   /* Now calculate the unit vector                                     */
   UVec.x = Vec.x / len;
   UVec.y = Vec.y / len;
   UVec.z = Vec.z / len;

   /* Calculate the vector from the point to an arbitrary point on the
      line (we'll choose End1)
   */
   PQVec.x = End1.x - Point.x;
   PQVec.y = End1.y - Point.y;
   PQVec.z = End1.z - Point.z;

   /* PRVect is the cross product of PQVect with the unit vector        */
   blCrossProd3(&PRVec, PQVec, UVec);

   /* The length we want is the length of this vector                   */
   return(blVecLen3(PRVec));
}

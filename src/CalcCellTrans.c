/************************************************************************/
/**

   \file       CalcCellTrans.c
   
   \version    V1.1
   \date       07.07.14
   \brief      Calculate offsets for creating a crystal lattice
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-1995
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
-  V1.0R 12.10.05 Original
-  V1.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Defines required for includes
*/

/************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Calculations
   #ROUTINE  blCalcCellTrans()
   Calculates the offsets to apply in X, Y and Z directions for creating
   a crystal lattice from the unit cell parameters.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "fsscanf.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void blCalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                        VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
   -----------------------------------------------------------------
*//**

   \param[in]     UnitCell       The unit cell dimensions
   \param[in]     CellAngles     The unit cell angles
   \param[out]    *xtrans        Translation to apply along X axis
   \param[out]    *ytrans        Translation to apply along Y axis
   \param[out]    *ztrans        Translation to apply along Z axis

   Calculates the offsets to apply in X, Y and Z directions for creating
   a crystal lattice from the unit cell parameters.

-  11.10.95 Original    By: ACRM, Based losely on code from Rasmol
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blCalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                     VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
{
   REAL lena, lenb, lenc,
        cosa, cosb, cosg, sing,
        tmpx, tmpy, tmpz, temp;
   
   lena =  UnitCell.x;                         /* A                     */
   lenb =  UnitCell.y;                         /* B                     */
   lenc = -1.0*UnitCell.z;                     /* C                     */

   cosa = (REAL)cos((double)(CellAngles.x));   /* Alpha                 */
   cosb = (REAL)cos((double)(CellAngles.y));   /* Beta                  */
   cosg = (REAL)cos((double)(CellAngles.z));   /* Gamma                 */
   sing = (REAL)sin((double)(CellAngles.z));   /* Gamma                 */

   temp = cosa*cosa + cosb*cosb + cosg*cosg - 2.0*cosa*cosb*cosg;
   tmpx = cosb; 
   tmpy = (cosa - cosb*cosg)/sing;
   tmpz = (REAL)sqrt((double)(1.0-temp))/sing;
   
   xtrans->x = lena;
   xtrans->y = (REAL)0.0;
   xtrans->z = (REAL)0.0;
   
   ytrans->x = lenb*cosg;
   ytrans->y = lenb*sing;
   ytrans->z = (REAL)0.0;
   
   ztrans->x = lenc*tmpx;
   ztrans->y = lenc*tmpy;
   ztrans->z = lenc*tmpz;
}



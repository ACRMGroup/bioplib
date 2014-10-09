/************************************************************************/
/**

   \file       angle.h
   
   \version    V1.10
   \date       14.08.14
   \brief      Include file for angle functions
   
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

**************************************************************************

   Revision History:
   =================
-  V1.5  27.03.95 
-  V1.6  09.07.96 Added TorToCoor()
-  V1.7  06.09.96 Added includes
-  V1.8  07.07.14 Use bl prefix for functions By: CTP
-  V1.9  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.10 14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP

*************************************************************************/
/* Includes
*/
#include "MathType.h"
#include "SysDefs.h"

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

#ifndef _ANGLE_H
#define _ANGLE_H

REAL blAngle(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
             REAL xk, REAL yk, REAL zk);
REAL blPhi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
           REAL xk, REAL yk, REAL zk, REAL xl, REAL yl, REAL zl);
REAL blSimpleangle(REAL ang);
REAL blTrueAngle(REAL opp, REAL adj);
BOOL blTorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, 
                 REAL bond, REAL theta, REAL torsion,
                 VEC3F *coords);



/************************************************************************/
/* Include deprecated functions                                         */
#define _ANGLE_H_DEPRECATED
#include "deprecated.h" 
/************************************************************************/


#endif

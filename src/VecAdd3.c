/************************************************************************/
/**

   \file       VecAdd3.c
   
   \version    V1.2
   \date       06.10.98
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996-8
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
/*>void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2)
   ----------------------------------------------
*//**

   \param[in]     In1       First vector
   \param[in]     In2       Second vector
   \param[out]    Out       Output vector

   Add 2 vectors
-  06.10.98 Original   By: ACRM
*/
void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   Out->x = In1.x + In2.x;
   Out->y = In1.y + In2.y;
   Out->z = In1.z + In2.z;
}



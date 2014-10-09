/************************************************************************/
/**

   \file       RotPDB.c
   
   \version    V1.1
   \date       August 1993
   \brief      Rotate a PDB linked list
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993
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

   Rotate a PDB linked list. Moves the structure to the origin first,
   applies the rotation and moves back from the origin.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/* Variables global to this file only
*/


/************************************************************************/
/*>void blRotatePDB(PDB *pdb, REAL matrix[3][3])
   ---------------------------------------------
*//**

   \param[in,out] *pdb          PDB linked list to rotate
   \param[in]     matrix        Rotation matrix

   Rotates a PDB linked list using ApplyMatrixPDB() which ignores 
   coordinates of 9999.0. The structure is moved to the origin, the 
   matrix is applied and the structure is moved back.

-  30.09.92 Original
-  01.10.92 Added check on NULL coordinates
-  22.07.93 Moves to origin first; calls ApplyMatrixPDB() to do the work
-  07.07.14 Renamed to blRotatePDB(). Use bl prefix for functions. By: CTP
*/
void blRotatePDB(PDB *pdb, REAL matrix[3][3])
{
   VEC3F CofG;
         
   blGetCofGPDB(pdb, &CofG);
   blOriginPDB(pdb);
   
   blApplyMatrixPDB(pdb, matrix);
   
   blTranslatePDB(pdb, CofG);
}


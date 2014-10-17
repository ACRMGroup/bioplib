/************************************************************************/
/**

   \file       ApMatPDB.c
   
   \version    V1.1
   \date       07.07.14
   \brief      
   
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Moving the structure
   #ROUTINE  blApplyMatrixPDB()
   Apply a rotation matrix to a PDB linked list.
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"
#include "pdb.h"
#include "matrix.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/* Variables global to this file only
*/

/************************************************************************/
/*>void blApplyMatrixPDB(PDB *pdb, REAL matrix[3][3])
   --------------------------------------------------
*//**

   \param[in,out] *pdb          PDB linked list
   \param[in]     matrix        Matrix to apply

   Apply a rotation matrix to a PDB linked list.

-  22.07.93 Original (old RotatePDB())   By: ACRM
-  07.07.14 Renamed to blApplyMatrixPDB() By: CTP
*/
void blApplyMatrixPDB(PDB  *pdb,
                      REAL matrix[3][3])
{
   PDB   *p;
   VEC3F incoords,
         outcoords;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x != 9999.0 && p->y != 9999.0 && p->z != 9999.0)
      {
         incoords.x = p->x;
         incoords.y = p->y;
         incoords.z = p->z;
         blMatMult3_33(incoords,matrix,&outcoords);
         p->x = outcoords.x;
         p->y = outcoords.y;
         p->z = outcoords.z;
      }
   }
}


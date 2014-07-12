/************************************************************************/
/**

   \file       OriginPDB.c
   
   \version    V1.2
   \date       07.07.14
   \brief      Move a PDB linked list to the origin
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-4
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
-  V1.0  01.10.92 Original
-  V1.1  22.02.94 Changed NULL check to any coordinate not 9999.0
-  V1.2  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#include <math.h>

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
/*>void blOriginPDB(PDB *pdb)
   --------------------------
*//**

   \param[in,out] *pdb    PDB linked list to move

   Moves a PDB linked list to the origin, ignoring NULL coordinates.

-  01.10.92 Original
-  22.02.94 Changed NULL check to any coordinate not 9999.0
-  11.03.94 Changed NULL check to >9998.0 Added cast to REAL
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blOriginPDB(PDB *pdb)
{
   PDB   *p;
   VEC3F cg;

   blGetCofGPDB(pdb,&cg);

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x < (REAL)9999.0 || p->y < (REAL)9999.0 || p->z < (REAL)9999.0)
      {
         p->x -= cg.x;
         p->y -= cg.y;
         p->z -= cg.z;
      }
   }
}


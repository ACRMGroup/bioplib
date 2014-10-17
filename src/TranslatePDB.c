/************************************************************************/
/**

   \file       TranslatePDB.c
   
   \version    V1.3
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-8
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
-  V1.1  01.03.94 Original
-  V1.2  27.02.98 Removed unreachable break from switch()
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Moving the structure
   #FUNCTION  blTranslatePDB()
   Translates a PDB linked list, ignoring null (9999.0) coordinates.
*/
/************************************************************************/
/* Includes
*/
#include "MathType.h"
#include "pdb.h"
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
/*>void blTranslatePDB(PDB *pdb,VEC3F tvect)
   -----------------------------------------
*//**

   \param[in,out] *pdb   PDB linked list to move
   \param[in]     tvect  Translation vector

   Translates a PDB linked list, ignoring null (9999.0) coordinates.

-  01.10.92 Original
-  11.03.94 Changed check on 9999.0 to >9998.0 and cast to REAL
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blTranslatePDB(PDB   *pdb,
                    VEC3F tvect)
{
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x < (REAL)9999.0 && p->y < (REAL)9999.0 && p->z < (REAL)9999.0)
      {
         p->x += tvect.x;
         p->y += tvect.y;
         p->z += tvect.z;
      }
   }
}


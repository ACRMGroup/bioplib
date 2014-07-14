/************************************************************************/
/**

   \file       CalcChiPDB.c
   
   \version    V1.2
   \date       07.07.14
   \brief      Perform calculations on PDB linked list
   
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
-  V1.0  22.02.94 Original
-  V1.1  27.02.98 Removed unreachable break from switch()
-  V1.2  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#include <math.h>

#include "MathType.h"
#include "macros.h"
#include "pdb.h"
#include "angle.h"

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
/*>REAL blCalcChi(PDB *pdb, int type)
   ----------------------------------
*//**

   \param[in]     *pdb     PDB linked list
   \param[in]     type     Torsion type (see below)
   \return                   Torsion angle

   Calculates a sidechain torsion angle from a pdb linked list. The atoms
   to be included in the calculation are specified by type.


         type     Atom names        Sequential atom numbers
         --------------------------------------------------
         0        N,  CA, CB, XG    (0 - 1 - 4 - 5)
         1        CA, CB, XG, XD    (1 - 4 - 5 - 6)
         2        CB, XG, XD, XE    (4 - 5 - 6 - 7)
         3        XG, XD, XE, XZ    (5 - 6 - 7 - 8)

-  13.05.92 Original
-  27.02.98 Removed unreachable break from switch()
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blCalcChi(PDB *pdb,
               int type)
{
   REAL  chi = 0.0;
   PDB   *one,
         *two,
         *three,
         *four;
   
   switch(type)
   {
   case 0:              /* N,  CA, CB, XG    (0 - 1 - 4 - 5)            */
      one   = blGetPDBByN(pdb, 0);
      two   = blGetPDBByN(pdb, 1);
      three = blGetPDBByN(pdb, 4);
      four  = blGetPDBByN(pdb, 5);
      break;
   case 1:              /* CA, CB, XG, XD    (1 - 4 - 5 - 6)            */
      one   = blGetPDBByN(pdb, 1);
      two   = blGetPDBByN(pdb, 4);
      three = blGetPDBByN(pdb, 5);
      four  = blGetPDBByN(pdb, 6);
      break;
   case 2:              /* CB, XG, XD, XE    (4 - 5 - 6 - 7)            */
      one   = blGetPDBByN(pdb, 4);
      two   = blGetPDBByN(pdb, 5);
      three = blGetPDBByN(pdb, 6);
      four  = blGetPDBByN(pdb, 7);
      break;
   case 3:              /* XG, XD, XE, XZ    (5 - 6 - 7 - 8)            */
      one   = blGetPDBByN(pdb, 5);
      two   = blGetPDBByN(pdb, 6);
      three = blGetPDBByN(pdb, 7);
      four  = blGetPDBByN(pdb, 8);
      break;
   default:
      return(chi);
   }
   
   /* Calculate the torsion angle                                       */
   chi = blphi(one->x,   one->y,   one->z,
               two->x,   two->y,   two->z,
               three->x, three->y, three->z,
               four->x,  four->y,  four->z);
   
   return(chi);
}


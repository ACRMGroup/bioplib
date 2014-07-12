/************************************************************************/
/**

   \file       RenumAtomsPDB.c
   
   \version    V1.3
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

**************************************************************************

   Revision History:
   =================
-  V1.1  01.03.94 Original
-  V1.2  27.02.98 Removed unreachable break from switch()
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#include "macros.h"
#include "pdb.h"

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
/*>void blRenumAtomsPDB(PDB *pdb)
   ------------------------------
*//**

   \param[in,out] *pdb   PDB linked list to renumber

   Renumber the atoms throughout a PDB linked list

-  01.08.93 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blRenumAtomsPDB(PDB *pdb)
{
   PDB *p;
   int i;

   for(p=pdb, i=1; p!=NULL; NEXT(p)) 
      p->atnum=i++;
}


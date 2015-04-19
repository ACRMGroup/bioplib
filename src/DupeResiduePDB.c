/************************************************************************/
/**

   \file       DupeResiduePDB.c
   
   \version    V1.2
   \date       07.07.14
   \brief      Create a new PDB linked list with a copy of a residue
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1996-2007
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
-  V1.0 27.08.96 Original from mutmodel  By: ACRM
-  V1.1 08.11.07 Initialize p and q; Moved into bioplib
-  V1.2 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blDupeResiduePDB()
   Makes a duplicate PDB linked list of just one residue
   Note that CONECT data will not be preserved since it would not
   be valid.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
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
/*>PDB *blDupeResiduePDB(PDB *in)
   ----------------------------
*//**

   \param[in]     *in     PDB linked list pointing to residue to duplicate
   \return                Duplicate linked list of residue at `in'
                          (NULL if allocation fails)

   Makes a duplicate PDB linked list of just the residue pointed to by
   `in'. Note that CONECT data will not be preserved since it would not
   be valid.

-  27.08.96 Original   By: ACRM
-  08.11.07 Initialize p and q
            Moved into bioplib
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDB *blDupeResiduePDB(PDB *in)
{
   PDB *pNext,
       *out = NULL,
       *p = NULL, 
       *q = NULL;
   
   /* Find the next residue                                             */
   pNext = blFindNextResidue(in);

   for(p=in; p!=pNext; NEXT(p))
   {
      if(out==NULL)
      {
         INIT(out, PDB);
         q=out;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      if(q==NULL)
      {
         FREELIST(out, PDB);
         return(NULL);
      }
      
      blCopyPDB(q, p);
      q->nConect=0;
   }
   
   return(out);
}



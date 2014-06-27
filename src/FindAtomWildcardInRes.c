/************************************************************************/
/**

   \file       FindAtomWildcardInRes.c
   
   \version    V1.0
   \date       27.08.96
   \brief      Find an atom within a residue allowing wild cards
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1996
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
-  V1.0  27.08.96 Original moved from mutmodel  By: ACRM

*************************************************************************/
/* Includes
*/
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
/*>PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern)
   ---------------------------------------------------
   Input:   PDB   *pdb      Pointer to start of a residue
            char  *pattern  Atom name pattern to find
   Returns: PDB   *         Pointer to requested atom or NULL if not
                            found.

   Finds an atom within the residue given as a PDB pointer. Allows 
   single character wildcards. Thus ?G? maybe used for any atom at the
   gamma position.

   Returns the first atom which matches

   27.08.96 Original   By: ACRM
*/
PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern)
{
   PDB  *p, *pNext;
   int  i;
   BOOL ok;
   
   pNext = FindNextResidue(pdb);
   
   for(p=pdb; p!=pNext; NEXT(p))
   {
      if(!strncmp(p->atnam,pattern,4))
         return(p);
      ok = TRUE;
      for(i=0; i<4; i++)
      {
         if((pattern[i] != p->atnam[i]) && pattern[i] != '?')
         {
            ok = FALSE;
            break;
         }
      }
      if(ok)
         return(p);
   }
   return(NULL);
}




/*************************************************************************

   Program:    
   File:       ExtractZonePDB.c
   
   Version:    V1.10R
   Date:       08.10.99
   Function:   PDB linked list manipulation
   
   Copyright:  (c) SciTech Software 1992-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
   V1.0  22.02.94 Original release
   V1.1  23.05.94 Added FindNextChainPDB()
   V1.2  05.10.94 KillSidechain() uses BOOL rather than int
   V1.3  24.07.95 Added TermPDB()
   V1.4  25.07.95 Added GetPDBChainLabels()
   V1.5  26.09.95 Fixed bug in TermPDB()
   V1.6  12.10.95 Added DupePDB(), CopyPDBCoords()
   V1.7  23.10.95 Moved FindResidueSpec() to ParseRes.c
   V1.8  10.01.96 Added ExtractZonePDB()
   V1.9  14.03.96 Added FindAtomInRes()
   V1.10 08.10.99 Initialised some variables

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"
#include "general.h"

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
/*>PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, 
                       char *insert1, char *chain2, int resnum2, 
                       char *insert2)
   -----------------------------------------------------------------------
   Input:   PDB    *inpdb   Input PDB linked list
            char   *chain1  Start residue chain name
            int    resnum1  Start residue number
            char   *insert1 Start residue insert code
            char   *chain2  End residue chain name
            int    resnum2  End residue number
            char   *insert2 End residue insert code
   Returns: PDB    *        PDB linked list of the region of interest.

   Reduces a PDB linked list to those residues within a specified zone.
   Note that the PDB linked list is duplicated before extraction so
   pointers do not match those in the input PDB linked list. Excess
   records in the new PDB linked list are freed.

   10.01.96 Original   By: ACRM
*/
PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                    char *chain2, int resnum2, char *insert2)
{
   PDB *pdb, *p, 
       *start = NULL, 
       *last  = NULL, 
       *prev  = NULL;

   /* Duplicate the PDB linked list                                     */
   if((pdb = DupePDB(inpdb))==NULL)
      return(NULL);

   /* Find the first residue in the PDB linked list                     */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    == resnum1)   &&
         (p->chain[0]  == chain1[0]) &&
         (p->insert[0] == insert1[0]))
      {
         start = p;
         break;
      }
      prev = p;
   }

   if(start==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Find the last residue in the PDB linked list                      */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    == resnum2)   &&
         (p->chain[0]  == chain2[0]) &&
         (p->insert[0] == insert2[0]))
      {
         last = p;
         break;
      }
   }

   if(last==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Step last onto the final atom in that residue                     */
   for(; last->next!=NULL; NEXT(last))
   {
      if((last->next->resnum    != resnum2)   ||
         (last->next->chain[0]  != chain2[0]) ||
         (last->next->insert[0] != insert2[0]))
      {
         break;
      }
   }

   /* Free linked list after 'last'                                     */
   if(last->next != NULL)
   {
      FREELIST(last->next, PDB);
      last->next = NULL;
   }
   
   /* Unlink 'start' from rest of linked list and free memory before 
      'start'
   */
   if(prev != NULL)
   {
      prev->next = NULL;
      FREELIST(pdb, PDB);
   }

   return(start);
}



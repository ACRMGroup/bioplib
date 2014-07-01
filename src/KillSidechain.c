/************************************************************************/
/**

   \file       KillSidechain.c
   
   \version    V1.10
   \date       08.10.99
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-6
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
-  V1.0  22.02.94 Original release
-  V1.1  23.05.94 Added FindNextChainPDB()
-  V1.2  05.10.94 KillSidechain() uses BOOL rather than int
-  V1.3  24.07.95 Added TermPDB()
-  V1.4  25.07.95 Added GetPDBChainLabels()
-  V1.5  26.09.95 Fixed bug in TermPDB()
-  V1.6  12.10.95 Added DupePDB(), CopyPDBCoords()
-  V1.7  23.10.95 Moved FindResidueSpec() to ParseRes.c
-  V1.8  10.01.96 Added ExtractZonePDB()
-  V1.9  14.03.96 Added FindAtomInRes()
-  V1.10 08.10.99 Initialised some variables

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

#include "SysDefs.h"
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
/*>BOOL KillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB)
   ----------------------------------------------------------
*//**

   \param[in]     *ResStart     Start of a residue in linked list
   \param[in]     *NextRes      Start of next residue
   \param[in]     doCB          Flag to kill CB as part of s/c
   \return                        Success?
   
   Kill a sidechain, by calls to KillPDB(). If doCB is set, will kill 
   the CB.
   N.B. At least 1 backbone atom must occur in the linked list before the
   sidechain.
   
-  12.05.92 Original
-  05.10.94 doCB is now a BOOL as is the return
*/
BOOL KillSidechain(PDB *ResStart,   /* Pointer to start of residue      */
                   PDB *NextRes,    /* Pointer to start if next residue */
                   BOOL doCB)       /* Flag to kill the CB              */
{
   PDB *p,
       *prev = NULL;
   
   for(p=ResStart; p && p!=NextRes; NEXT(p))
   {
      if(strcmp(p->atnam, "N   ") &&
         strcmp(p->atnam, "CA  ") &&
         strcmp(p->atnam, "C   ") &&
         strcmp(p->atnam, "O   ") &&
         strcmp(p->atnam, "CB  "))
      {
         if(prev == NULL) return(FALSE); /* No b/b atom before s/c      */

         /* KillPDB() returns the next in list, so exit if list ended   */
         if(KillPDB(p, prev) == NULL) break;
         p = prev;
      }
      
      /* Kill the CB if required                                        */
      if(doCB && !strcmp(p->atnam,"CB  "))
      {
         if(prev == NULL) return(FALSE);  /* No b/b atom before s/c     */

         /* KillPDB() returns the next in list, so exit if list ended   */
         if(KillPDB(p, prev) == NULL) break;
         p = prev;
      }
      
      prev = p;
   }
   return(TRUE);
}


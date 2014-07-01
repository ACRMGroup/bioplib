/************************************************************************/
/**

   \file       KillPDB.c
   
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
/*>PDB *KillPDB(PDB *pdb, PDB *prev)
   ---------------------------------
*//**

   \param[in]     *pdb    Pointer to item in PDB linked list to be removed
   \param[in]     *prev   Pointer to previous item in linked list
   \return                Next item in PDB linked list

   Kill an item in the PDB linked list and re-link correctly. Returns the
   next item in the list, so will be NULL when the last item in the list
   is killed.

-  12.05.92 Original
-  11.03.94 Now handles prev==NULL to delete first item in a list
*/
PDB *KillPDB(PDB *pdb,              /* Pointer to record to kill        */
             PDB *prev)             /* Pointer to previous record       */
{
   PDB *p;

/* Old action was just to return if prev==NULL
   if(prev == NULL) return(NULL);
*/
   p = pdb->next;

   if(prev!=NULL)
      prev->next = pdb->next;       /* Relink the list                  */
   free(pdb);                       /* Free the item                    */

   return(p);
}


/************************************************************************/
/**

   \file       FindRawAtomInRes.c
   
   \version    V1.13
   \date       07.07.14
   \brief      PDB linked list manipulation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2014
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
-  V1.11 28.02.01 Added FindRawAtomInRes()
-  V1.12 03.06.05 Compares 4 rather than 5 characters
-  V1.13 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Searching the PDB linked list        
   #FUNCTION  blFindRawAtomInRes()
   Searches the raw atom name (atnam_raw) field of the current residue
   for the specified atom name
*/
/************************************************************************/
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
/*>PDB *blFindRawAtomInRes(PDB *pdb, char *atnam_in)
   -------------------------------------------------
*//**

   \param[in]     *pdb         The beginning of a residue in a PDB 
                               linked list
   \param[in]     *atnam_in    An atom name to search for (doesn't need
                               to be space-padded)
   \return                     Pointer to required atom, NULL if not
                               found

   Searches the raw atom name (atnam_raw) field for the specified atom
   name

-  28.02.01 Original based on FindAtomInRes()  By: ACRM
-  03.06.05 Now compares 4 characters rather than 5
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDB *blFindRawAtomInRes(PDB *pdb, char *atnam_in)
{
   PDB *end,
       *p;

   char atnam[8];
   
   /* First copy the specified atom name and pad to 5 chars             */
   strcpy(atnam,atnam_in);
   blPadterm(atnam,5);
   
   /* Find the end of this residue                                      */
   end = blFindNextResidue(pdb);
   
   /* Search for the required atom                                      */
   for(p=pdb; p!=end; NEXT(p))
   {
      if(!strncmp(p->atnam_raw,atnam,4))
         return(p);
   }
   
   return(NULL);
}

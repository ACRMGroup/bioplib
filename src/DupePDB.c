/************************************************************************/
/**

   \file       DupePDB.c
   
   \version    V1.12
   \date       19.04.15
   \brief      PDB linked list manipulation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2015
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
-  V1.11 07.07.14 Use bl prefix for functions By: CTP
-  V1.12 19.04.15 Added call to blCopyConect()   By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blDupePDB()
   Duplicates a PDB linked list. CONECT data are updated to point within
   the new list.
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
/*>PDB *blDupePDB(PDB *in)
   -----------------------
*//**

   \param[in]     *in     Input PDB linked list
   \return                Duplicated PDB linked list
                          (NULL on allocation failure)

   Duplicates a PDB linked list. Allocates new linked list with identical
   data. CONECT data are updated to point within the new list.

-  11.10.95 Original   By: ACRM
-  08.10.99 Initialise q to NULL
-  07.07.14 Use bl prefix for functions By: CTP
-  19.04.15 Added call to blCopyConect()   By: ACRM
*/
PDB *blDupePDB(PDB *in)
{
   PDB *out = NULL,
       *p, *q = NULL;

   for(p=in; p!=NULL; NEXT(p))
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
   }
   
   if(!blCopyConects(out, in))
   {
      FREELIST(out, PDB);
      return(NULL);
   }

   return(out);
}



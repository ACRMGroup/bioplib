/************************************************************************/
/**

   \file       TermPDB.c
   
   \version    V1.12
   \date       07.07.14
   \brief      
   
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
-  V1.11 04.02.14 Use CHAINMATCH macro. By: CTP
-  V1.12 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #ROUTINE  blTermPDB()
   Terminate a PDB linked list after length residues, returning a pointer
   to the next residue. 
*/
/************************************************************************/
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
/*>PDB *blTermPDB(PDB *pdb, int length)
   ----------------------------------
*//**

   \param[in]     *pdb         PDB linked list
   \param[in]     length       Number of residues after which to terminate
   \return                     Pointer to next residue after terminated
                               list. NULL if not enough residues in linked
                               list.

   Terminate a PDB linked list after length residues, returning a pointer
   to the next residue. 

   Note that the number of residues may cross chain boundaries.

-  06.07.95 Original    By: ACRM
-  26.09.95 Corrected update of resnum etc to use p-> not pdb-> (!!)
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDB *blTermPDB(PDB *pdb, int length)
{
   int  resnum,
        count;
   char insert,
        chain[8];
   PDB  *p,
        *prev = NULL;
   
   resnum = pdb->resnum;
   insert = pdb->insert[0];
   strcpy(chain,pdb->chain);

   for(p=pdb, count=1; p!=NULL; NEXT(p))
   {
      if((p->resnum    != resnum) ||
         !CHAINMATCH(p->chain,chain) ||
         (p->insert[0] != insert))
      {
         if(++count > length)
         {
            prev->next = NULL;
            return(p);
         }

         resnum = p->resnum;
         insert = p->insert[0];
         strcpy(chain,p->chain);
      }
      prev = p;
   }

   return(NULL);
}



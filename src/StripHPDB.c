/************************************************************************/
/**

   \file       StripHPDB.c
   
   \version    V1.7
   \date       16.10.96
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
-  V1.0  01.03.90 Original   By: ACRM
-  V1.1  28.03.90 Modified to match new version of pdb.h
-  V1.2  24.05.90 Fixed so the variables passed in as sel[] don't 
                  *have* to be 4 chars.
-  V1.3  17.05.93 Modified for book. Returns BOOL.
-  V1.4  09.07.93 Modified to return PDB pointer. Changed allocation 
                  scheme. Changed back to sel[] variables *must* be 4
                  chars.
-  V1.5  01.11.94 Added HStripPDB()
-  V1.6  26.07.95 Removed unused variables
-  V1.7  16.10.96 Added SelectCaPDB()

*************************************************************************/
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
/*>PDB *StripHPDB(PDB *pdbin, int *natom)
   --------------------------------------
*//**

   \param[in]     *pdbin      Input list
   \param[out]    *natom      Number of atoms kept
   \return                    Output list

   Take a PDB linked list and returns the PDB list minus hydrogens

   N.B. The routine is non-destructive; i.e. the original PDB linked 
        list is intact after the selection process

-  01.11.94 Original based on SelAtomsPDB()   By: ACRM
-  26.07.95 Removed unused variables
*/
PDB *StripHPDB(PDB *pdbin, int *natom)
{
   PDB   *pdbout  = NULL,
         *p,
         *q;
    
   *natom = 0;
   
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!= NULL; NEXT(p))
   {
      if(p->atnam[0] != 'H')
      {
         /* Allocate a new entry                                        */
         if(pdbout==NULL)
         {
            INIT(pdbout, PDB);
            q = pdbout;
         }
         else
         {
            ALLOCNEXT(q, PDB);
         }
         
         /* If failed, free anything allocated and return               */
         if(q==NULL)
         {
            if(pdbout != NULL) FREELIST(pdbout,PDB);
            *natom = 0;
            return(NULL);
         }
         
         /* Increment atom count                                        */
         (*natom)++;
         
         /* Copy the record to the output list (sets ->next to NULL)    */
         CopyPDB(q, p);
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}



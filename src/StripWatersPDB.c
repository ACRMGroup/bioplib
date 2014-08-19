/************************************************************************/
/**

   \file       StripWatersPDB.c
   
   \version    V1.2
   \date       19.08.14
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2008-2014
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
-  V1.0  30.04.08 Original based on StripHPDB()   By: ACRM
-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  19.08.14 Renamed function blStripWatersPDB() to 
                  blStripWatersPDBAsCopy() By: CTP

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
/*>PDB *blStripWatersPDBAsCopy(PDB *pdbin, int *natom)
   ---------------------------------------------------
*//**

   \param[in]     *pdbin      Input list
   \param[out]    *natom      Number of atoms kept
   \return                    Output list

   Take a PDB linked list and returns the PDB list minus waters

   N.B. The routine is non-destructive; i.e. the original PDB linked 
        list is intact after the selection process

-  30.04.08 Original based on StripHPDB()   By: ACRM
-  19.08.14 Renamed function to blStripWatersPDBAsCopy() By: CTP
*/
PDB *blStripWatersPDBAsCopy(PDB *pdbin, int *natom)
{
   PDB   *pdbout  = NULL,
         *p,
         *q;
    
   *natom = 0;
   
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!= NULL; NEXT(p))
   {
      if(!ISWATER(p))
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
         blCopyPDB(q, p);
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}

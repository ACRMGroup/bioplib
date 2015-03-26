/************************************************************************/
/**

   \file       GetPDBChainAsCopy.c
   
   \version    V1.0
   \date       26.03.15
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2015
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
-  V1.0  26.03.16 Original By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Extracting data
   #FUNCTION  blGetPDBChainAsCopy()
   Extracts a specified chain from a PDB linked list allocating a new
   list containing only that chain. The original list is unchanged.
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
/*>PDB *blGetPDBChainAsCopy(PDB *pdbin, char *chain)
   -------------------------------------------------
*//**
   \param[in]    *pdbin     PDB linked list
   \param[in]    *chain     Chain label
   \return                  PDB linked list for requested chain

   Extracts a specified chain from a PDB linked list allocating a new
   list containing only that chain. The original list is unchanged.

-  26.03.15  Original   By: ACRM
*/
PDB *blGetPDBChainAsCopy(PDB *pdbin, char *chain)
{
   PDB *pdbout  = NULL,
       *p,
       *q;
    
   /* Step through the input PDB linked list                            */
   for(p=pdbin; p!=NULL; NEXT(p))
   {
      if(CHAINMATCH(p->chain, chain))
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
            FREELIST(pdbout,PDB);
            return(NULL);
         }
        
         /* Copy the record to the output list (sets ->next to NULL)    */
         blCopyPDB(q, p);
      }
   }

   /* Return pointer to start of output list                            */
   return(pdbout);
}



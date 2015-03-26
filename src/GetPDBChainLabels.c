/************************************************************************/
/**

   \file       GetPDBChainLabels.c
   
   \version    V1.13
   \date       26.03.15
   \brief      
   
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
-  V1.11 25.03.14 Deprecated GetPDBChainLabels() and added replacement 
                  function blGetPDBChainLabels() By: CTP
-  V1.12 31.07.14 Moved GetPDBChainLabels() to dperecated.c By: CTP
-  V1.13 26.03.15 blGetPDBChainLabels() no longer assumes all instances
                  of one label occur together  By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Extracting data
   #FUNCTION  blGetPDBChainLabels()
   Scans a PDB linked list for chain names. Allocates memory for an 
   array of strings containing these labels which is returned together 
   with the number of chains found.
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
/*>char **blGetPDBChainLabels(PDB *pdb, int *nChains)
   --------------------------------------------------
*//**

   \param[in]     *pdb        PDB linked list
   \param[out]    nChains     Number of chains found.
   \return                    Allocated array of strings containing chain
                              labels. NULL if unable to allocate memory.

   Scans a PDB linked list for chain names. Allocates memory for an 
   array of strings containing these labels which is returned together 
   with the number of chains found.

   N.B. You must free the allocated memory for both the array of chains 
        and for each individual chain label when you've finished with it!

-  25.03.14 Original based on GetPDBChainLabels(). By: CTP
-  26.03.15 No longer assumes chain labels do not appear again later
            (i.e. chains L,H,Y,L only returns L,H,Y)  By: ACRM
*/
char **blGetPDBChainLabels(PDB *pdb, int *nChains)
{
   char **chains = NULL,
        prevChain[8];
   PDB  *p;

   /* Initialize previous chain and zero chain count                    */
   prevChain[0] = '\0';
   *nChains = 0;
   
   /* Just return if linked list is NULL                                */
   if(pdb == NULL)
      return(NULL);
   
   /* Allocate memory for first chain                                   */
   if((chains = (char **)malloc(sizeof(char *)))==NULL)
      return(NULL);
   if((chains[0] = (char *)malloc(8 * sizeof(char)))==NULL)
      return(NULL);
   *nChains = 1;      
   
   /* Copy in the first chain                                           */
   strncpy(chains[0], pdb->chain, 8);  
   strncpy(prevChain, pdb->chain, 8);

   /* Run through the pdb linked list                                   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(!CHAINMATCH(p->chain, prevChain))
      {
         BOOL found = FALSE;
         int  chainNum;
         
         strncpy(prevChain, p->chain, 8);

         /* See if this chain has appeared before                       */
         for(chainNum=0; chainNum<(*nChains); chainNum++)
         {
            if(CHAINMATCH(p->chain, chains[chainNum]))
            {
               found = TRUE;
               break;
            }
         }
         /* If it hasn't appeared before, allocate memory to store it
            and copy in the new label
         */
         if(!found)
         {
            if((chains = 
                (char **)realloc(chains, (*nChains + 1)*sizeof(char *)))
               == NULL)
               return(NULL);

            if((chains[*nChains] = 
                (char *)malloc(8 * sizeof(char))) == NULL)
               return(NULL);

            strncpy(chains[*nChains],p->chain, 8);  
            (*nChains)++;
         }
      }
   }

   return(chains);
}

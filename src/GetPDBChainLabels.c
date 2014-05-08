/*************************************************************************

   Program:    
   File:       GetPDBChainLabels.c
   
   Version:    V1.11
   Date:       25.03.14
   Function:   
   
   Copyright:  (c) SciTech Software 1992-6
   Author:     Dr. Andrew C. R. Martin
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@bioinf.org.uk
               
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
   V1.11 25.03.14 Deprecated GetPDBChainLabels() and added replacement 
                  function blGetPDBChainLabels() By: CTP

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
/*>char *GetPDBChainLabels(PDB *pdb)
   ---------------------------------
   Input:   PDB    *pdb      PDB linked list
   Returns: char   *         Allocated string containing chain labels
                             NULL if unable to allocate memory

   Scans a PDB linked list for chain names. Allocates memory for a 
   string containing these labels which is returned.

   N.B. You must free the allocated memory when you've finished with it!

   25.07.95 Original    By: ACRM
   25.03.14 Added deprecated message. By: CTP
   07.05.14 Use DEPRECATED() macro. By: CTP
*/
char *GetPDBChainLabels(PDB *pdb)
{
   char *chains;
   int  nchains   = 0,
        maxchains = 16;
   PDB  *p;

   DEPRECATED("GetPDBChainLabels()","blGetPDBChainLabels()");
   
   /* Just return if linked list is NULL                                */
   if(pdb==NULL)
      return(NULL);

   /* Allocate a chunk for storing the chains                           */
   if((chains = (char *)malloc(maxchains * sizeof(char)))==NULL)
      return(NULL);

   /* Set up first chain label                                          */
   chains[nchains] = pdb->chain[0];

   /* Run through the linked list                                       */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(p->chain[0] != chains[nchains])
      {
         /* Increment chain count and reallocate memory if needed       */
         if(++nchains == maxchains)
         {
            maxchains += 16;
            if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
               return(NULL);
         }
         /* Store this new chain label                                  */
         chains[nchains] = p->chain[0];
      }
   }

   /* Increment chain count and reallocate memory if needed             */
   if(++nchains == maxchains)
   {
      maxchains += 16;
      if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
         return(NULL);
   }

   /* Terminate the chain list with a NUL character                     */
   chains[nchains] = '\0';

   return(chains);
}


/************************************************************************/
/*>char **blGetPDBChainLabels(PDB *pdb, int *nchains)
   --------------------------------------------------
   Input:   PDB    *pdb      PDB linked list
   Returns: int    *nchains  Number of chains found.
            char   **        Allocated array of strings containing chain
                             labels. NULL if unable to allocate memory.

   Scans a PDB linked list for chain names. Allocates memory for an 
   array of strings containing these labels which is returned together 
   with the number of chains found.

   N.B. You must free the allocated memory for both the array of chains 
        and for each individual chain label when you've finished with it!

   25.03.14 Original based on GetPDBChainLabels(). By: CTP
   
*/
char **blGetPDBChainLabels(PDB *pdb, int *nchains)
{
   char **chains;
   PDB  *p;
   
   /* Zero Chain Count                                                  */
   *nchains = 0;
   
   /* Just return if linked list is NULL                                */
   if(pdb == NULL)
   {
      return(NULL);
   }
   
   /* Get First Chain                                                   */
   /* allocate memory */
   if((chains = (char **)malloc(sizeof(char *)))==NULL)
   {
      return(NULL);
   }
   if((chains[0] = (char *)malloc(8 * sizeof(char)))==NULL)
   {
      return(NULL);
   }
   
   /* get first chain and set count */
   strcpy(chains[0],pdb->chain);  
   *nchains = 1;      

   /* Run through the pdb linked list                                   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(!CHAINMATCH(chains[*nchains - 1],p->chain))
      {
         /* allocate/reallocate memory */
         if( (chains = (char **)realloc( chains, 
                                         (*nchains + 1) * sizeof(char *) )
             ) == NULL )
         {
            return(NULL);
         }
         if((chains[*nchains] = (char *)malloc(8 * sizeof(char))) == NULL)
         {
            return(NULL);
         }

         /* get chain and increment count */
         strcpy(chains[*nchains],p->chain);  
         *nchains += 1;
      }
   }

   return(chains);
}

/************************************************************************/
/**

   \file       GetPDBChainLabels.c
   
   \version    V1.12
   \date       31.07.14
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
-  V1.11 25.03.14 Deprecated GetPDBChainLabels() and added replacement 
                  function blGetPDBChainLabels() By: CTP
-  V1.12 31.07.14 Moved GetPDBChainLabels() to dperecated.c By: CTP

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
/*>char **blGetPDBChainLabels(PDB *pdb, int *nchains)
   --------------------------------------------------
*//**

   \param[in]     *pdb        PDB linked list
   \param[out]    nchains     Number of chains found.
   \return                    Allocated array of strings containing chain
                              labels. NULL if unable to allocate memory.

   Scans a PDB linked list for chain names. Allocates memory for an 
   array of strings containing these labels which is returned together 
   with the number of chains found.

   N.B. You must free the allocated memory for both the array of chains 
        and for each individual chain label when you've finished with it!

-  25.03.14 Original based on GetPDBChainLabels(). By: CTP
   
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

/************************************************************************/
/**

   \file       SplitSeq.c
   
   \version    V1.11
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2000
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
-  V1.0  29.09.92 Original
-  V1.1  07.06.93 Corrected allocation
-  V1.2  18.06.93 Handles multi-chains and skips NTER and CTER residues.
                  Added SplitSeq()
-  V1.3  09.07.93 SplitSeq() cleans up properly if allocation failed
-  V1.4  11.05.94 Added TrueSeqLen()
-  V1.5  13.05.94 Fixed bug in PDB2Seq().
                  Added KnownSeqLen().
-  V1.6  07.09.94 Fixed allocation bug in SplitSeq()
-  V1.7  19.07.95 Added check for ATOM records
-  V1.8  24.01.96 Fixed bug when no ATOM records in linked list
                  Returns a blank string
-  V1.9  26.08.97 Renamed DoPDB2Seq() with handling of Asx/Glx and
                  protein-only. Added macros to recreate the
                  old PDB2Seq() interface and similar new calls
-  V1.10 02.10.00 Added NoX option
-  V1.11 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP blSplitSeq()
   Splits a sequence stored as a linear array with each chain separated
   by a * into an array of sequences. Returns the number of chains
   found.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

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
/*>int blSplitSeq(char *LinearSeq, char **seqs)
   --------------------------------------------
*//**

   \param[in]     *LinearSeq   Array containing sequence with chains
                               terminated by *'s
   \param[out]    **seqs       Allocated set of character arrays 
                               containing one chain per array
   \return                        Number of chains found

   Splits a sequence stored as a linear array with each chain separated
   by a * into an array of sequences. Returns the number of chains
   found.
   
-  18.06.93 Original    By: ACRM
-  09.07.93 Cleans up properly of allocation failed
-  07.09.94 Sequence space was being allocated one too small
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blSplitSeq(char *LinearSeq, char **seqs)
{
   char  *ptr,
         *star;
   int   seqlen   = strlen(LinearSeq),
         NSeq     = 0;
   
   ptr = LinearSeq;
   
   while(ptr-LinearSeq < seqlen)
   {
      star = strchr(ptr,'*');
      if(star != NULL)
      {
         *star = '\0';
         seqs[NSeq] = (char *)malloc((1+strlen(ptr)) * sizeof(char));
         if(seqs[NSeq] == NULL)
         {
            int i;
            for(i=0; i<=NSeq; i++) if(seqs[i] != NULL) free(seqs[i]);
            return(0);
         }
         strcpy(seqs[NSeq],ptr);
         NSeq++;
         ptr = star+1;
      }
      else
      {
         if(ptr < LinearSeq+seqlen)
         {
            seqs[NSeq] = (char *)malloc((1+strlen(ptr)) * sizeof(char));
            if(seqs[NSeq] == NULL)
            {
               int i;
               for(i=0; i<=NSeq; i++) if(seqs[i] != NULL) free(seqs[i]);
               return(0);
            }
            strcpy(seqs[NSeq],ptr);
            NSeq++;
         }
         break;
      }
   }
   
   return(NSeq);
}


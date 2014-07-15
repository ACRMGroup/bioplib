/************************************************************************/
/**

   \file       TrueSeqLen.c
   
   \version    V1.11
   \date       07.07.14
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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
/* Includes
*/
#include "deprecated.h"

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
/*>int blTrueSeqLen(char *sequence)
   --------------------------------
*//**

   \param[in]     *sequence    A sequence containing deletions
   \return                        Length without deletions

   Scans a 1-letter code sequence and calculate the length without
   `-' or ` ' residues

-  14.04.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP

*/
int blTrueSeqLen(char *sequence)
{
   int length = 0,
       i = 0;
   
   for(i=0; sequence[i]; i++)
   {
      if(sequence[i] != '-' && sequence[i] != ' ')
         length++;
   }
   
   return(length);
}


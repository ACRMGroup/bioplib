/************************************************************************/
/**

   \file       WritePIR.c
   
   \version    V1.0
   \date       11.06.16
   \brief      Writes a PIR sequence file
   
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
-  V1.0  11.06.15 

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP File IO
   PIR and FASTA writer
*/
/************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "SysDefs.h"
#include "macros.h"
#include "seq.h"

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
/*>void blWriteOneStringPIR(FILE *out, char *label, char *title, 
                            char *sequence, char **chainLabels, 
                            BOOL ByChain, BOOL doFasta)
   --------------------------------------------------------------------
*//**

   \param[in]  *out      File pointer
   \param[in]  *label    Sequence label
   \param[in]  *title    Sequence title
   \param[in]  *sequence Sequence (1-letter code) with chains separated 
                         by *
   \param[in]  **chainLabels  Chain labels (may be set to NULL unless
                         ByChain is set)
   \param[in]  ByChain   Print a separate header for each chain
   \param[in]  doFasta   Output FASTA format instead of PIR

   Writes a PIR sequence file from a 1-letter code sequence. Multiple
   chains are split with '*'. If ByChain is set the the chainLabels array
   must be non-NULL and contains labels for each chain
   Adds a terminating * if required.

-  10.05.94 Original    By: ACRM
-  22.08.97 Can now handle chains separately
-  26.08.97 If chains are handled separately, don't bother writing out
            an empty chain
-  10.08.98 Basically a total rewrite to fix a bug which caused the 
            header not to be printed with -c -p for a chain after one
            which was non-protein. Much simplified the code by printing
            the header at the beginning of a chain rather than end
            of previous chain.
-  18.10.00 Added code to write FASTA as well
-  11.06.15 Reset count=0 after a * - tidies up the output!
*/
void blWriteOneStringPIR(FILE *out, char *label, char *title, 
                         char *sequence, char **chainLabels, 
                         BOOL ByChain, BOOL doFasta)
{
   int  i, 
        count,
        chaincount = 0;
   BOOL GotStar = FALSE,
        DoneHeader = FALSE,
        Printed;
   char chn[blMAXCHAINLABEL],
        outlabel[blMAXPIRLABEL];

   strcpy(outlabel,label);

   /* If we are not going by chain create a label                       */
   if(!ByChain)
   {
      strcpy(outlabel, "PDBPIR");
      if(label[0])
         strcpy(outlabel, label);
   }

   /* If we are not doing it by-chain, then simple display the label
      and header now otherwise set flag to say we haven't done header
   */
   if(ByChain)
   {
      DoneHeader = FALSE;
   }
   else
   {
      if(doFasta)
      {
         fprintf(out,">%s\n",outlabel);
      }
      else
      {
         fprintf(out,">P1;%s\n",outlabel);
         fprintf(out,"Sequence extracted from PDB file - %s\n",
                 (title[0]?title:"By pdb2pir"));
      }
      DoneHeader = TRUE;
   }

   /* Start by setting flag to say nothing has been printed             */
   Printed = FALSE;
   
   /* Loop through the sequence with i. Use count to count number of
      residues printed on a line

      11.06.15 Fixed initialization of count to zero rather than 1
   */
   for(i=0,count=0; sequence[i]; i++)
   {
      if(ByChain)
      {
         /* If we don't have a header, then print one                   */
         if((sequence[i] != '*') && (!DoneHeader))
         {
            if(label[0])
            {
               strcpy(outlabel,label);
               strcpy(chn, chainLabels[chaincount]);
               strcat(outlabel,chn);
            }
            else
            {
               sprintf(outlabel,"Chain%s",chainLabels[chaincount]);
            }
            if(doFasta)
            {
               fprintf(out,">%s\n",outlabel);
            }
            else
            {
               fprintf(out,">P1;%s\n",outlabel);
               fprintf(out,"Sequence extracted from PDB file - %s\n",
                       (title[0]?title:"By pdb2pir"));
            }
            
            DoneHeader = TRUE;
            count = 0;
         }
      }
      
      if(count==30)
      {
         fputc('\n',out);
         count = 0;
      }

      if(sequence[i] == '*')
      {
         chaincount++;
         if(Printed)
         {
            if(doFasta)
            {
               fprintf(out,"\n");
            }
            else
            {
               fprintf(out,"*\n");
            }
            
            GotStar = TRUE;
         }
         
         DoneHeader = FALSE;
         Printed    = FALSE;
         count      = 0;          /* 11.06.15                           */
      }
      else
      {
         fputc(sequence[i],out);
         GotStar = FALSE;
         Printed = TRUE;
      }

      if(Printed)
         count++;
   }

   if(!GotStar)
   {
      if(doFasta)
      {
         fprintf(out,"\n");
      }
      else
      {
         fprintf(out,"*\n");
      }
   }
}



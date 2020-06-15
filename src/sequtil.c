/************************************************************************/
/**

   \file       sequtil.c
   
   \version    V1.1
   \date       10.11.17   
   \brief      Six-frame translation and FASTA handling
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2017-2019
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
-  V1.0  10.11.17 Original  By: ACRM
-  V1.1  01.11.19 Moved MAXBUFF into here from sequtil.h

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP Sequence Utilities

   #FUNCTION blSixFTBest()
   Performs a 6FT and find the longest translation starting from the
   beginning of the DNA or a Met. Returns malloc'd memory for the 
   protein sequence. Optionally outputs the ORF that was translated.

   #FUNCTION blReverseComplement()
   Malloc's a reverse complement sequence

   #FUNCTION blTranslateFrame()
   Translates a given frame of DNA filling in the protein sequence for
   the complete DNA - stop codons are indicated with *

   #FUNCTION blFindLongestTranslation()
   Finds the longest protein sequence within the long sequence where
   each section is separated with a *

   #FUNCTION blWriteFASTA()
   Writes a sequence in FASTA format

   #FUNCTION blReadFASTA()
   Each call reads another sequence from a FASTA file allocating
   memory for the sequence. Note that this is not re-entrant - it 
   uses a static buffer to look ahead to the next line.

   #FUNCTION blReadFASTAExtBuffer()
   Each call reads another sequence from a FASTA file allocating
   memory for the sequence.

*/

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "sequtil.h"
#include "macros.h"
#include "general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256

/************************************************************************/
/* Globals
*/

/************************************************************************/
/*>char *blSixFTBest(char *inDna, char *orf)
   -----------------------------------------
*//**
   \param[in]    inDna    DNA sequence
   \param[out]   orf      DNA of the ORF (if NULL, doesn't provide this)
   \return                Malloc'd protein sequence

   Performs a 6FT and find the longest translation starting from the
   beginning of the DNA or a Met. Returns malloc'd memory for the 
   protein sequence. Optionally outputs the ORF that was translated.

-  10.11.17 Original   By: ACRM
*/
char *blSixFTBest(char *inDna, char *orf)
{
   char *bestProtSeq = NULL,
        *protSeq     = NULL,
        *dna[2];
   int  bestLen      = 0,
        i, frame;

   /* For working copy of protein sequence                              */
   if((protSeq = (char *)malloc((2+strlen(inDna)/3) * sizeof(char)))
      == NULL)
      return(NULL);

   /* For best protein sequence used as output                          */
   if((bestProtSeq = (char *)malloc((2+strlen(inDna)/3) * sizeof(char)))
      == NULL)
   {
      free(protSeq);
      return(NULL);
   }

   /* Put the DNA and RC DNA into an array                              */
   dna[0] = inDna;
   dna[1] = blReverseComplement(inDna);

   if(dna[1] == NULL)
   {
      free(protSeq);
      free(bestProtSeq);
      return(NULL);
   }
   
   /* For the DNA and RC DNA                                            */
   for(i=0; i<2; i++)
   {
      /* For each reading frame                                         */
      for(frame=0; frame<3; frame++)
      {
         int offset,
             protLen;

         blTranslateFrame(dna[i], frame, protSeq);
         offset = blFindLongestTranslation(protSeq, &protLen);
         if(protLen > bestLen)
         {
            bestLen = protLen;
            strncpy(bestProtSeq, protSeq+offset, protLen);
            bestProtSeq[protLen] = '\0';
            
            if(orf != NULL)
            {
               strncpy(orf, dna[i]+(offset*3)+frame, protLen*3);
               orf[protLen*3] = '\0';
            }
            
         }
      }
   }

   FREE(dna[1]);
   
   return(bestProtSeq);
}


/************************************************************************/
/*>char *blReverseComplement(char *dna)
   ------------------------------------
*//**
   \param[in]    dna   DNA sequence
   \return             Malloc'd reverse complement DNA sequence

   Malloc's a reverse complement sequence

-  10.11.17 Original   By: ACRM
*/
char *blReverseComplement(char *dna)
{
   char *rcDNA = NULL;
   int  i, dnaLen;
   
   dnaLen = strlen(dna);

   if((rcDNA = (char *)malloc((1+dnaLen)*sizeof(char)))==NULL)
      return(NULL);

   for(i=0; i<dnaLen; i++)
   {
      switch(dna[dnaLen-1-i])
      {
      case 'a':
      case 'A':
         rcDNA[i] = 't';
         break;
      case 't':
      case 'u':
      case 'T':
         rcDNA[i] = 'a';
         break;
      case 'c':
      case 'C':
         rcDNA[i] = 'g';
         break;
      case 'g':
      case 'G':
         rcDNA[i] = 'c';
         break;
      default:
         rcDNA[i] = 'n';
      }
   }
   rcDNA[dnaLen] = '\0';
   
   return(rcDNA);
}


/************************************************************************/
/*>void blTranslateFrame(char *dna, int frame, char *protein)
   ----------------------------------------------------------
*//**
   \param[in]    dna      DNA sequence
   \param[in]    frame    Reading frame (0..2)
   \param[out]   protein  Protein sequence with * for stop codons

   Translates a given frame of DNA filling in the protein sequence for
   the complete DNA - stop codons are indicated with *

-  10.11.17 Original   By: ACRM
*/
void blTranslateFrame(char *dna, int frame, char *protein)
{
   int i, j, k, dnaLen;
   
   static char *codons[] = {"aaa", "aac", "aag", "aat",
                            "aca", "acc", "acg", "act",
                            "aga", "agc", "agg", "agt",
                            "ata", "atc", "atg", "att",

                            "caa", "cac", "cag", "cat",
                            "cca", "ccc", "ccg", "cct",
                            "cga", "cgc", "cgg", "cgt",
                            "cta", "ctc", "ctg", "ctt",

                            "gaa", "gac", "gag", "gat",
                            "gca", "gcc", "gcg", "gct",
                            "gga", "ggc", "ggg", "ggt",
                            "gta", "gtc", "gtg", "gtt",
                          
                            "taa", "tac", "tag", "tat",
                            "tca", "tcc", "tcg", "tct",
                            "tga", "tgc", "tgg", "tgt",
                            "tta", "ttc", "ttg", "ttt", NULL };
   
   static char aas[] = {'K', 'N', 'K', 'N',
                        'T', 'T', 'T', 'T',
                        'R', 'S', 'R', 'S',
                        'I', 'I', 'M', 'I',

                        'Q', 'H', 'Q', 'H',
                        'P', 'P', 'P', 'P',
                        'R', 'R', 'R', 'R',
                        'L', 'L', 'L', 'L',

                        'E', 'D', 'E', 'D',
                        'A', 'A', 'A', 'A',
                        'G', 'G', 'G', 'G',
                        'V', 'V', 'V', 'V',
                 
                        '*', 'Y', '*', 'Y',
                        'S', 'S', 'S', 'S',
                        '*', 'C', 'W', 'C',
                        'L', 'F', 'L', 'F', '\0' };
   j=0;k=0;
   dnaLen = strlen(dna);

   for(i=frame; i<dnaLen; i+=3)
   {
      protein[k] = 'X';
      
      if(strlen(dna+i) >= 3)
      {
         for(j=0; codons[j] != NULL; j++)
         {
            if(!strncmp(dna+i, codons[j], 3))
            {
               protein[k] = aas[j];
               break;
            }
         }
         k++;
      }
      
   }

   protein[k] = '\0';
}


/************************************************************************/
/*>int blFindLongestTranslation(char *protSeq, int *protLen)
   -------------------------------------------------------
*//**
   \param[in]    protSeq    Protein sequence with stops indicated
   \param[out]   protLen    Length of the longest section of protein
                            starting from the start of the sequence or
                            from a Met
   \return                  Offset to the start of the longest sequence

   Finds the longest protein sequence within the long sequence where
   each section is separated with a *

-  10.11.17 Original   By: ACRM
*/
int blFindLongestTranslation(char *protSeq, int *protLen)
{
   int seqLen, offset, pLen, bestLen = 0, bestOffset = 0;
   seqLen = strlen(protSeq);
   *protLen = 0;
   
   for(offset=0; offset<seqLen; offset++)
   {
      if((offset==0) || (protSeq[offset] == 'M'))
      {
         int i;
         pLen = 0;
         for(i=offset; (i<seqLen) && (protSeq[i] != '*'); i++)
         {
            pLen++;
         }
         if(pLen > bestLen)
         {
            bestLen = pLen;
            bestOffset = offset;
            *protLen = pLen;
         }
      }
   }

   return(bestOffset);
}


/************************************************************************/
/*>void blWriteFASTA(FILE *out, char *header, char *sequence, 
                int width, BOOL pad)
   ---------------------------------------------------------
*//**
   \param[in]   out       Output file pointer
   \param[in]   header    FASTA header
   \param[in]   sequence  Sequence to print
   \param[in]   width     Width of each line
   \param[in]   pad       Print a space every 10 characters

   Writes a sequence in FASTA format

-  10.11.17 Original   By: ACRM
*/
void blWriteFASTA(FILE *out, char *header, char *sequence, 
                  int width, BOOL pad)
{
   int seqLen, i, padPos;
   
   fprintf(out, "%s\n", header);
   seqLen = strlen(sequence);

   padPos = 0;
   for(i=0; i<seqLen; i++)
   {
      if(((i%width)==0) && (i!=0))
      {
         fputc('\n', out);
         padPos = 0;
      }
      else if(pad && (padPos!=0) && ((padPos%10)==0))
      {
         fputc(' ', out);
      }
      fputc(sequence[i], out);
      padPos++;
   }
   fputc('\n', out);
}


/************************************************************************/
/*>char *blReadFASTA(FILE *in, char *header, int size)
   ---------------------------------------------------
*//**
   \param[in]    in     Input file pointer
   \param[out]   header The FASTA header
   \param[in]    size   The maximum size of the FASTA header  
   \return              Pointer to allocated memory for sequence

   Each call reads another sequence from a FASTA file allocating
   memory for the sequence. Note that this is not re-entrant - it 
   uses a static buffer to look ahead to the next line.

-  10.11.17 Original   By: ACRM
*/
char *blReadFASTA(FILE *in, char *header, int size)
{
   static char buffer[MAXBUFF];
   
   return(blReadFASTAExtBuffer(in, header, size, buffer, MAXBUFF));
}


/************************************************************************/
/*>char *blReadFASTAExtBuffer(FILE *in, char *header, int headerSize, 
                              char *buffer, int bufferSize)
   ------------------------------------------------------------------
*//**
   \param[in]     in         Input file pointer
   \param[out]    header     The FASTA header
   \param[in]     size       The maximum size of the FASTA header  
   \param[in,out] buffer     Buffer to look ahead to next line
   \param[in]     bufferSize Size of the lookahead buffer
   \return                   Pointer to allocated memory for sequence

   Each call reads another sequence from a FASTA file allocating
   memory for the sequence.

-  10.11.17 Original   By: ACRM
*/
char *blReadFASTAExtBuffer(FILE *in, char *header, int headerSize, 
                           char *buffer, int bufferSize)
{
   char *sequence  = NULL,
        *seqBuffer = NULL;

   /* Copy the existing buffer - which should be the next header        */
   strncpy(header, buffer, headerSize);

   while(fgets(buffer, bufferSize, in))
   {
      TERMINATE(buffer);
      
      if(buffer[0] == '>')  /* A header                                 */
      {
         if(seqBuffer == NULL) /* Only on first sequence                */
         {
            strncpy(header, buffer, headerSize);
         }
         else
         {
            break;
         }
      }
      else
      {
         if((seqBuffer = blStrcatalloc(seqBuffer, buffer)) == NULL)
            return(NULL);
      }
   }
   
   sequence = blRemoveSpaces(seqBuffer);
   free(seqBuffer);

   return(sequence);
}



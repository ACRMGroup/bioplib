/************************************************************************/
/**

   Program:    
   \file       sequtil.h
   
   \version    V1.0
   \date       10.11.17   
   \brief      Sequence utilities
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2017
   \author     Dr. Andrew C. R. Martin
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
   - V1.0   10.11.17  Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256

/************************************************************************/
/* Globals
*/

char *blSixFTBest(char *dna, char *orf);
char *blReverseComplement(char *dna);
int blFindLongestTranslation(char *protSeq, int *protLen);
void blWriteFASTA(FILE *out, char *header, char *sequence, 
                int width, BOOL pad);
char *blReadFASTA(FILE *in, char *header, int size);
char *blReadFASTAExtBuffer(FILE *in, char *header, int headerSize, 
                         char *buffer, int bufferSize);
char *blRemoveSpaces(char *inText);
void blTranslateFrame(char *dna, int frame, char *protein);

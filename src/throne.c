/*************************************************************************

   Program:    
   File:       throne.c
   
   Version:    V1.2R
   Date:       25.07.95
   Function:   Convert between 1 and 3 letter aa codes
   
   Copyright:  (c) SciTech Software 1993-5
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
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
   V1.0  29.09.92 Original    By: ACRM   
   V1.1  11.03.94 Added PCA, ASX and GLX to translation table.
                  PCA translates to E
                  Added routines to handle asx/glx
   V1.2  25.07.95 handles nucleic acids
                  Sets the gBioplibSeqNucleicAcid flag if it's a 
                  nucleic acid.

*************************************************************************/
/* Includes
*/
#include <string.h>
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define NUMAAKNOWN 29

/************************************************************************/
/* Globals
*/

/* N.B. The order in sTab1[] and sTab3[] must be the same and they must
   end with X/UNK.
   Also, nucleic acids must come *after* amino acids.
*/
static char sTab1[]    = {'A','C','D','E','F',
                          'G','H','I','K','L',
                          'M','N','P','Q','R',
                          'S','T','V','W','Y',
                          'E','B','Z',
                          'A','T','C','G','U',
                          'X'
                         };
static char sTab3[][8] = {"ALA ","CYS ","ASP ","GLU ","PHE ",
                          "GLY ","HIS ","ILE ","LYS ","LEU ",
                          "MET ","ASN ","PRO ","GLN ","ARG ",
                          "SER ","THR ","VAL ","TRP ","TYR ",
                          "PCA ","ASX ","GLX ",
                          "  A ","  T ","  C ","  G ","  U ",
                          "UNK "
                         };

BOOL gBioplibSeqNucleicAcid = FALSE;

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>char throne(char *three)
   ------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as X
   
   29.09.92 Original    By: ACRM
   11.03.94 Modified to handle ASX and GLX in the tables
   25.07.95 Added handling of gBioplibSeqNucleicAcid
*/
char throne(char *three)
{
   int j;

   if(three[0] == ' ' && three[1] == ' ')
      gBioplibSeqNucleicAcid = TRUE;
   else
      gBioplibSeqNucleicAcid = FALSE;

   if(three[2] == 'X')
      return('X');

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}


/************************************************************************/
/*>char thronex(char *three)
   -------------------------
   Input:   char  *three    Three letter code
   Returns: char            One letter code

   Converts 3-letter code to 1-letter code.
   Handles ASX and GLX as B and Z.
   
   29.09.92 Original    By: ACRM
   25.07.95 Added handling of gBioplibSeqNucleicAcid
*/
char thronex(char *three)
{
   int j;

   if(three[0] == ' ' && three[1] == ' ')
      gBioplibSeqNucleicAcid = TRUE;
   else
      gBioplibSeqNucleicAcid = FALSE;

   for(j=0;j<NUMAAKNOWN;j++)
      if(!strncmp(sTab3[j],three,3)) return(sTab1[j]);

   /* Only get here if the three letter code was not found              */
   return('X');
}


/************************************************************************/
/*>char *onethr(char one)
   ----------------------
   Input:   char  one     One letter code
   Returns: char  *       Three letter code (padded to 4 chars with a 
                          space)

   Converts 1-letter code to 3-letter code (actually as 4 chars).

   07.06.93 Original    By: ACRM
   25.07.95 If the gBioplibSeqNucleicAcid flag is set, assumes nucleic
            acids rather than amino acids
*/
char *onethr(char one)
{
   int j;

   if(gBioplibSeqNucleicAcid) /* Work from end of table                 */
   {
      for(j=NUMAAKNOWN-1;j>=0;j++)
         if(sTab1[j] == one) return(sTab3[j]);
   }
   else                       /* Work from start of table               */
   {
      for(j=0;j<NUMAAKNOWN;j++)
         if(sTab1[j] == one) return(sTab3[j]);
   }

   /* Only get here if the one letter code was not found                */
   return(sTab3[NUMAAKNOWN-1]);
}


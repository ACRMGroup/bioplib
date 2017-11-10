/************************************************************************/
/**

   \file       stringutil.c
   
   \version    V1.1
   \date       10.11.17
   \brief      String utilities
   
   \copyright  (c) UCL / Dr. Andrew C.R. Martin, 2015-2017
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

   See documentation for details

**************************************************************************

   Revision History:
   =================
-  V1.0  28.04.15 Original
-  V1.1  10.11.17 Added blRemoveSpaces()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP String handling

   #FUNCTION blCollapseSpaces()
   Takes a string and collapses multiple spaces down to a single space
   Equivalent of perl 's/\s+/ /g'

   #FUNCTION blStrdup()
   Duplicates a string, allocating memory for it.
   An implementation of SVr4, 4.3BSD, POSIX.1-2001 strdup() which is 
   not standard ANSI C

   #FUNCTION blRemoveSpaces()
   Strips all whitespace out of a string. Allocates a new string.
*/

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXWORD 8

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>char *blCollapseSpaces(char *inText)
   ------------------------------------
*//**
   \param[in]     *inText   Input text string
   \return                  malloc()'d string with multiple spaces 
                            collapsed

   Takes a string and collapses multiple spaces down to a single space
   Equivalent of perl 's/\s+/ /g'
   The input string is unmodified and malloc()s the output.

   28.04.15  Original   By: ACRM
*/
char *blCollapseSpaces(char *inText)
{
   int  nchar = 0;
   char *ch, *chp, *chq,
        *outText = NULL;

   if(inText==NULL)
      return(NULL);
   
   /* Get length of input string                                        */
   nchar = strlen(inText) + 1;

   /* Allocate new space                                                */
   if((outText=(char *)malloc(nchar * sizeof(char)))==NULL)
      return(NULL);

   /* Copy characters skipping repeated spaces                          */
   chp=NULL;
   chq=outText;
   for(ch=inText; *ch!='\0'; ch++)
   {
      if(*ch == '\t') *ch = ' '; /* Tab to space                        */
      if((chp   == NULL) ||      /* 1st char                            */
         (*chp != ' ')   ||      /* Prev char not \s                    */
         (*ch  != ' '))          /* Curr char not \s                    */
      {
         *chq = *ch;
         chq++;
      }
      chp=ch;
   }
   *chq = '\0';

   return(outText);
}


/************************************************************************/
/*>char *blStrdup(char *instr)
   ---------------------------
*//**
   \param[in]   instr    A string
   \return               A malloc'd copy of the string

   An implementation of SVr4, 4.3BSD, POSIX.1-2001 strdup() which is 
   not standard ANSI C

-  12.05.15  Original   By: ACRM
*/
char *blStrdup(char *instr)
{
   int  len     = strlen(instr)+1;
   char *outstr = NULL;

   if((outstr = (char*)malloc(len*sizeof(char)))!=NULL)
      strcpy(outstr, instr);
   
   return(outstr);
}


/************************************************************************/
/*>char *blRemoveSpaces(char *inText)
   --------------------------------
*//**
   \param[in]    inText  Input string
   \return               Malloc'd string without spaces

   Allocates a string and copies the input to it skipping whitespace.

-  10.11.17 Original   By: ACRM
*/
char *blRemoveSpaces(char *inText)
{
   int  nchar = 0;
   char *chIn, 
        *chOut,
        *outText = NULL;

   if(inText==NULL)
      return(NULL);

   /* Count the non-space characters                                    */
   for(chIn=inText; *chIn!='\0'; chIn++)
   {
      if((*chIn != '\t') && 
         (*chIn != '\n')  &&
         (*chIn != '\r')  &&
         (*chIn != ' '))
      {
         nchar++;
      }
   }
   nchar++;
   
   /* Allocate new space                                                */
   if((outText=(char *)malloc(nchar * sizeof(char)))==NULL)
      return(NULL);

   /* Copy characters skipping repeated spaces                          */
   chOut = outText;
   for(chIn=inText; *chIn!='\0'; chIn++)
   {
      if((*chIn != '\t') && 
         (*chIn != '\n') &&
         (*chIn != '\r') &&
         (*chIn != ' '))
      {
         *chOut = *chIn;
         chOut++;
      }
   }
   *chOut = '\0';

   return(outText);
}



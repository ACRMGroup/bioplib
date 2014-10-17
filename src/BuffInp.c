/************************************************************************/
/**

   \file       BuffInp.c
   
   \version    V1.2
   \date       07.07.14
   \brief      Read from a file a line at a time, allowing one to probe
               ahead and look at the contants of the next line without
               removing it from the input stream.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2014
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
-  V1.0  08.03.94 Original
-  V1.1  11.03.94 Changes to ReadBufferedFile() and ProbeBufferedFile()
                  to RePrompt() if reading from stdin when we get a 
                  blank line
-  V1.2  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP File IO
   #ROUTINE  blOpenBufferedFile()
   Open a file for buffered input. This allows probe-ahead to look at the
   contents of the next line without removing it from the input stream.

   #ROUTINE  blReadBufferedFile()
   Reads a line from a buffered file (like fgets()).
   Blank lines in the file will be skipped.

   #ROUTINE  blProbeBufferedFile()
   Read the next line from a buffered file without removing it from
   the input stream. Repeated calls will thus return the same string.
   The next call to ReadBufferedFile will also output the same string,
   but will remove the line from the input stream.
   Blank lines in the file will be skipped.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "BuffInp.h"
#include "WindIO.h"

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
/*>INBUFFER blOpenBufferedFile(char *filename, int maxstr)
   -------------------------------------------------------
*//**

   \param[in]     *filename   File name
   \param[in]     maxstr      Max string length in file for buffering
   \return                    Pointer to a buffered file

   Open a file for buffered input. This allows probe-ahead to look at the
   contents of the next line without removing it from the input stream.

-  28.02.94 Original    By: ACRM
-  03.03.94 If filename is NULL, make file stdin
-  07.07.14 Use bl prefix for functions By: CTP
*/
INBUFFER *blOpenBufferedFile(char *filename, int maxstr)
{
   FILE     *fp;
   INBUFFER *BuffStruc = NULL;
   
   if(filename != NULL)
      fp=fopen(filename,"r");
   else
      fp=stdin;

   if(fp!=NULL)
   {
      if((BuffStruc = (INBUFFER *)malloc(sizeof(INBUFFER)))!=NULL)
      {
         BuffStruc->fp         = fp;
         BuffStruc->nlines     = 0;
         BuffStruc->maxstr     = maxstr;
         if((BuffStruc->buffer = 
             (char *)malloc(maxstr * sizeof(char))) != NULL)
            return(BuffStruc);
      }
      fclose(fp);
   }
   return(NULL);
}

/************************************************************************/
/*>BOOL blReadBufferedFile(INBUFFER *bfp, char *string, int length)
   ----------------------------------------------------------------
*//**

   \param[in]     *bfp      Pointer to a buffered file structure
   \param[in]     length    Size of output string
   \param[out]    *string   Output string read from file
   \return                       TRUE: Successful read
                               FALSE: End of file (or error)

   Reads a line from a buffered file (like fgets()).
   Blank lines in the file will be skipped.

-  28.02.94 Original    By: ACRM
-  07.03.94 Added code to skip blank lines
-  10.03.94 Added call to RePrompt() if we're reading from stdin and
            we get a blank line
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blReadBufferedFile(INBUFFER *bfp, char *string, int length)
{
   int bufflen = 0;
   
   if(bfp == NULL)
      return(FALSE);

   if(bfp->nlines == 0)
   {
      while(bufflen==0)
      {
         if(fgets(string,length,bfp->fp))
         {
            TERMINATE(string);
            bufflen=strlen(string);
         }
         else
         {
            return(FALSE);
         }

         if(!bufflen && bfp->fp == stdin)
            blRePrompt();
      }
   }
   else
   {
      strncpy(string,bfp->buffer,length);
      (bfp->nlines)--;
   }

   return(TRUE);
}

/************************************************************************/
/*>BOOL blProbeBufferedFile(INBUFFER *bfp, char *string, int length)
   -----------------------------------------------------------------
*//**

   \param[in]     *bfp      Pointer to a buffered file structure
   \param[in]     length    Size of output string
   \param[out]    *string   Output string read from file
   \return                       TRUE: Successful read
                               FALSE: End of file (or error)

   Read the next line from a buffered file without removing it from
   the input stream. Repeated calls will thus return the same string.
   The next call to ReadBufferedFile will also output the same string,
   but will remove the line from the input stream.
   Blank lines in the file will be skipped.

-  28.02.94 Original    By: ACRM
-  07.03.94 Added code to skip blank lines
-  10.03.94 Added call to RePrompt() after blank line
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blProbeBufferedFile(INBUFFER *bfp, char *string, int length)
{
   int bufflen = 0;
   
   if(bfp==NULL)
      return(FALSE);

   while(bufflen==0)
   {
      if(bfp->nlines == 0)
      {
         if(fgets(bfp->buffer,bfp->maxstr,bfp->fp))
         {
            TERMINATE(bfp->buffer);

            bufflen = strlen(bfp->buffer);
            if(bufflen)   (bfp->nlines)++;
         }
         else
         {
            return(FALSE);
         }

         if(!bufflen && bfp->fp == stdin)
            blRePrompt();
      }
   }

   strncpy(string,bfp->buffer,length);
   
   return(TRUE);
}


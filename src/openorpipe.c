/*************************************************************************

   Program:    
   File:       openorpipe.c
   
   Version:    V1.3
   Date:       18.08.98
   Function:   Open a file for writing unless the filename starts with
               a | in which case open as a pipe
   
   Copyright:  (c) SciTech Software 1997-8
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      amartin@stagleys.demon.co.uk
               
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
   V1.0  26.05.97 Original   By: ACRM
   V1.1  26.06.97 Added calls to signal()
   V1.2  27.02.98 Uses port.h
   V1.3  18.08.98 Added cast to popen() for SunOS

*************************************************************************/
/* Includes
*/
#include "port.h"    /* Required before stdio.h                         */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
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
/*>FILE *OpenOrPipe(char *filename)
   --------------------------------
   Input:   char  *filename     A file or pipe to be opened
   Returns: FILE  *             A file pointer

   Opens a file for writing unless the filename begins with a | in which
   case it is opened as a pipe.

   Broken pipe signals are ignored.

   26.05.97 Original   By: ACRM
   26.06.97 Added call to signal()
   18.08.98 Added case to popen() for SunOS
*/
FILE *OpenOrPipe(char *filename)
{
   char *fnam;
   
   KILLLEADSPACES(fnam, filename);
   if(fnam[0] == '|')
   {
      signal(SIGPIPE, SIG_IGN);
      fnam++;
      KILLLEADSPACES(fnam, fnam);
      return((FILE *)popen(fnam, "w"));
   }
   else
   {
      return(fopen(fnam, "w"));
   }
}

/************************************************************************/
/*>int CloseOrPipe(FILE *fp)
   -------------------------
   Input:   FILE  *fp        File pointer to be closed
   Returns: int              Error code (as for fclose())

   Attempts to close a file pointer as a pipe. If it isn't associated 
   with a pipe (i.e. popen returns (-1)), tries again to close it as
   a normal file.

   26.05.97 Original   By: ACRM
   26.06.97 Added call to signal()
*/
int CloseOrPipe(FILE *fp)
{
   int ret;
   
   if((ret=pclose(fp)) == (-1))
      return(fclose(fp));

   signal(SIGPIPE, SIG_DFL);
   return(ret);
}


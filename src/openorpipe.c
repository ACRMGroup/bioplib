/************************************************************************/
/**

   \file       openorpipe.c
   
   \version    V1.9
   \date       07.07.14
   \brief      Open a file for writing unless the filename starts with
               a | in which case open as a pipe
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1997-2014
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
-  V1.0  26.05.97 Original   By: ACRM
-  V1.1  26.06.97 Added calls to signal()
-  V1.2  27.02.98 Uses port.h
-  V1.3  18.08.98 Added cast to popen() for SunOS
-  V1.4  28.01.04 Added NOPIPE define. Allows compilation on systems
                  which don't support unix pipes
-  V1.5  03.02.06 Added prototypes for popen() and pclose()
-  V1.6  29.06.07 popen() and pclose() prototypes now skipped for MAC OSX
                  which defines them differently
-  V1.7  17.03.09 popen() prototype now skipped for Windows.
-  V1.8  02.04.09 Clean compile with NOPIPE defined
-  V1.9  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#ifndef NOPIPE
#include "port.h"    /* Required before stdio.h                         */
#endif

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
#if !defined(__APPLE__) && !defined(MS_WINDOWS)
FILE *popen(char *, char *);
#endif
#ifndef __APPLE__
int  pclose(FILE *);
#endif

/************************************************************************/
/*>FILE *blOpenOrPipe(char *filename)
   ----------------------------------
*//**

   \param[in]     *filename     A file or pipe to be opened
   \return                        A file pointer

   Opens a file for writing unless the filename begins with a | in which
   case it is opened as a pipe.

   Broken pipe signals are ignored.

-  26.05.97 Original   By: ACRM
-  26.06.97 Added call to signal()
-  18.08.98 Added case to popen() for SunOS
-  28.01.05 Added NOPIPE define
-  07.07.14 Use bl prefix for functions By: CTP
*/
FILE *blOpenOrPipe(char *filename)
{
   char *fnam;
   
   KILLLEADSPACES(fnam, filename);
#ifdef NOPIPE
   return(fopen(fnam, "w"));
#else
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
#endif
}

/************************************************************************/
/*>int blCloseOrPipe(FILE *fp)
   ---------------------------
*//**

   \param[in]     *fp        File pointer to be closed
   \return                      Error code (as for fclose())

   Attempts to close a file pointer as a pipe. If it isn't associated 
   with a pipe (i.e. popen returns (-1)), tries again to close it as
   a normal file.

-  26.05.97 Original   By: ACRM
-  26.06.97 Added call to signal()
-  28.01.05 Added NOPIPE define
-  02.04.09 Moved 'int ret' to be in the #else
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blCloseOrPipe(FILE *fp)
{
#ifdef NOPIPE
   return(fclose(fp));
#else
   int ret;

   if((ret=pclose(fp)) == (-1))
      return(fclose(fp));

   signal(SIGPIPE, SIG_DFL);
   return(ret);
#endif
}


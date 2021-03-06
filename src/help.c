/************************************************************************/
/**

   \file       help.c
   
   \version    V1.4
   \date       07.07.14
   \brief      Provides a simple file-based command line help utility
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2014
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

   A pair of routines for handling a help facility. All screen/keyboard
   I/O is via screen() and GetKybdString() routines.

**************************************************************************

   Usage:
   ======

\code
   blDoHelp(string,helpfilename)
\endcode

   This is the normal entry point and should
   be supplied with the complete command 
   including the word `help'. If this is the
   only word given, the routine will prompt
   with Help> and give help on each word 
   typed until return is hit to exit help. 
   If help followed by a keyword is given, 
   help only on that topic will be supplied.

\code
   blHelp(string,helpfilename)
\endcode

   Generates help from helpfilename on the 
   topic named by string. If this is `help' 
   or `?', available topics will be listed.

\code
   blHelp(NULL,"CLOSE")
\endcode

   Used to close the help file

   Under Unix, the environment variable, HELPDIR should be set to
   specify the directory in which help files are stored.
   Under operating systems such as VAX/VMS and AmigaDOS which support
   the assign command, the HELPDIR: assign should be set up for the help
   directory.

**************************************************************************

   Revision History:
   =================
-  V1.0  29.09.92 Original
-  V1.1  12.08.93 Returns correctly if help file not found
-  V1.2  04.01.94 Changed DoHelp() to fix problem with compilers which
                  don't let you write to strings defined in double 
                  inverted commas and never assigned to a variable.
                  Added getenv call for Unix support
-  V1.3  18.01.95 Help() changed to call OpenFile() rather than handling 
                  alternative directory internally. Consequently assign 
                  or envvar is called HELPDIR.
-  V1.4  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP User interaction
   #FUNCTION  blHelp()
   Generates help from a help file on the topic named by string. If 
   this is `help' or `?', available topics will be listed.

   #FUNCTION  blDoHelp()
   Handles help facility. This is the normal entry point.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "macros.h"
#include "WindIO.h"
#include "parse.h"
#include "help.h"
#include "general.h"

/************************************************************************/
/* Defines and macros
*/
#define BUFFLEN  240       /* Buffer for reading from file              */
#define HELPENV "HELPDIR"  /* The help directory (variable) for unix or
                              assign name for VMS/AmigaDOS              */

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void blHelp(char *string, char *HelpFile)
   -----------------------------------------
*//**

   \param[in]     *string     Topic on which to provide help. If "help" or
                              "?" then list all topics
   \param[in]     *HelpFile   Name of help file, or "CLOSE" to close the
                              help file

   Generates help from a help file on the topic named by string. If 
   this is `help' or `?', available topics will be listed.

   Only one help file may be open at a time. Once a file is open, no check
   is made on the HelpFile string to check that you still require the same
   help file. You must first close the file before changing to a new file.

   Calling:
      Help(NULL,"CLOSE")
   will close the help file.

   The directory specified by the environment variable or assign
   specified by the #define HELPENV (normally HELPDIR) will be searched 
   for the help file if not found in the local directory.
   
-  25.09.92 Original
-  28.09.92 Added HELP
-  29.09.92 Changed to file-based version
-  07.10.92 Added paging support
-  12.08.93 Correctly returns if help file not found
-  05.01.94 Handles getenv() call for Unix
-  11.03.94 Resets FirstCall to TRUE when file is closed
-  18.01.95 Calls OpenFile() rather than handling alternative directory
            internally. Consequently assign or envvar is called HELPDIR.
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blHelp(char *string,
            char *HelpFile)
{
   int         nletters,
               buffpos,
               i,
               Found;
   static int  FirstCall = TRUE;
   char        FileBuff[BUFFLEN],
               buffer[80],
               *ptr;
   static FILE *fp       = NULL;
   
   if(!strcmp(HelpFile,"CLOSE"))
   {
      if(fp)
         fclose(fp);
      FirstCall = TRUE;
      return;
   }
   
   /* If first call, open file                                          */
   if(FirstCall)
   {
      BOOL NoEnv;

      if((fp=blOpenFile(HelpFile, HELPENV, "r", &NoEnv))==NULL)
      {
         blScreen("   Error==> Unable to open help file.\n");
         sprintf(FileBuff,"            The %s environment variable \
or assign has not been set.\n",HELPENV);
         blScreen(FileBuff);
         return;
      }

      FirstCall = FALSE;
   }
   
   /* Rewind the file                                                   */
   rewind(fp);

   /* If asking from general help, display known commands               */
   if(blMatch(string,"HELP",&nletters) || string[0] == '?')
   {
      /* Search the file for keywords, echoing them to the screen       */
      buffpos = 0;

      while(fgets(FileBuff, BUFFLEN, fp))
      {
         TERMINATE(FileBuff);
         if(FileBuff[0] == '#')
         {
            if(buffpos + strlen(FileBuff) > 58)
            {
               buffer[buffpos] = '\0';

               blScreen("   ");
               blScreen(buffer);
               blScreen("\n");
               buffpos = 0;
            }
            for(i=1; i<strlen(FileBuff); i++)
            {
               buffer[buffpos++] = FileBuff[i];
            }
            buffer[buffpos++] = ' ';
         }
      }

      /* Print the last line                                            */
      if(buffpos != 0)
      {
         buffer[buffpos] = '\0';

         blScreen("   ");
         blScreen(buffer);
         blScreen("\n");
      }
   }
   else  /* Asking for help on a specific subject                       */
   {
      Found = FALSE;
      blPagingOn();
      
      while(fgets(FileBuff, BUFFLEN, fp))
      {
         TERMINATE(FileBuff);
         if(FileBuff[0] == '#')
         {
            ptr = FileBuff+1;
            UPPER(ptr);
            if(blMatch(string,ptr,&nletters))
            {
               Found = TRUE;
               while(fgets(FileBuff, BUFFLEN, fp))
               {
                  TERMINATE(FileBuff);
                  if(FileBuff[0] == '#') break;

                  blScreen(FileBuff);
                  blScreen("\n");
               }
            }
         }
      }
      if(!Found)
      {
         blScreen("   Sorry, no help on '");
         blScreen(string);
         blScreen("'\n");
      }
      blPagingOff();
   }
}

/************************************************************************/
/*>void blDoHelp(char *string, char *HelpFile)
   -------------------------------------------
*//**

   \param[in]     *string   String on which to give help, must include
                            the word "help". If given on its own, sits
                            in a loop prompting for help.
   \param[in]     *HelpFile The help file to search

   Handles help facility.
   This is the normal entry point and should be supplied with the 
   complete command including the word `help'. 

   e.g. DoHelp("help","foo.hlp");
   or   DoHelp("help bar","foo.hlp");
   
   If help is the only word given, the routine will prompt with Help> 
   and give help on each word typed until return is hit to exit help. If 
   help followed by a keyword is given, help only on that topic will be 
   supplied.
   
-  25.09.92 Original
-  28.09.92 Changed to call Help("Help")
-  02.10.92 Added GetKybdString() rather than gets() 
-  04.01.94 Changed to fix problem with compilers which
            don't let you write to strings defined in double 
            inverted commas and never assigned to a variable
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blDoHelp(char *string,
              char *HelpFile)
{
   int   i;
   char  *str,
         buffer[160];
   
   /* Trim any trailing spaces                                          */
   strcpy(buffer,string);   /* 04.01.94: Put string in buffer           */
   for(i=strlen(buffer)-1;buffer[i]==' '||buffer[i]=='\t';i--);
   buffer[++i]='\0';
   
   if((str=strchr(buffer,' '))!=NULL) /* See if a space is in the string*/
   {
      /* Yes, set pointer to the position after the space 
         (i.e. the keyword)
      */
      str++;
   }
   else
   {
      /* No keyword was specified, so give help on help                 */
      blHelp("Help",HelpFile);
   }
   
   
   if(str)  /* If specified, just give help on the keyword              */
   {
      blHelp(str,HelpFile);
   }
   else     /* Sit in a loop handling each keyword                      */
   {
      for(;;)
      {
         blPrompt("Help");
         blGetKybdString(buffer, 160);

         TERMINATE(buffer);
         if(buffer[0])
            blHelp(buffer,HelpFile);
         else
            break;
      }
   }
}

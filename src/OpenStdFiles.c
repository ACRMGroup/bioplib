/************************************************************************/
/**

   \file       OpenStdFiles.c
   
   \version    V1.22
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1991-2014
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
-  V1.1  08.02.91 Added KillLine()
-  V1.2  10.02.91 Added setextn() and index()
-  V1.3  20.03.91 Added Word()
-  V1.4  28.05.92 ANSIed
-  V1.5  22.06.92 Added tab check to Word(). Improved setextn().
                  Added WordN(). Documented other routines.
-  V1.6  27.07.93 Corrected fsscanf() for double precision
-  V1.7  07.10.93 Checks made on case before toupper()/tolower()
                  for SysV compatibility. Also index() becomes
                  chindex()
-  V1.8  18.03.94 getc() -> fgetc()
-  V1.9  11.05.94 Added GetFilestem(), upstrcmp(), upstrncmp() &
                  GetWord()
-  V1.10 24.08.94 Added OpenStdFiles()
-  V1.11 08.03.95 Corrected OpenFile() for non-UNIX
-  V1.12 09.03.95 Added check on non-NULL filename in OpenFile()
-  V1.13 17.07.95 Added countchar()
-  V1.14 18.10.95 Moved YorN() to WindIO.c
-  V1.15 06.11.95 Added StoreString(), InStringList() and FreeStringList()
-  V1.16 22.11.95 Moved ftostr() to generam.c
-  V1.17 15.12.95 Added QueryStrStr()
-  V1.18 18.12.95 OpenStdFiles() treats filename of - as stdin/stdout
-  V1.19 05.02.96 OpenStdFiles() allows NULL pointers instead if filenames
-  V1.20 18.09.96 Added padchar()
-  V1.21 18.06.02 Added string.h
-  V1.22 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP File IO
   #FUNCTION  blOpenStdFiles()
   Open the files if specified. Does not modify the file handles if 
   files are not specified. Typically used to open files for input and
   output using stdin and stdout if files not given.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "SysDefs.h"

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
/*>BOOL blOpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
   -----------------------------------------------------------------------
*//**

   \param[in]     *infile     Input filename
   \param[in]     *outfile    Output filename
   \param[out]    **in        Input file pointer
   \param[out]    **out       Output file pointer
   \return                    Success?

   Open the files if specified. In and out are not modified if files
   are not specified.

-  29.06.94 Original    By: ACRM
-  24.08.94 Name changed from OpenFiles() and placed in gen lib.
-  18.12.95 Now treats a filename of - as stdin/stdout
-  05.02.96 Allows NULL pointers
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blOpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
{
   if(infile!=NULL && infile[0] && strcmp(infile,"-"))
   {
      if((*in = fopen(infile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",infile);
         return(FALSE);
      }
   }
      
   if(outfile!=NULL && outfile[0] && strcmp(outfile,"-"))
   {
      if((*out = fopen(outfile,"w"))==NULL)
      {
         fprintf(stderr,"Unable to open output file: %s\n",outfile);
         return(FALSE);
      }
   }
   
   return(TRUE);
}



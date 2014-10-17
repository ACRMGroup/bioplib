/************************************************************************/
/**

   \file       wrapprint.c
   
   \version    V1.2
   \date       07.07.14
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, University of Reading, 2002-14
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
-  V1.2  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP String handling
   #ROUTINE  blWrapString()
   Wraps a string with double inverted commas if it contains spaces
   and escapes any contained double inverted commas with a backslash.

   #ROUTINE blWrapPrint()
   Wraps a string with double inverted commas if it contains spaces
   and escapes any contained double inverted commas with a backslash.
   Then prints the string
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SysDefs.h"
#include "general.h"

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
/*>BOOL blWrapString(char *in, char *out, int maxlen)
   --------------------------------------------------
*//**

   \param[in]     *in      Input string
   \param[in]     maxlen   Max length of output string
   \param[out]    *out     Output wrapped string
   \return                 Output string was long enough

   Wraps a string with double inverted commas if it contains spaces
   and escapes any contained double inverted commas with a backslash.
   If the output string wasn't big enough, then the routine just
   returns FALSE without copying anything into the output string.

-  30.05.02 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blWrapString(char *in, char *out, int maxlen)
{
   int  len,
        ndic,
        i;
   char *chp;
   BOOL hasWhitespace = FALSE;
   
   out[0] = '\0';

   /* Get the length of the input string                                */
   len = strlen(in);
   
   /* If the string is blank, just return two double inverted commas    */
   if(!len)
   {
      if(maxlen<3)
         return(FALSE);
      strcpy(out, "\"\"");
      return(TRUE);
   }

   /* Increment len to account for the terminating NUL                  */
   len++;
   
   /* See if the string contains spaces - if so we need to wrap the
      string with double inverted commas. Set the flag to say there
      is whitespace and increment the string length by two to
      account for the two double inverted commas.
   */
   if(strchr(in,' ') || strchr(in,'\t'))
   {
      hasWhitespace = TRUE;
      len += 2;
   }
   
   /* See if the string has double inverted commas - if so we need
      to escape them, so increment the string length
   */
   ndic = blCountchar(in, '"');
   len += ndic;

   /* Check that there is space for our padded string                   */
   if(len > maxlen)
      return(FALSE);
   
   /* Now create the output string                                      */
   i=0;
   
   if(hasWhitespace)
   {
      /* Put in the leading double inverted commas                      */
      out[i++] = '"';
   }
   
   /* Copy the string one character at a time, escaping the double 
      inverted commas
   */
   for(chp=in; *chp; chp++)
   {
      if(*chp == '"')
      {
         out[i++] = '\\';
         out[i++] = '"';
      }
      else
      {
         out[i++] = *chp;
      }
   }

   if(hasWhitespace)
   {
      /* Put in the traling double inverted commas                      */
      out[i++] = '"';
   }

   /* Terminate the string                                              */
   out[i] = '\0';
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL blWrapPrint(FILE *out, char *string)
   -----------------------------------------
*//**

   \param[in]     *out     Output file pointer
   \param[in]     *string  String to be printed
   \return                     OK?

   Wraps a string with double inverted commas if it contains spaces
   and escapes any contained double inverted commas with a backslash.
   Allocates memory for temporary storage of the wrapped string.
   Returns FALSE if this memory allocation failed.

-  30.05.02 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blWrapPrint(FILE *out, char *string)
{
   int len;
   char *buffer;
   
   len = 2 * strlen(string);
   if(len < 8)
      len=8;
   
   if((buffer=(char *)malloc(len*sizeof(char)))==NULL)
      return(FALSE);
   
   if(!blWrapString(string, buffer, len))
   {
      free(buffer);
      return(FALSE);
   }
   
   fprintf(out, "%s ", buffer);
   free(buffer);
   
   return(TRUE);
}


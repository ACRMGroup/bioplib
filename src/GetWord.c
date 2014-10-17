/************************************************************************/
/**

   \file       GetWord.c
   
   \version    V2.2
   \date       03.08.14
   \brief      Get a space delimited word from a string
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1995-2014
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
-  V1.0  02.03.99 Original   By: ACRM
-  V2.0  10.06.99 Complete rewrite to allow escaping of characters
-  V2.1  07.07.14 Use bl prefix for functions By: CTP
-  V2.2  08.03.14 Made doGetWord() a static function. By CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP String handling

   #ROUTINE  blGetWord()
   Reads a whitespace/comma delimted word out of buffer into word.

   #ROUTINE  blGetWordNC()
   Reads a whitespace delimted word out of buffer into word. Commas
   are treated just like normal characters.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "macros.h"
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
static char *doGetWord(char *buffer, char *word, int maxlen, BOOL comma);


/************************************************************************/
/*>static char *doGetWord(char *buffer, char *word, int maxlen,
                          BOOL comma)
   ------------------------------------------------------------
*//**

   \param[in]     *buffer     Input buffer to read words from
   \param[in]     maxlen      Max length of output word
   \param[in]     comma       Treat commas like white space?
   \param[out]    *word       Word read from buffer
   \return                        Pointer to start of next word in buffer
                                or NULL

   This code is designed to be called from GetWord() or GetWordNC()

   Reads a whitespace delimted word out of buffer into word. If comma is
   TRUE, then commas are treated just like white space, otherwise they
   are treated like normal characters.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

-  10.06.99 Original   By: ACRM (based on code from Bioplib)
-  03.08.14 Made static By: CTP
*/
static char *doGetWord(char *buffer, char *word, int maxlen, BOOL comma)
{
   int  i, j;
   BOOL dic    = FALSE,
        escape = FALSE;
   char *chp;
   
   /* Decrement maxlen so we can terminate correctly                    */
   maxlen--;
   
   /* Check validity of passed pointers                                 */
   if(word==NULL)
      return(NULL);

   word[0] = '\0';
   if(buffer==NULL)
      return(NULL);
   
   KILLLEADSPACES(chp, buffer);

   /* Run through each character in the input buffer                    */
   for(i=0, j=0; chp[i]; i++)
   {
      switch(chp[i])
      {
      case '\\':
         /* Use backslash as an escape character. If we've just had an
            escape, then simply store it
         */
         if(escape)
         {
            escape = FALSE;
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            escape = TRUE;
         }
         break;
      case '\"':
         /* Double inverted commas enclose strings containing white space
            If we've just had an escape then handle as a normal character,
            otherwise, toggle the dic flag
         */
         if(escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            TOGGLE(dic);
         }
         escape = FALSE;
         break;
      case ',':
         /* A comma is handled as white space or a normal character,
            depending on the comma flag
         */
         if(!comma)   /* Treat as default                               */
         {
            if(j<maxlen)
               word[j++] = chp[i];
            escape = FALSE;
            break;
         }
         /* Otherwise, if comma is true, just fall through to treat it
            like whitespace
         */
      case ' ':
      case '\t':
         /* If we are in double inverted commas or last char was an escape
            just handle as a normal character
         */
         if(dic || escape)
         {
            if(j<maxlen)
               word[j++] = chp[i];
         }
         else
         {
            /* Otherwise, this terminates the word, so terminate, move 
               the pointer on and return
            */
            word[j] = '\0';
            chp += i;
            KILLLEADSPACES(chp, chp);
            if(comma)
            {
               /* If we are handling commas as whitespace, then k
                  the comma if found      
               */
               if(*chp == ',') chp++;
            }
            if(*chp == '\0') chp = NULL;
            return(chp);
         }
         escape = FALSE;
         break;
      default:
         /* A normal character, copy it across                          */
         if(j<maxlen)
            word[j++] = chp[i];
         escape = FALSE;
      }
   }

   word[j] = '\0';
   return(NULL);
}

/************************************************************************/
/*>char *blGetWord(char *buffer, char *word, int maxlen)
   -----------------------------------------------------
*//**

   \param[in]     *buffer     Input buffer to read words from
   \param[in]     maxlen      Max length of output word
   \param[out]    *word       Word read from buffer
   \return                        Pointer to start of next word in buffer
                                or NULL

   This code is a wrapper to doGetWord()

   Reads a whitespace/comma delimted word out of buffer into word.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

-  10.06.99 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blGetWord(char *buffer, char *word, int maxlen)
{
   return(doGetWord(buffer, word, maxlen, TRUE));
}

/************************************************************************/
/*>char *blGetWordNC(char *buffer, char *word, int maxlen)
   -------------------------------------------------------
*//**

   \param[in]     *buffer     Input buffer to read words from
   \param[in]     maxlen      Max length of output word
   \param[out]    *word       Word read from buffer
   \return                    Pointer to start of next word in buffer
                              or NULL

   This code is a wrapper to doGetWord()

   Reads a whitespace delimted word out of buffer into word. Commas
   are treated just like normal characters.

   Words containing white space may be wrapped in double inverted commas.
   A \ is used as an escape character and maybe used to escape *any*
   following character. In particular:
      "\\" -> '\'     To get a backslash
      "\ " -> ' '     To get a hard whitespace (alternatively wrap the
                      string in double inverted commas)
      "\"" -> '"'     To get a double inverted comma

-  10.06.99 Original By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blGetWordNC(char *buffer, char *word, int maxlen)
{
   return(doGetWord(buffer, word, maxlen, FALSE));
}


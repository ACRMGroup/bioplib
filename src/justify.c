/************************************************************************/
/**

   \file       justify.c
   
   \version    V1.1
   \date       18.06.02
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, University of Reading, 2002
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
-  V1.0  30.05.02 Original
-  V1.1  18.06.02 Added string.h

*************************************************************************/
/* Includes
*/
#include <string.h>

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
/*>void RightJustify(char *string)
   -------------------------------
   I/O:     char    *string           A string 

   Right justifies a string in place

   30.05.02 Original   By: ACRM
*/
void RightJustify(char *string)
{
   int len, dec;

   len = strlen(string);
   len--;
   
   if(len)
   {
      while(string[len] == ' ')
      {
         for(dec = len; dec; dec--)
         {
            string[dec] = string[dec-1];
         }
         string[0] = ' ';
      }
   }
}


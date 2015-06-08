/************************************************************************/
/**

   \file       StringCat.c
   
   \version    V1.1
   \date       04.06.15
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2015
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
-  V1.0  26.03.15 Original   By:ACRM
-  V1.1  04.06.15 Fixed bug in counting number of characters to copy

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP String handling
   #FUNCTION  blStrncat()
   Like strncat() but takes the max number of characters that the
   output string can hold rather than the maximum number of characterss
   to be appended.
*/
/************************************************************************/
/* Includes
*/
#include <ctype.h>
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
/*>char *blStrncat(char *out, const char *in, size_t len)
   ------------------------------------------------------
*//**

   \param[in]      *in     Input string to be appended
   \param[in]      len     Length of output string
   \param[in,out]  *out    String that we are appending to

   A simpler version of strncat. strncat takes the max number of chars
   to be appended whereas this takes the max number of chars that 'out'
   can hold.

-  16.01.15  Original   By: ACRM
-  04.06.15  Fixed - was subtracting the input string length as well!  
*/
char *blStrncat(char *out, const char *in, size_t len)
{
   int lenOut, cpLen;
   
   lenOut = strlen(out);
   cpLen  = len - lenOut;
   
   strncat(out, in, cpLen);
   return(out);
}

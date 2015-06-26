/************************************************************************/
/**

   \file       strcatalloc.c
   
   \version    V1.3
   \date       26.06.15
   
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
-  V1.0  22.05.99 Original   By: ACRM
-  V1.1  11.07.00 Check that realloc succeeded
-  V1.2  07.07.14 Use bl prefix for functions By: CTP
-  V1.3  26.06.15 Corrected checks on input strings being NULL  By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP String handling

   #FUNCTION  blStrcatalloc()
   Like strcat() but uses a realloc() on instr to make space available.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
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
/*>char *blStrcatalloc(char *instr, char *catstr)
   ----------------------------------------------
*//**

   \param[in]     *instr    String to append to
   \param[in]     *catstr   String to append
   \return                  realloc'd version of instr with catstr
                            appended

   Like strcat() but uses a realloc() on instr to make space available.

-  22.05.99 Original   By: ACRM
-  16.06.99 Initialise outstr to NULL
-  25.08.99 Fixed bug where testing for NULL outstr instead of catstr
-  11.07.00 Check that realloc succeeded
-  07.07.14 Use bl prefix for functions By: CTP
-  26.06.15 Corrected checks on instr and catstr being null  By: ACRM
*/
char *blStrcatalloc(char *instr, char *catstr)
{
   int  totLen;
   char *outstr = NULL;
   
   totLen = ((instr==NULL)  ? 0 : strlen(instr)) + 
            ((catstr==NULL) ? 0 : strlen(catstr));
   if((outstr = realloc(instr, totLen+1))!=NULL)
   {
      /* If the input string was NULL, outstr will have been allocated
         using the equivalent of malloc() so must be initialized before
         concatenating
      */
      if(instr==NULL)
         outstr[0] = '\0';
      /* If the additional string was non-NULL then add it on           */
      if(catstr!=NULL)
         strcat(outstr, catstr);
   }
   
   return(outstr);
}

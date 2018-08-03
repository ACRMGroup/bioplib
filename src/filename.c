/************************************************************************/
/**

   \file       filename.c
   
   \version    V1.0
   \date       12.03.15
   \brief      Extract parts of filename
   
   \copyright  (c) Dr. Andrew C. R. Martin 2015
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
-  V1.0  12.03.15 Original

*************************************************************************/
/* Doxygen
   -------
   #GROUP     General Programming
   #SUBGROUP  Filenames
   #FUNCTION  blCheckProgName()
   Tests if a program name matches the specified name
*/
/************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdlib.h>
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

/************************************************************************/
/*>BOOL *blCheckProgName(char *progname, char *expected)
   -----------------------------------------------------
*//**

   \param[in]     *filename    A PDB filename containing a PDB code
   \return                     Does the program name match?
                               
   This routine strips the path and checks a program name.

-  12.03.15 Original    By: ACRM
*/
BOOL blCheckProgName(char *progname, char *expected)
{
   char *chp;
   
   /* Is the whole thing the expected program name                      */
   if(!strcmp(progname, expected))
      return(TRUE);

   /* Othewise, find the last / in the progname and return FALSE if it's
      not found
   */
   if((chp = strrchr(progname, '/'))==NULL)
      return(FALSE);

   /* Bump the pointer so it is the first character after the / and 
      compare with the expected name
   */
   chp++;
   if(!strcmp(chp, expected))
      return(TRUE);
   
   /* No match so the names differ                                      */
   return(FALSE);
}


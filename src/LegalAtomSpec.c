/************************************************************************/
/**

   \file       LegalAtomSpec.c
   
   \version    V1.7
   \date       11.10.99
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
-  V1.0  01.03.94 Original
-  V1.1  07.07.95 Now non-destructive
-  V1.2  17.07.95 Now checks that a number was specified as part of the
                  spec. and returns a BOOL
-  V1.3  23.10.95 Moved FindResidueSpec() from PDBList.c
-  V1.4  08.02.96 Added FindResidue() and changed FindResidueSpec() to
                  use it
-  V1.5  23.07.96 Added AtomNameMatch() and LegalAtomSpec()
-  V1.6  18.03.98 Added option to include a . to separate chain and 
                  residue number so numeric chain names can be used
-  V1.7  11.10.99 Allow a . to be used to start a number (such that the
                  default blank chain name is used). Allows negative 
                  residue numbers

*************************************************************************/
/* Includes
*/
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
/*>BOOL LegalAtomSpec(char *spec)
   ------------------------------
*//**

   Partner routine for AtomNameMatch(). Checks whether a wildcard
   specfication is legal (i.e. will not return an error when used
   with AtomNameMatch()).

   The only thing which is not legal is characters following a *

-  23.07.96 Original   By: ACRM
*/
BOOL LegalAtomSpec(char *spec)
{
   char *chp;
   
   for(chp=spec; *chp; chp++)
   {
      if(*chp == '\\')
      {
         chp++;
      }
      else if(*chp == '*')
      {
         chp++;
         if(*chp && *chp != ' ')
            return(FALSE);
      }
   }
   return(TRUE);
}


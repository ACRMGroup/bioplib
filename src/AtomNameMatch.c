/************************************************************************/
/**

   \file       AtomNameMatch.c
   
   \version    V1.8
   \date       07.07.14
   \brief      Tests for matching atom names with wild cards
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-9
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
-  V1.8  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Miscellaneous functions
   #ROUTINE  blAtomNameMatch()
   Tests whether an atom name matches an atom name specification.
   ? or % is used to match a single character
   * is used to match any trailing characters; it may not be used for
   leading characters or in the middle of a specification (e.g. *B*,
   C*2 are both illegal).
   Wildcards may be escaped with a backslash.

   #ROUTINE  blAtomNameRawMatch()
   Tests whether an atom name matches an atom name specification
   having been given a 'raw' atom name rather than the 
   massaged one. i.e. " CA " is C-alpha, "CA  " is Calcium
   Normally it checks against the second character onwards unless the
   spec starts with a < in which case it checks from the beginning of
   the string.
*/
/************************************************************************/
/* Includes
*/
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "SysDefs.h"
#include "pdb.h"

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
/*>BOOL blAtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn)
   --------------------------------------------------------------
*//**

   \param[in]     *atnam      The atom name to test
   \param[in]     *spec       The atom specification
   \param[in,out] *ErrorWarn  On input, if TRUE, this routine will
                              indicate errors.
                              On output, indicates whether there
                              was an error.
                              Note that you must be careful to supply
                              an lvalue here, you can't just use TRUE
                              or FALSE since it's modified on return.
                              NULL is allowed if you don't care about
                              errors.

   Tests whether an atom name matches an atom name specification.
   ? or % is used to match a single character
   * is used to match any trailing characters; it may not be used for
   leading characters or in the middle of a specification (e.g. *B*,
   C*2 are both illegal).
   Wildcards may be escaped with a backslash.

   For example: C* matches all carbon atoms,
                O5\* matches an atom called O5*
                ?B* matches all beta atoms

-  23.07.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blAtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   char *specp,
        *atnamp;
   
   /* Step through the specification and the atom name                  */
   for(specp=spec, atnamp = atnam; *specp; specp++, atnamp++)
   {
      switch(*specp)
      {
      case '\\':
         /* If the specification has a \ then we are escaping the next
            character, so just step on to that character
         */
         specp++;
         break;
      case '?':
         /* A query in the specification matches anything, so just
            continue
         */
         continue;
      case '*':
         /* Matches the rest of the string                              */
         if(ErrorWarn != NULL)
         {
            /* Check that there aren't any illegal characters following */
            if(*(specp+1) && *(specp+1) != ' ')
            {
               if(*ErrorWarn)
               {
                  fprintf(stderr,"Error in atom wildcard: %s\n",spec);
               }
               *ErrorWarn = TRUE;
            }
            else
            {
               *ErrorWarn = FALSE;
            }
         }
         return(TRUE);
      default:
         break;
      }

      /* If there is a mismatch return FALSE                            */
      if(*specp != *atnamp)
      {
         if(ErrorWarn != NULL)
            *ErrorWarn = FALSE;
         return(FALSE);
      }

      /* 07.06.05 If both specifications have ended with a space of 
         end of string then return TRUE. Fixed for if the atnam is
         shorter (after moving the alternate atom indicator into its
         own field)
      */
      if((*specp == ' ') && ((*atnamp == ' ') || (*atnamp == '\0')))
      {
         if(ErrorWarn != NULL)
            *ErrorWarn = FALSE;
         return(TRUE);
      }
   }

   /* There have been no errors and we don't need the error flag again  */
   if(ErrorWarn != NULL)
      *ErrorWarn = FALSE;

   /* The specification has run out, see if there are any atom characters
      left
   */
   if(*atnamp && *atnamp!=' ')
      return(FALSE);

   /* Both have ended OK, so the names match                            */
   return(TRUE);
}


/************************************************************************/
/*>BOOL blAtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn)
   -----------------------------------------------------------------
*//**

   \param[in]     *atnam      The atom name to check
   \param[in]     *spec       The atom specification
   \param[in,out] *ErrorWarn  On input, if TRUE, this routine will
                              indicate errors.
                              On output, indicates whether there
                              was an error.
                              Note that you must be careful to supply
                              an lvalue here, you can't just use TRUE
                              or FALSE since it's modified on return.
                              NULL is allowed if you don't care about
                              errors.

   Tests whether an atom name matches an atom name specification.

   This version should be given the raw atom name rather than the 
   massaged one. i.e. " CA " is C-alpha, "CA  " is Calcium

   Normally it checks against the second character onwards unless the
   spec starts with a < in which case it checks from the beginning of
   the string

   Written as a wrapper to AtomNameMatch()

-  15.02.01 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blAtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   /* If atom spec starts with a < then just bump the spec pointer, 
      otherwise bump the atom name pointer since we will look from the 
      second character of the atom name
   */
   if(*spec == '<')
   {
      spec++;
   }
   else
   {
      atnam++;
   }

   return(blAtomNameMatch(atnam, spec, ErrorWarn));
}

#ifdef TEST_MAIN
int main(int argc, char **argv)
{
   char spec[8], atnam[8];
   
   strcpy(atnam, " CA*");
   printf("Atom name '%s':\n", atnam);

   strcpy(spec,"CA");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"<CA");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C*");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"CA*");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"CA?");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C\\*");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));
   
   strcpy(spec,"C?");
   printf("'%s' matches? %s\n", spec, (blAtomNameRawMatch(atnam, spec, NULL)?"YES":"NO"));

   return(0);
}
#endif

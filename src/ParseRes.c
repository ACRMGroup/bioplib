/************************************************************************/
/**

   \file       ParseRes.c
   
   \version    V1.14
   \date       10.03.15
   \brief      Parse a residue specification
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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
-  V1.8  29.09.05 Moved ParseResSpec() into DoParseResSpec() with extra
                  param and added wrappers for ParseResSpec() and 
                  ParseResSpecNoUpper()  (Changes by Tony Lewis) By: TL
-  V1.9  05.01.12 Default behaviour of ParseResSpec() is now not to
                  upcase the chain label - there are now too many PDB
                  entries with lower case chain names for this to be
                  sensible.   By: ACRM
-  V1.10 12.10.12 insert is now a properly terminated string when there is
                  no insert
-  V1.11 28.08.13 chain is now a properly terminated string
-  V1.12 26.02.14 Parsing handles multi-letter chains. By: CTP
-  V1.13 07.07.14 Use bl prefix for functions By: CTP
-  V1.14 10.03.15 Added blPrintResSpecHelp()  By: ACRM
                  Removed blParseResSpecNoUpper() since blParseResSpec() 
                  now does this

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Miscellaneous functions
   #FUNCTION  blParseResSpec()
   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert code are not up-cased

   #FUNCTION  blDoParseResSpec()
   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Gives control over up-casing

   #FUNCTION  blPrintResSpecHelp()
   Prints a help message on the resdidue specfication format to make
   help messages more consistent

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
/*>BOOL blParseResSpec(char *spec, char *chain, int *resnum, char *insert)
   -----------------------------------------------------------------------
*//**

   \param[in]     *spec    Residue specification
   \param[out]    *chain   Chain label
   \param[out]    *resnum  Residue number
   \param[out]    *insert  Insert label
   \return                   Success?

   Note that chain and insert must be arrays of at least 2 characters,
   not character pointers

   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert are optional and will
   be set to spaces if not specified. Converts the resiude specification
   to upper case before processing.
   
   Moved the code that was here to a new function, DoParseResSpec()
   and made this function just call that new function.  See
   DoParseResSpec()'s comments for notes on previous changes.  This
   move is to allow the underlying function to have an extra parameter
   to specify whether or not the residue specification should be upper
   cased (without affecting code that calls this function).

-  29.09.05 Original   By: TL
-  05.01.12 Now behaves the same as ParseResSpecNoUpper(). There are now
            too many PDB files with lower case chain names (e.g. 1gav,
            3n9r, etc.) for the old default behaviour or up-casing 
            everything.   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blParseResSpec(char *spec, char *chain, int *resnum, char *insert)
{
   return blDoParseResSpec(spec, chain, resnum, insert, FALSE);
}


/************************************************************************/
/*>BOOL blDoParseResSpec(char *inSpec, char *chain, int *resnum, 
                         char *insert, BOOL uppercaseresspec)
   -------------------------------------------------------------
*//**

   \param[in]     *inSpec              Residue specification
   \param[out]    *chain               Chain label
   \param[out]    *resnum              Residue number
   \param[out]    *insert              Insert label
   \param[in]     uppercaseresspec     Convert spec to upper case.
   \return                             Success?

   Note that chain and insert must be arrays of at least 2 characters,
   not character pointers

   Splits up a residue specification of the form 
         [c][.]num[i]
   into chain, resnum and insert. Chain and insert are optional and will
   be set to spaces if not specified. If uppercaseresspec equals TRUE,
   the spec is upper cased before processing
   
   Multi-letter chain IDs can be parsed. Additionally, chain IDs with 
   numerical characters can be parsed if a period is used to separate the 
   chain from the residue number.
   
-  21.07.93 Original    By: ACRM
-  17.07.95 Added BOOL return
-  18.03.98 Added option to include a . to separate chain and residue
            number so numeric chain names can be used
-  29.09.05 Moved this code to from ParseResSpec() to DoParseResSpec()
            and made that function just call this new function.
            This move is to allow this underlying function to have an
            extra parameter to specify whether or not the residue
            specification should be upper cased (without affecting code
            that calls the old function). By: TL
-  12.10.12 insert is now a properly terminated string when there is
            no insert
-  28.08.12 chain  is now a properly terminated string
            The input specification is now copied so that actual strings
            can be passed into the routine as opposed to string delimited
            variables. This also removes the need for restoring the 
            string which has now been removed
-  26.02.14 Parsing handles multi-letter chains and numerical chain IDs.
            The "Extract chain from spec" section was re-written.
            If the period separator between the chain id and the residue
            number is absent then the chain id is set from any non-numeric
            lead characters. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blDoParseResSpec(char *inSpec, char *chain, int *resnum,
                      char *insert, BOOL uppercaseresspec)
{
   char  *ptr,
         *ptr2,
         spec[64];
   BOOL  /* DoRestore = FALSE, */
         retval    = TRUE,
         chain_found = FALSE;
   int   i;

   strncpy(spec, inSpec, 64);

   /* 11.10.99 Default resnum of 0                                      */
   *resnum = 0;

   /* Upper case the residue specification if it has been requested     */
   if (uppercaseresspec == TRUE)
   {
      UPPER(spec);   
   }
   KILLLEADSPACES(ptr, spec);
     
   /* Extract chain from spec.                   Added 26.02.14 By: CTP */

   /* Extract chain from spec (dot format)                              */
   for(ptr2=ptr,i=0;*ptr2;ptr2++,i++)
   {
      if(*ptr2 == '.')
      {
         /* set chain */
         if(i > 0)
         {
            strncpy(chain,ptr,i);
            chain[i] = '\0';
         }
         else 
         {
            strcpy(chain," ");
         }
         
         chain_found = TRUE;
         ptr = ptr2 + 1; /* update start point */
         break;
      }
   }
   
   /* Extract chain from spec (non-numeric lead characters)             */
   if(chain_found == FALSE)
   {
      for(ptr2=ptr,i=0;*ptr2;ptr2++,i++)
      {
         if(!isdigit(*ptr2) && (*ptr2 != '-'))
         {
            chain[i]    = *ptr2;
            chain[i+1]  = '\0';
            chain_found = TRUE;
            ptr = ptr2 + 1; /* update start point */
         }
         else
         {
            break;
         }
      }
   }
   
   /* Extract chain from spec (set chain to space)                      */
   if(chain_found == FALSE)
   {
      strcpy(chain," ");
   }

   
   /* Extract insert from spec                                          */
   insert[0] = ' ';
   insert[1] = '\0';  /* Added 12.10.12                                 */
   
   for(ptr2 = ptr; *ptr2; ptr2++)
   {
      /* 11.10.99 Now also checks that it isn't a - as the first 
         character 
      */
      if(!isdigit(*ptr2) && ((ptr2!=ptr)||(*ptr2 != '-')))
      {
         insert[0] = *ptr2;
         insert[1] = '\0';
         *ptr2   = '\0';
/*         DoRestore = TRUE; */
         break;
      }
   }
   
   /* Extract residue number from spec                                  */
   if(sscanf(ptr,"%d",resnum) == 0)
      retval = FALSE;

/*   if(DoRestore) */
/*   { */
      /* V1.1: Restore the original string                              */
/*      *ptr2 = *insert; */
/*   } */

   return(retval);
}


/************************************************************************/
/*>void blPrintResSpecHelp(FILE *fp)
   ---------------------------------
   \param[in]   *fp   File pointer

   Simply prints a help message fo how to use a residue specifier. Makes
   help messages from programs more consistent.

   10.03.15  Original   By: ACRM
*/
void blPrintResSpecHelp(FILE *fp)
{
   fprintf(fp,"resspec is a residue specification of the form \
[c[.]]nnn[i] where c is\n");
   fprintf(fp,"a (multi-character) chain label optionally followed by \
a '.' (required\n");
   fprintf(fp,"if the chain label is numeric), nnn is a residue number \
and i is an \n");
   fprintf(fp,"optional insert code.\n");
}


/*************************************************************************

   Program:    
   File:       FindResidue.c
   
   Version:    V1.9
   Date:       25.02.14
   Function:   Parse a residue specification
   
   Copyright:  (c) SciTech Software 1993-2014
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
   V1.0  01.03.94 Original
   V1.1  07.07.95 Now non-destructive
   V1.2  17.07.95 Now checks that a number was specified as part of the
                  spec. and returns a BOOL
   V1.3  23.10.95 Moved FindResidueSpec() from PDBList.c
   V1.4  08.02.96 Added FindResidue() and changed FindResidueSpec() to
                  use it
   V1.5  23.07.96 Added AtomNameMatch() and LegalAtomSpec()
   V1.6  18.03.98 Added option to include a . to separate chain and 
                  residue number so numeric chain names can be used
   V1.7  11.10.99 Allow a . to be used to start a number (such that the
                  default blank chain name is used). Allows negative 
                  residue numbers
   V1.8  24.02.14 Added BiopFindResidue(). By: CTP
   V1.9  25.02.14 Added error message for FindResidue(). By: CTP

*************************************************************************/
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
/*>PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list.
  Uses char for string and insert.

  06.02.96 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindResidue(). By: CTP
  25.02.14 Added error message. By: CTP
*/
PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
{
   char chain_a[2]  = " ",
        insert_a[2] = " ";
   
#ifdef BIOPLIB_CHECK
   fprintf(stderr, 
           "This code uses FindResidue() which is now deprecated!\n");
#endif

   chain_a[0]  = chain;
   insert_a[0] = insert;
   
   return(BiopFindResidue(pdb, chain_a, resnum, insert_a));
}

/************************************************************************/
/*>PDB *BiopFindResidue(PDB *pdb, char *chain, int resnum, char *insert)
  --------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list.
  Uses string for chain and insert.

  24.02.14 Original   By: CTP
*/
PDB *BiopFindResidue(PDB *pdb, char *chain, int resnum, char *insert)
{
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum == resnum) &&
         !strcmp(p->insert,insert) &&
         CHAINMATCH(p->chain,chain))
         return(p);
   }

   return(NULL);
}

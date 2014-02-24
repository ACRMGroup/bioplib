/*************************************************************************

   Program:    
   File:       ParseResidueSpec.c
   
   Version:    V1.9
   Date:       24.02.14
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
   V1.8  15.08.13 FindResidueSpec() modified as chain and insert now need
                  to be arrays
   V1.9  24.02.14 Now calls BiopFindResidue() By: CTP

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
/*>PDB *FindResidueSpec(PDB *pdb, char *resspec)
   ---------------------------------------------
   Input:   PDB   *pdb      PDB linked list
            char  *resspec  Residue specification
   Returns: PDB   *         Pointer to first atom of specified residue
                            (NULL if not found).

   Search a PDB linked list for a specified residue (given as
   [chain][.]num[insert])

   08.08.95 Original    By: ACRM
   08.02.96 Now calls FindResidue() to do the actual work
   15.08.13 chain[] and insert[] are now arrays because of changes to
            ParseResSpec()
   24.02.14 Now calls BiopFindResidue() By: CTP
*/
PDB *FindResidueSpec(PDB *pdb, char *resspec)
{
   char chain[8],
        insert[8];
   int  resnum;

   if(ParseResSpec(resspec, chain, &resnum, insert))
      return(BiopFindResidue(pdb, chain, resnum, insert));
   
   return(NULL);
}



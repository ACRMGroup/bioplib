/*************************************************************************

   Program:    
   File:       FindHetatmResidue.c
   
   Version:    V1.2
   Date:       25.02.14
   Function:   Parse a residue specification
   
   Copyright:  2011-2014
   Author:     Dr. Andrew C. R. Martin
   EMail:      martin@biochem.ucl.ac.uk
               
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
   V1.0  26.10.11 Original based on FindResidue.c
   V1.1  24.02.14 Added BiopFindHetatmResidue() By: CTP
   V1.2  25.02.14 Added error message for FindHetatmResidue(). By: CTP

*************************************************************************/
/* Includes
*/
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "pdb.h"
#include "SysDefs.h"
#include "MathType.h"
#include "macros.h"


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
/*>PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list, but
  requires the residue is a HETATM record.
  Uses char for chain and insert.

  26.10.11 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindHetatmResidue(). By: CTP
  25.02.14 Added error message. By: CTP
*/
PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
{
   char chain_a[2]  = " ",
        insert_a[2] = " ";
   
#ifdef BIOPLIB_CHECK
   fprintf(stderr, 
         "This code uses FindHetatmResidue() which is now deprecated!\n");
#endif

   chain_a[0]  = chain;
   insert_a[0] = insert;
   
   return(BiopFindHetatmResidue(pdb, chain_a, resnum, insert_a));
}

/************************************************************************/
/*>PDB *BiopFindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
  ------------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list, but
  requires the residue is a HETATM record
  Uses string for chain and insert.

  24.02.14 Original. By: CTP
*/
PDB *BiopFindHetatmResidue(PDB *pdb, char *chain, int resnum, char *insert)
{
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((!strncmp(p->record_type,"HETATM",6)) &&
      	 (p->resnum    == resnum) &&
         !strncmp(p->insert,insert,1) &&
         CHAINMATCH(p->chain,chain))
	 {
           return(p);
	 }
   }
   return(NULL);
}


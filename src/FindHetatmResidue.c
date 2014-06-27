/************************************************************************/
/**

   \file       FindHetatmResidue.c
   
   \version    V1.2
   \date       25.02.14
   \brief      Parse a residue specification
   
   \copyright  2011-2014
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
-  V1.0  26.10.11 Original based on FindResidue.c
-  V1.1  24.02.14 Added BiopFindHetatmResidue() By: CTP
-  V1.2  25.02.14 Added error message for FindHetatmResidue(). By: CTP
-  V1.3  07.05.14 Moved FindHetatmResidue() to deprecated.h By: CTP

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


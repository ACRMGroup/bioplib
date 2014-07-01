/************************************************************************/
/**

   \file       FindHetatmResidueSpec.c
   
   \version    V1.2
   \date       24.02.14
   \brief      Parse a residue specification
   
   \copyright  (c) 2011-2014
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
-  V1.0  26.10.11 Original based on FindResidueSpec.c
-  V1.1  28.08.13 Mofified for new ParseResSpec that terminates strings
-  V1.2  24.02.14 Now calls BiopFindResidue(). By: CTP

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
/*>PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec)
   ---------------------------------------------------
*//**

   \param[in]     *pdb      PDB linked list
   \param[in]     *resspec  Residue specification
   \return                  Pointer to first atom of specified residue
                            (NULL if not found).

   Search a PDB linked list for a specified residue (given as
   [chain]num[insert]) but limits search to HETATM residues

-  26.10.11 Original    By: ACRM
-  28.08.13 Mofified for new ParseResSpec that terminates strings
-  24.02.14 Now calls BiopFindResidue(). By: CTP
*/
PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec)
{
   char chain[8],
        insert[8];
   int  resnum;

   if(ParseResSpec(resspec, chain, &resnum, insert))    
      return(BiopFindHetatmResidue(pdb, chain, resnum, insert));
   
   return(NULL);
}


/************************************************************************/
/**

   \file       FindResidueSpec.c
   
   \version    V1.10
   \date       07.07.14
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
-  V1.8  15.08.13 FindResidueSpec() modified as chain and insert now need
                  to be arrays
-  V1.9  24.02.14 Now calls BiopFindResidue() By: CTP
-  V1.10 07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Searching the PDB linked list        
   #FUNCTION  blFindResidueSpec()
   Search a PDB linked list for a specified residue (given as
   [chain][.]num[insert])
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
/*>PDB *blFindResidueSpec(PDB *pdb, char *resspec)
   -----------------------------------------------
*//**

   \param[in]     *pdb      PDB linked list
   \param[in]     *resspec  Residue specification
   \return                     Pointer to first atom of specified residue
                            (NULL if not found).

   Search a PDB linked list for a specified residue (given as
   [chain][.]num[insert])

-  08.08.95 Original    By: ACRM
-  08.02.96 Now calls FindResidue() to do the actual work
-  15.08.13 chain[] and insert[] are now arrays because of changes to
            ParseResSpec()
-  24.02.14 Now calls BiopFindResidue() By: CTP
-  07.07.14 Now calls blFindResidue() Use bl prefix for functions By: CTP
*/
PDB *blFindResidueSpec(PDB *pdb, char *resspec)
{
   char chain[8],
        insert[8];
   int  resnum;

   if(blParseResSpec(resspec, chain, &resnum, insert))
      return(blFindResidue(pdb, chain, resnum, insert));
   
   return(NULL);
}



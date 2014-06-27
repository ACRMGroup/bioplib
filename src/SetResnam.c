/************************************************************************/
/**

   \file       SetResnam.c
   
   \version    V1.2
   \date       01.03.94
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
-  V1.1  01.03.94
-  V1.2  27.02.98 Removed unreachable break from switch()

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

#include "pdb.h"
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
/*>void SetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum,
                  char *insert, char *chain)
   ---------------------------------------------------------------------
   I/O:    PDB  *ResStart   Pointer to start of residue (linked list)
   Input:  PDB  *NextRes    Pointer to start of next residue
           char *resnam     Residue name to set
           int  resnum      Residue number to set
           char *insert     Insert label to set
           char *chain      Chain label to set

   Change the residue name, number, insert and chain for an amino acid.

   12.05.92 Original
*/
void SetResnam(PDB  *ResStart,
               PDB  *NextRes,
               char *resnam,
               int  resnum,
               char *insert,
               char *chain)
{
   PDB *p;
   
   for(p=ResStart; p && p!=NextRes; NEXT(p))
   {
      strcpy(p->resnam, resnam);
      strcpy(p->insert, insert);
      strcpy(p->chain,  chain);
      p->resnum = resnum;
   }
}


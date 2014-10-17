/************************************************************************/
/**

   \file       SetResnam.c
   
   \version    V1.3
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin, University of Reading, 
               2002-2014
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
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #ROUTINE  blSetResnam()
   Change the residue name, number, insert and chain for an amino acid.
*/
/************************************************************************/
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
/*>void blSetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum,
                    char *insert, char *chain)
   -----------------------------------------------------------------------
*//**

   \param[in,out] *ResStart   Pointer to start of residue (linked list)
   \param[in]     *NextRes    Pointer to start of next residue
   \param[in]     *resnam     Residue name to set
   \param[in]     resnum      Residue number to set
   \param[in]     *insert     Insert label to set
   \param[in]     *chain      Chain label to set

   Change the residue name, number, insert and chain for an amino acid.

-  12.05.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blSetResnam(PDB  *ResStart,
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


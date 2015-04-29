/************************************************************************/
/**

   \file       RenumAtomsPDB.c
   
   \version    V1.5
   \date       29.04.15
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2015
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
-  V1.1  01.03.94 Original
-  V1.2  27.02.98 Removed unreachable break from switch()
-  V1.3  07.07.14 Use bl prefix for functions By: CTP
-  V1.4  23.02.15 Properly handled TER card numbering and now takes an
                  additional parameter - an offset for the start of
                  numbering
-  V1.5  29.04.15 Increment atom number at start of HETATM records.
                  By: CTP
*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blRenumAtomsPDB()
   Renumber the atoms throughout a PDB linked list
*/
/************************************************************************/
/* Includes
*/
#include "macros.h"
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
/*>void blRenumAtomsPDB(PDB *pdb)
   ------------------------------
*//**

   \param[in,out] *pdb   PDB linked list to renumber
   \param[in]     offset Number for the first atom

   Renumber the atoms throughout a PDB linked list

-  01.08.93 Original
-  07.07.14 Use bl prefix for functions By: CTP
-  23.02.15 More intelligent version that allows for TER records which 
            are also numbered. Added offset parameter.
-  29.04.15 Increment atom number at start of HETATM records.  By: CTP
*/
void blRenumAtomsPDB(PDB *pdb, int offset)
{
   PDB *p, 
       *prev = NULL;
   int i;

   i=offset;

   for(p=pdb; p!=NULL; NEXT(p)) 
   {
      /* If there is a prevous atom and the chain of this atom is 
         different from that of the previous atom, bump the atom
         counter
         or
         If the current atom is a hetatm and the previous atom is a 
         standard atom then bump the atom counter.
      */
      if((prev != NULL) && 
         ( (!CHAINMATCH(p->chain, prev->chain) ) ||
           ( (!strncmp(p->record_type,    "HETATM", 6)) && 
             (!strncmp(prev->record_type, "ATOM  ", 6)) )))
      {
         i++;
      }
      p->atnum=i++;
      prev=p;
   }
}


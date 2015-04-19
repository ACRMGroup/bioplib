/************************************************************************/
/**

   \file       BuildAtomNeighbourPDBList.c
   
   \version    V1.4
   \date       19.08.14
   \brief      Build a new PDB linked list containing atos within a given
               distance of a specified residue
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1996-2014
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
-  V1.0  27.08.96 Original   By: ACRM
-  V1.1  17.11.05 Fixed freed memory access
-  V1.2  08.11.07 Moved out of mutmodel. Added check on memory allocation.
                  Uses occ rather than BVal as a flag
-  V1.3  07.07.14 Use bl prefix for functions By: CTP
-  V1.4  19.08.14 Renamed blBuildAtomNeighbourPDBListAsCopy to 
                  blBuildAtomNeighbourPDBListAsCopy() By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blBuildAtomNeighbourPDBListAsCopy()
   Builds a PDB linked list of atoms neighbouring those in a specified
   residue. 
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
/*>PDB *blBuildAtomNeighbourPDBListAsCopy(PDB *pdb, PDB *pRes,
                                          REAL NeighbDist)
   -----------------------------------------------------------
*//**

   \param[in]     *pdb        PDB linked list of whole structure
   \param[in]     *pRes       Pointer to start of residue of interest (may
                              be in a separate linked list providing it's
                              in the same coordinate frame)
   \param[in]     NeighbDist  Cutoff neighbour distance
   \return                    PDB linked list of atoms within cutoff
                              distance of the residue of interest.
                              (NULL if allocations failed)

   Builds a PDB linked list of atoms neighbouring those in a specified
   residue. The input list is unmodified.

-  27.08.96 Original   By: ACRM
-  17.11.05 Fixed freed memory access
-  08.11.07 Moved out of mutmodel. Added check on memory allocation.
            Uses occ rather than BVal as a flag
            MOVED FROM MUTMODEL
-  07.07.14 Use bl prefix for functions By: CTP
-  19.08.14 Renamed function to blBuildAtomNeighbourPDBListAsCopy() 
            By: CTP
*/
PDB *blBuildAtomNeighbourPDBListAsCopy(PDB *pdb, PDB *pRes, 
                                       REAL NeighbDist)
{
   PDB  *pdbN  = NULL,
        *prev  = NULL,
        *p, *q,
        *pNext = NULL;
   REAL DCutSq = NeighbDist * NeighbDist;
   
   /* First, simply duplicate the PDB linked list                       */
   if(!(pdbN = blDupePDB(pdb)))
      return(NULL);

   /* We'll use the Occ column as a flag                                */
   for(p=pdbN; p!=NULL; NEXT(p))
      p->occ = 0.0;
   
   /* Find the atom after the residue in which we are interested        */
   pNext = blFindNextResidue(pRes);

   /* Look at each atom in our residue in turn flagging atoms in the
      duplicate list which are in range
   */
   for(p=pRes; p!=pNext; NEXT(p))
   {
      for(q=pdbN; q!=NULL; NEXT(q))
      {
         if(DISTSQ(p,q) <= DCutSq)
            q->occ = (REAL)1.0;
      }
   }
   
   /* Now run through the linked list removing any items which do not
      have the Occ flag set
   */
   for(p=pdbN, prev=NULL; p!=NULL;)
   {
      if(p->occ < (REAL)0.5)
      {
         if(prev==NULL)
         {
            pdbN = p->next;
            free(p);
            p=pdbN;
            continue;
         }
         else
         {
            prev->next = p->next;
            free(p);
            p=prev->next;   /* 17.11.05              Fixed from p->next */
            continue;
         }
      }
      prev = p;
      NEXT(p);
   }

   /* Return the reduced list                                           */
   return(pdbN);
}



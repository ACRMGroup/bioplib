/************************************************************************/
/**

   \file       KillPDB.c
   
   \version    V1.12
   \date       30.09.17
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2017
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
-  V1.0  22.02.94 Original release
-  V1.1  23.05.94 Added FindNextChainPDB()
-  V1.2  05.10.94 KillSidechain() uses BOOL rather than int
-  V1.3  24.07.95 Added TermPDB()
-  V1.4  25.07.95 Added GetPDBChainLabels()
-  V1.5  26.09.95 Fixed bug in TermPDB()
-  V1.6  12.10.95 Added DupePDB(), CopyPDBCoords()
-  V1.7  23.10.95 Moved FindResidueSpec() to ParseRes.c
-  V1.8  10.01.96 Added ExtractZonePDB()
-  V1.9  14.03.96 Added FindAtomInRes()
-  V1.10 08.10.99 Initialised some variables
-  V1.11 07.07.14 Use bl prefix for functions By: CTP
-  V1.12 30.09.17 Added vlDeleteResiduePDB() By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION blKillPDB()
   Remove an item in the PDB linked list and re-link correctly. Generally
   better to use blDeleteAtomPDB()

   #FUNCTION blDeleteAtomPDB()
   Delete an atom from the linked list re-linking the list and returning
   the new start of the list (in case the first atom has been deleted)

   #FUNCTION blDeleteAtomRangePDB()
   Deletes a range of atoms from the linked list re-linking the list and 
   returning the new start of the list (in case the first atom has been 
   deleted)

   #FUNCTION blDeleteResiduePDB()
   Deletes a residue from the linked list re-linking the list. Returns
   a pointer to the next residue and updates the start of the list if
   it's the first residue that's been deleted
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>

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
/*>PDB *blDeleteAtomRangePDB(PDB *pdb, PDB *start, PDB *stop)
   ----------------------------------------------------------
*//**
   \param[in]    *pdb    Start of PDB linked list
   \param[in]    *start  First atom to delete
   \param[in]    *stop   Atom after the one to be deleted
   \return               New start of linked list

   Deletes a range of atoms from the linked list re-linking the list and 
   returning the new start of the list (in case the first atom has been 
   deleted)

-  17.03.15  Original  By: ACRM
*/
PDB *blDeleteAtomRangePDB(PDB *pdb, PDB *start, PDB *stop)
{
   PDB  *p,
        *prev = NULL;
   BOOL found = FALSE;
 
   /* Find the atom previous to start                                   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p==start)
      {
         found = TRUE;
         break;
      }
      
      prev = p;
   }
   if(!found)
      return(pdb);
   
   /* Now step from start to stop deleting atoms                        */
   for(p=start; p!=stop;)
   {
      PDB *next;
      next=blKillPDB(p, prev);
      p=next;
   }

   if(prev==NULL)
      return(stop);
   return(pdb);
}


/************************************************************************/
/*>PDB *blDeleteAtomPDB(PDB *pdb, PDB *atom)
   -----------------------------------------
*//**
   \param[in]     *pdb    Start of PDB linked list
   \param[in]     *atom   Atom to delete
   \return                New start of PDB linked list

   Deletes an atom from the PDB linked list. Should be called as 
   pdb=blDeleteAtomPDB(pdb, atom);
   to allow for the first atom in the linked list being deleted.

   Returns NULL of all atoms have been deleted. Returns the input
   pdb linked list unmodified if the atom isn't found.

-  17.03.15  Original   By: ACRM
*/
PDB *blDeleteAtomPDB(PDB *pdb, PDB *atom)
{
   PDB *p, 
       *next,
       *prev=NULL;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p==atom)
      {
         next = blKillPDB(atom, prev);
         /* Killed atom from start of linked list                       */
         if(prev==NULL)
            return(next);
         return(pdb);
      }
      
      prev=p;
   }

   return(pdb);
}

/************************************************************************/
/*>PDB *blKillPDB(PDB *pdb, PDB *prev)
   -----------------------------------
*//**

   \param[in]     *pdb    Pointer to item in PDB linked list to be removed
   \param[in]     *prev   Pointer to previous item in linked list
   \return                Next item in PDB linked list

   Kill an item in the PDB linked list and re-link correctly. Returns the
   next item in the list, so will be NULL when the last item in the list
   is killed.

-  12.05.92 Original
-  11.03.94 Now handles prev==NULL to delete first item in a list
-  07.07.14 Use bl prefix for functions By: CTP
-  16.03.15 Checks and removes any CONECT data  By: ACRM
*/
PDB *blKillPDB(PDB *pdb,              /* Pointer to record to kill      */
               PDB *prev)             /* Pointer to previous record     */
{
   PDB *next;

   next = pdb->next;

   blDeleteAtomConects(pdb);
   
   if(prev!=NULL)
      prev->next = next;            /* Relink the list                  */
   free(pdb);                       /* Free the item                    */

   return(next);
}

/************************************************************************/
/*>PDB *blDeleteResiduePDB(PDB **pPDB, PDB *res)
   ---------------------------------------------
*//**

   \param[in,out] *pPDB   Pointer to pointer to start of PDB linked list
   \param[in]     *res    Pointer to residue to be removed
   \return                Next residue in PDB linked list

   Remove a residue from the PDB linked list and re-link correctly.
   Returns the next item in the list, so will be NULL when the last 
   residue has been deleted.

-  30.09.17 Original
*/
PDB *blDeleteResiduePDB(PDB **pPDB, PDB *res)
{
   PDB *prevAtom,
       *nextRes, 
       *p;
   
   if(*pPDB == NULL)  return(NULL);

   /* Find the previous atom                                            */
   if(res == *pPDB)
   {
      prevAtom = NULL;
   }
   else
   {
      for(prevAtom=*pPDB; 
          ((prevAtom != NULL) && (prevAtom->next != res));
          NEXT(prevAtom));
   }

   /* Find the next residue                                             */
   nextRes = blFindNextResidue(res);
   
   /* Link the previous atom to the next residue                        */
   if(prevAtom == NULL)    /* Start of linked list                      */
   {
      *pPDB = nextRes;
   }
   else                    /* Elsewhere in the linked list              */
   {
      prevAtom->next = nextRes;
   }
   
   for(p=res; p!=nextRes; NEXT(p))
   {
      blDeleteAtomConects(p);
      free(p);
   }
   
   return(nextRes);
}


/************************************************************************/
/**

   \file       StructurePDB.c
   
   \version    V1.3
   \date       07.07.14
   \brief      Build a structured PDB representation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2009-14
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
-  V1.0   24.11.09  Original
-  V1.1   19.05.10  PDBRESIDUE and PDBCHAIN are doubly linked lists
-  V1.2   04.02.14  Use CHAINMATCH By: CTP
-  V1.3   07.07.14  Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blAllocPDBStructure()
   Takes a PDB linked list and converts it into a hierarchical structure
   of chains, residues and atoms

   #FUNCTION  blFreePDBStructure()
   Frees memory used by the hierarchical description of a PDB structure.
   Note this does not free the underlying PDB linked list

   #FUNCTION  blFindNextChain()
   Takes a PDB linked list and find the start of the next chain. This is
   similar to another Bioplib routine (blFindNextChainPDB()) which 
   terminates the first chain, but this routines doesn't terminate.
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
/*>PDBSTRUCT *blAllocPDBStructure(PDB *pdb)
   ----------------------------------------
*//**

   \param[in]     *pdb    PDB linked list
   \return                PDBSTRUCT structure containing chain
                          linked list

   Takes a PDB linked list and converts it into a hierarchical structure
   of chains, residues and atoms

-  24.11.09  Original   By: ACRM
-  19.05.10  Correctly do doubly linked lists
-  02.06.10  Fixed some initializations!
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDBSTRUCT *blAllocPDBStructure(PDB *pdb)
{
   PDB        *pdbc, *pdbr;
   PDB        *nextchain, *nextres;
   PDBSTRUCT  *pdbstruct = NULL;
   PDBCHAIN   *chain = NULL;
   PDBRESIDUE *residue = NULL;
   
   
   /* Initialization                                                    */
   if((pdbstruct = (PDBSTRUCT *)malloc(sizeof(PDBSTRUCT)))==NULL)
      return(NULL);
   pdbstruct->pdb = pdb;
   pdbstruct->chains = NULL;
   pdbstruct->extras = NULL;             /* 02.06.10                    */
   
   /* Build the chain list                                              */
   for(pdbc=pdb; pdbc!=NULL; pdbc=nextchain)
   {
      nextchain = blFindNextChain(pdbc);

      if(pdbstruct->chains == NULL)
      {
         INITPREV(pdbstruct->chains, PDBCHAIN);
         chain = pdbstruct->chains;
      }
      else
      {
         ALLOCNEXTPREV(chain, PDBCHAIN);
      }
      if(chain == NULL)
      {
         blFreePDBStructure(pdbstruct);
         return(NULL);
      }
      
      chain->start    = pdbc;
      chain->stop     = nextchain;
      chain->residues = NULL;              /* 02.06.10                  */
      chain->extras   = NULL;              /* 02.06.10                  */
      strcpy(chain->chain, pdbc->chain);
   }
   
   /* Build the residue list for each chain                             */
   for(chain = pdbstruct->chains; chain!=NULL; NEXT(chain))
   {
      for(pdbr = chain->start; pdbr != chain->stop; pdbr=nextres)
      {
         char resid[24];
         int  i, j;
         
         nextres = blFindNextResidue(pdbr);
         
         if(chain->residues == NULL)
         {
            INITPREV(chain->residues, PDBRESIDUE);
            residue = chain->residues;
         }
         else
         {
            ALLOCNEXTPREV(residue, PDBRESIDUE);
         }
         if(residue == NULL)
         {
            blFreePDBStructure(pdbstruct);
            return(NULL);
         }
         
         residue->start = pdbr;
         residue->stop  = nextres;
         strcpy(residue->chain, pdbr->chain);
         strcpy(residue->insert, pdbr->insert);
         strcpy(residue->resnam, pdbr->resnam);
         residue->resnum = pdbr->resnum;
         residue->extras = NULL;            /* 02.06.10                 */

         /* Build a residue label                                       */
         sprintf(resid, "%s.%d%s", 
                 residue->chain, residue->resnum, residue->insert);
         for(i=0, j=0; resid[i]; i++)
         {
            /* Copy the character if it's not a space, but only copy a 
               full stop if it's not going to be the first character
            */
            if(resid[i] != ' ' && (resid[i] != '.' || j>0))
            {
               residue->resid[j++] = resid[i];
            }
         }
         residue->resid[j++] = '\0';
      }
   }
   
   return(pdbstruct);
}

/************************************************************************/
/*>void blFreePDBStructure(PDBSTRUCT *pdbstruct)
   ---------------------------------------------
*//**

   \param[in]     *pdbstruct       PDBSTRUCT structure containing chain
                                   linked list

   Frees memory used by the hierarchical description of a PDB structure.
   Note this does not free the underlying PDB linked list

-  24.11.09  Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blFreePDBStructure(PDBSTRUCT *pdbstruct)
{
   PDBCHAIN *chain = NULL;

   if(pdbstruct == NULL)
      return;
   
   for(chain = pdbstruct->chains; chain!=NULL; NEXT(chain))
   {
      FREELIST(chain->residues, PDBRESIDUE);
   }
   
   FREELIST(pdbstruct->chains, PDBCHAIN);
   free(pdbstruct);
}


/************************************************************************/
/*>PDB *blFindNextChain(PDB *pdb)
   ------------------------------
*//**

   \param[in]     *pdb    PDB linked list
   \return                Pointer to start of next chain

   Takes a PDB linked list and find the start of the next chain. This is
   similar to another Bioplib routine (blFindNextChainPDB()) which 
   terminates the first chain, but this routines doesn't terminate.

-  24.11.09  Original   By: ACRM
-  04.02.14 Use CHAINMATCH By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDB *blFindNextChain(PDB *pdb)
{
   PDB  *p, *ret = NULL;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->next == NULL) || !CHAINMATCH(p->next->chain,pdb->chain))
      {
         ret = p->next;
         break;
      }
   }

   return(ret);
}

/************************************************************************/
#ifdef DEMO
#include <stdio.h>
int main(int argc, char **argv)
{
   FILE *fp = NULL;
   int natoms;
   PDB *pdb, *p;
   PDBSTRUCT *pdbs;
   PDBCHAIN *chain;
   PDBRESIDUE *residue;
   
   
   fp   = fopen(argv[1], "r");
   pdb  = blReadPDB(fp, &natoms);
   pdbs = blAllocPDBStructure(pdb);

   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      fprintf(stdout,"COMMENT: Chain %s\n", chain->chain);
      for(residue=chain->residues; residue!=NULL; NEXT(residue))
      {
         fprintf(stdout, "COMMENT: Residue %s\n", residue->resid);
         for(p=residue->start; p!=residue->stop; NEXT(p))
         {
            blWritePDBRecord(stdout, p);
         }
      }
   }

   for(chain=pdbs->chains; chain!=NULL; NEXT(chain))
   {
      fprintf(stdout,"COMMENT: Chain %s\n", chain->chain);
      for(p=chain->start; p!=chain->stop; NEXT(p))
      {
         blWritePDBRecord(stdout, p);
      }
   }

   blFreePDBStructure(pdbs);

   return(0);
}

#endif
   

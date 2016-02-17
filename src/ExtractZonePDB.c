/************************************************************************/
/**

   \file       ExtractZonePDB.c
   
   \version    V1.19
   \date       17.02.16
   \brief      PDB linked list manipulation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2016
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
-  V1.11 22.03.05 Extracted range is limited by specified residues
-  V1.12 22.03.06 Modified ExtractZonePDB() to allow non-exact ranges
-  V1.13 29.10.10 Fixed bug when end of zone was last residue in a chain
-  V1.14 04.02.14 Use CHAINMATCH By: CTP
-  V1.15 07.07.14 Use bl prefix for functions By: CTP
-  V1.16 19.08.14 Renamed function to blExtractZonePDBAsCopy() By: CTP
-  V1.17 07.10.15 Added function blExtractNotZonePDBAsCopy() and
                  blExtractNotZonePDBAsCopy() By: ACRM
-  V1.18 08.12.15 Modified blExtractZoneAsCopy() to check for exact 
                  match first

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Searching the PDB linked list
   #FUNCTION  blExtractZonePDBAsCopy()
   Reduces a PDB linked list to those residues within a specified zone
   forming a new linked list. Uses separate chain, residue number and 
   insert rather than residue specifications.

   #FUNCTION  blExtractZoneSpecPDBAsCopy()
   Reduces a PDB linked list to those residues within a specified zone
   forming a new linked list. Uses residue specifications ([c]nnn[i])
   rather than separate chain, residue number and insert.

   #FUNCTION  blExtractNotZonePDBAsCopy()
   Reduces a PDB linked list to those residues outside a specified zone
   forming a new linked list. Uses separate chain, residue number and 
   insert rather than residue specifications.

   #FUNCTION  blExtractNotZoneSpecPDBAsCopy()
   Reduces a PDB linked list to those residues outside a specified zone
   forming a new linked list. Uses residue specifications ([c]nnn[i])
   rather than separate chain, residue number and insert.
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"
#include "general.h"

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
/*>PDB *blExtractZonePDBAsCopy(PDB *inpdb, char *chain1, int resnum1, 
                               char *insert1, char *chain2, int resnum2, 
                               char *insert2)
   -----------------------------------------------------------------------
*//**

   \param[in]     *inpdb   Input PDB linked list
   \param[in]     *chain1  Start residue chain name
   \param[in]     resnum1  Start residue number
   \param[in]     *insert1 Start residue insert code
   \param[in]     *chain2  End residue chain name
   \param[in]     resnum2  End residue number
   \param[in]     *insert2 End residue insert code
   \return                 PDB linked list of the region of interest.

   Reduces a PDB linked list to those residues within a specified zone.
   Note that the PDB linked list is duplicated before extraction so
   pointers do not match those in the input PDB linked list. Excess
   records in the new PDB linked list are freed.

-  10.01.96 Original   By: ACRM
-  22.03.06 Modified to allow non-exact zones. i.e. the extracted zone
            will be the widest subset of the specified zone. So, if
            you specifiy 30-35Z and the PDB file only has 30-35B
            then that will be extracted.
-  29.10.10 Fixed extraction where end of zone matched last residue in
            a chain
-  04.02.14 Use CHAINMATCH By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  19.08.14 Renamed function to blExtractZonePDBAsCopy() By: CTP
-  08.12.15 Now checks for an exact match first to deal with 6INS that
            uses the insert code for something different! By: ACRM
-  17.02.16 Added fix on check for exact match - prev wasn't being set
            correctly if the match was to the first residue in the file

***    TODO - This doesn't deal with CONECT information properly!      ***
*/
PDB *blExtractZonePDBAsCopy(PDB *inpdb, 
                            char *chain1, int resnum1, char *insert1, 
                            char *chain2, int resnum2, char *insert2)
{
   PDB *pdb, *p, 
       *start = NULL, 
       *last  = NULL, 
       *prev  = NULL;

   /* Duplicate the PDB linked list                                     */
   if((pdb = blDupePDB(inpdb))==NULL)
      return(NULL);

   /* Try to find the exact residue specification, setting start to the
      start of that residue and prev to the last atom of the previous
      residue
   */
   if((start=blFindResidue(pdb, chain1, resnum1, insert1))!=NULL)
   {
      /* Find the previous atom                                         */
      if(start==pdb)
      {
         prev=NULL;
      }
      else
      {
         for(prev=pdb; prev->next!=start; NEXT(prev)) ;
      }
      
   }
   else
   {
      /* We didn't find an exact match so find the residue that would
         follow that residue specification.
         Again prev will point to the last atom before the first atom in
         the zone and start will point to the first atom in the zone
      */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(CHAINMATCH(p->chain,chain1) &&
            ((p->resnum > resnum1) ||
             ((p->resnum == resnum1) &&
              (p->insert[0] >= insert1[0]))))
         {
            start = p;
            break;
         }
         prev = p;
      }
   }

   /* If we didn't find anything then exit                              */
   if(start==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }


   /* See if we have an exact match to the residue requested
      last will be the last atom in that residue
   */
   if((last=blFindResidue(pdb, chain2, resnum2, insert2))!=NULL)
   {
      PDB *nextres = blFindNextResidue(last);

      /* Step to the end of this residue                                */
      for(; last->next!=nextres; NEXT(last)) ;
   }
   else
   {
      /* There was no exact match so find the one after that residue
         specification. 
         Again, last will point to the last atom in the zone

         29.10.10 Also breaks out if chain1 and chain2 are the same but
         we've now come to a different chain. This fixes a bug where
         the code wouldn't break out if resnum2 was the last residue in
         a chain By: ACRM
      */
      for(p=start; p!=NULL; NEXT(p))
      {
         if(CHAINMATCH(p->chain,chain2) && /* If same chain and...      */
            ((p->resnum > resnum2) ||      /* Residue number exceeded or*/
             ((p->resnum == resnum2) &&    /* Resnum same and...        */
              (p->insert[0] > insert2[0]))))  /* insert exceeded        */
         {
            break;
         }
         /* Both zone ends are in the same chain so, if we got here we have 
            found the right chain. If the current chain is now different 
            from chain2, then we've gone off the end of the chain
         */
         if( CHAINMATCH(chain1,chain2) &&
             !CHAINMATCH(p->chain,chain2))
         {
            break;
         }
         last = p;
      }
   }

   if(last==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Free linked list after 'last'                                     */
   if(last->next != NULL)
   {
      FREELIST(last->next, PDB);
      last->next = NULL;
   }
   
   /* Unlink 'start' from rest of linked list and free memory before 
      'start'
   */
   if(prev != NULL)
   {
      prev->next = NULL;
      FREELIST(pdb, PDB);
   }

   return(start);
}


/************************************************************************/
/*>PDB *blExtractZoneSpecPDBAsCopy(PDB *pdb, char *firstRes, 
                                   char *lastRes)
   ---------------------------------------------------------
*//**
   \param[in]   pdb       PDB linked list
   \param[in]   firstRes  Residue spec ([chain]resnum[insert])
   \param[in]   lastRes   Residue spec ([chain]resnum[insert])

   Extracts a zone from a PDB linked list, making a copy of the original
   list.

-  08.10.14  Original   By: ACRM
*/
PDB *blExtractZoneSpecPDBAsCopy(PDB *pdb, char *firstRes, char *lastRes)
{
   char chain1[8],  chain2[8],
        insert1[8], insert2[8];
   int  resnum1,    resnum2;
   PDB  *zone = NULL;

   if(blParseResSpec(firstRes, chain1, &resnum1, insert1) &&
      blParseResSpec(lastRes,  chain2, &resnum2, insert2))
   {
      zone = blExtractZonePDBAsCopy(pdb, 
                                    chain1, resnum1, insert1,
                                    chain2, resnum2, insert2);
   }
   return(zone);
}

/************************************************************************/
/*>PDB *blExtractNotZonePDBAsCopy(PDB *inpdb, char *chain1, int resnum1, 
                                  char *insert1, char *chain2, int resnum2, 
                                  char *insert2)
   -----------------------------------------------------------------------
*//**

   \param[in]     *inpdb   Input PDB linked list
   \param[in]     *chain1  Start residue chain name
   \param[in]     resnum1  Start residue number
   \param[in]     *insert1 Start residue insert code
   \param[in]     *chain2  End residue chain name
   \param[in]     resnum2  End residue number
   \param[in]     *insert2 End residue insert code
   \return                 PDB linked list of the region of interest.

   Reduces a PDB linked list to those residues outside a specified zone.
   Note that the PDB linked list is duplicated before extraction so
   pointers do not match those in the input PDB linked list. Excess
   records in the new PDB linked list are freed.

-  07.10.15 Original   By: ACRM
-  08.10.15 Fixed return(FALSE) to return(NULL)
*/
PDB *blExtractNotZonePDBAsCopy(PDB *inpdb, 
                               char *chain1, int resnum1, char *insert1, 
                               char *chain2, int resnum2, char *insert2)
{
   PDB  *pdb   = NULL,
        *start = NULL, 
        *stop  = NULL,
        *p, *q;

   /* Find the start and stop of the zone in the old linked list        */
   if((start = blFindResidue(inpdb, chain1, resnum1, insert1)) == NULL)
      return(NULL);
   if((stop  = blFindResidue(inpdb, chain2, resnum2, insert2)) == NULL)
      return(NULL);
   stop  = blFindNextResidue(stop);
   
   for(p=inpdb; p!=start; NEXT(p))
   {
      if(pdb == NULL)
      {
         INIT(pdb, PDB);
         q = pdb;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      if(q==NULL)
      {
         FREELIST(pdb, PDB);
         return(NULL);
      }
      
      blCopyPDB(q, p);
   }

   for(p=stop; p!=NULL; NEXT(p))
   {
      if(pdb == NULL)
      {
         INIT(pdb, PDB);
         q = pdb;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      if(q==NULL)
      {
         FREELIST(pdb, PDB);
         return(NULL);
      }
      
      blCopyPDB(q, p);
   }

   blCopyConects(pdb, inpdb);

   return(pdb);
}


/************************************************************************/
/*>PDB *blExtractNotZoneSpecPDBAsCopy(PDB *pdb, char *firstRes, 
                                      char *lastRes)
   ------------------------------------------------------------
*//**
   \param[in]   pdb       PDB linked list
   \param[in]   firstRes  Residue spec ([chain]resnum[insert])
   \param[in]   lastRes   Residue spec ([chain]resnum[insert])

   Extracts atoms outside a zone from a PDB linked list, making a copy 
   of the original list.

-  07.10.15  Original   By: ACRM
*/
PDB *blExtractNotZoneSpecPDBAsCopy(PDB *pdb, char *firstRes, 
                                   char *lastRes)
{
   char chain1[8],  chain2[8],
        insert1[8], insert2[8];
   int  resnum1,    resnum2;
   PDB  *zone = NULL;

   if(blParseResSpec(firstRes, chain1, &resnum1, insert1) &&
      blParseResSpec(lastRes,  chain2, &resnum2, insert2))
   {
      zone = blExtractNotZonePDBAsCopy(pdb, 
                                       chain1, resnum1, insert1,
                                       chain2, resnum2, insert2);
   }
   return(zone);
}



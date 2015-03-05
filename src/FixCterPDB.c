/************************************************************************/
/**

   \file       FixCterPDB.c
   
   \version    V1.10
   \date       05.03.15
   \brief      Routine to add C-terminal oxygens.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2015
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

   BOOL FixCterPDB(PDB *pdb, int style)

   Renames Cter oxygens in the required style (CTER_STYLE_STD,
   CTER_STYLE_GROMOS or CTER_STYLE_CHARMM) and generates a second
   oxygen if required. If the style is CHARMM an extra CTER residue
   will be created. Input may be in any of the 3 styles.

**************************************************************************

   Revision History:
   =================
-  V1.0  24.08.94 Original    By: ACRM
-  V1.1  05.10.94 Removed unused variables
-  V1.2  12.11.96 If any of the antecedant coordinates are undefined, set
                  the terminal oxygen to NULL coordinates
-  V1.3  13.11.96 Also checks for missing CA,C and O1 records
-  V1.4  06.02.03 Handles atnam_raw
-  V1.5  03.06.05 Handles altpos
-  V1.6  04.02.14 Use CHAINMATCH By: CTP
-  V1.7  07.07.14 Use bl prefix for functions By: CTP
-  V1.8  17.02.15 Added check on memory allocation   By: ACRM
-  V1.9  25.02.15 Ensures terminal oxygens are only added to amino acids 
                  not HETATM groups
-  V1.10 05.03.15 Replaced blFindEndPDB() with blFindNextResidue()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Modifying the structure
   #FUNCTION  blFixCterPDB()
   Renames C-ter atoms in required style and calls CalcCterCoords()
   as required to calculate coordinates ans splices them in.
   The input PDB linked list may have standard, CHARMM or GROMOS
   style.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF      160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void StandardiseCTers(PDB *pdb);
static BOOL SetCterStyle(PDB *start, PDB *end, int style);

/************************************************************************/
/*>BOOL blFixCterPDB(PDB *pdb, int style)
   --------------------------------------
*//**

   \param[in,out] *pdb   PDB linked list to modify
   \param[in]     style  CTER_STYLE_STD, CTER_STYLE_GROMOS or 
                         CTER_STYLE_CHARMM
   \return               Memory allocation OK?

   Renames C-ter atoms in required style and calls CalcCterCoords()
   as required to calculate coordinates ans splices them in.
   The input PDB linked list may have standard, CHARMM or GROMOS
   style.

-  24.08.94 Original    By: ACRM
-  05.10.94 Removed unused variables
-  12.11.96 If any of the antecedant coordinates are undefined, set
            the terminal oxygen to NULL coordinates
-  13.11.96 Added check on finding CA,C and O1
            If no C or O1 present, OXT not added
-  06.02.03 Handles atnam_raw
-  03.06.05 Handles altpos
-  04.02.14 Use CHAINMATCH By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  17.02.15 Checks INIT() succeeded  By: ACRM
-  25.02.15 Ensures terminal oxygens are only added to amino acids not
            HETATM groups
-  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
*/
BOOL blFixCterPDB(PDB *pdb, int style)
{
   PDB   *p,
         *start,
         *end,
         *CA   = NULL,
         *C    = NULL,
         *O1   = NULL,
         *O2   = NULL;
      
   /* Remove CTERs - we'll add them back later if required              */
   StandardiseCTers(pdb);   
   
   /* Step through the linked list                                      */
   start = pdb;
   while(start != NULL)
   {
      end = blFindNextResidue(start);
      
      /* Skip if this is not an amino acid                              */
      if(strncmp(start->record_type, "ATOM  ", 6))
      {
         start = end;
         continue;
      }

      if(end==NULL || !CHAINMATCH(end->chain,start->chain))
      {
         /* We're in a c-ter residue; find the atoms                    */
         CA = C = O1 = O2 = NULL;
         
         for(p=start; p!=end; NEXT(p))
         {
            if(!strncmp(p->atnam,"CA  ",4))
               CA = p;
            if(!strncmp(p->atnam,"C   ",4))
               C  = p;
            if(!strncmp(p->atnam,"O   ",4))
               O1 = p;
            if(!strncmp(p->atnam,"OXT ",4))
               O2 = p;
         }

         /* If C and O1 are missing, return                             */
         if(C==NULL && O1==NULL)
            return(TRUE);
         
         /* If O2 is missing, generate coordinates                      */
         if(O2 == NULL)
         {
            INIT(O2, PDB);
            if(O2==NULL)                     /* 17.02.15                */
               return(FALSE);
            
            CLEAR_PDB(O2);
            if(O1 != NULL) 
               blCopyPDB(O2, O1);
            else if(CA != NULL)
               blCopyPDB(O2, CA);
            else if(C != NULL)
               blCopyPDB(O2, C);
               
            O2->bval = 20.0;
            strcpy(O2->atnam,"OXT ");
            strcpy(O2->atnam_raw," OXT");
            O2->altpos = ' ';         /* 03.06.05                       */
            
            if(CA==NULL || C==NULL || O1==NULL ||
               ((CA->x == (REAL)9999.0) && 
                (CA->y == (REAL)9999.0) && 
                (CA->z == (REAL)9999.0)) ||
               ((C->x  == (REAL)9999.0) && 
                (C->y  == (REAL)9999.0) && 
                (C->z  == (REAL)9999.0)) ||
               ((O1->x == (REAL)9999.0) && 
                (O1->y == (REAL)9999.0) && 
                (O1->z == (REAL)9999.0)))
            {
               O2->x = O2->y = O2->z = (REAL)9999.0;
               O2->occ = O2->bval = (REAL)0.0;
            }
            else
            {
               if(!blCalcCterCoords(O2,CA,C,O1))
                  return(FALSE);
            }

            /* Splice O2 into the PDB linked list after O1              */
            if(O1 != NULL)
            {
               O2->next = O1->next;
               O1->next = O2;
            }
            else if(C != NULL)
            {
               O2->next = C->next;
               C->next = O2;
            }
            else if(CA != NULL)
            {
               O2->next = CA->next;
               CA->next = O2;
            }
         }

         /* Now set the required style                                  */
         SetCterStyle(start, end, style);
      }

      start = end;
   }
   
   return(TRUE);
}

/************************************************************************/
/*>static void StandardiseCTers(PDB *pdb)
   --------------------------------------
*//**

   \param[in,out] *pdb      PDB linked list to standardise

   Runs through a PDB linked list and corrects the C-terminal residues
   to the standard PDB style. Removes CTER residue names, but does
   not generate any missing Oxygens.

-  24.08.94 Original    By: ACRM
-  06.02.03 Handles atnam_raw
-  03.06.05 Handles altpos
-  04.02.14 Use CHAINMATCH By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
*/
static void StandardiseCTers(PDB *pdb)
{
   PDB *start,
       *end,
       *O1,
       *O2,
       *p,
       *prev = NULL;

   start = pdb;
   while(start != NULL)
   {
      end = blFindNextResidue(start);

      /* If it's a c-terminal residue or the one before CTER,
         fix the names of the oxygens
      */
      if(end==NULL                            || 
         !CHAINMATCH(end->chain,start->chain) ||
         !strncmp(end->resnam,"CTER",4))
      {
         O1 = O2 = NULL;
         
         for(p=start; p!=end; NEXT(p))
         {
            if(!strncmp(p->atnam,"O   ",4) ||
               !strncmp(p->atnam,"O1  ",4) ||
               !strncmp(p->atnam,"OT1 ",4))
               O1 = p;
            if(!strncmp(p->atnam,"OXT ",4) ||
               !strncmp(p->atnam,"O2  ",4) ||
               !strncmp(p->atnam,"OT2 ",4))
               O2 = p;
         }

         if(O1 != NULL)
         {
            strcpy(O1->atnam,"O   ");
            strcpy(O1->atnam_raw," O  ");
            O1->altpos = ' ';         /* 03.06.05                       */
         }
         
         if(O2 != NULL)
         {
            strcpy(O2->atnam,"OXT ");
            strcpy(O2->atnam_raw," OXT");
            O2->altpos = ' ';         /* 03.06.05                       */
         }
         
      }

      /* If it's a CTER residue, change it's name and number to that
         of the previous residue
      */
      if(!strncmp(start->resnam,"CTER",4))
      {
         if(prev != NULL)
         {
            for(p=start; p!=end; NEXT(p))
            {
               strcpy(p->resnam, prev->resnam);
               strcpy(p->insert, prev->insert);
               p->resnum = prev->resnum;
            }
         }
      }
      
      prev = start;
      start = end;
   }
   
}

/************************************************************************/
/*>static BOOL SetCterStyle(PDB *start, PDB *end, int style)
   ---------------------------------------------------------
*//**

   \param[in,out] *start     Start of PDB linked list for c-ter residue
   \param[in]     *end       End of PDB linked list for c-ter residue
   \param[in]     style      CTER_STYLE_STD, CTER_STYLE_GROMOS or 
                             CTER_STYLE_CHARMM
   \return                     Success?

   For a c-terminal or CTER residue, sets the required style for the
   two oxygens.

-  24.08.94 Original    By: ACRM
-  06.02.03 Handles atnam_raw
-  03.06.05 Handles altpos
*/
static BOOL SetCterStyle(PDB *start, PDB *end, int style)
{
   PDB *O1,
       *O2,
       *O2prev,
       *p;

   O1 = O2 = O2prev = NULL;
         
   for(p=start; p!=end; NEXT(p))
   {
      if(!strncmp(p->atnam,"O   ",4))
         O1 = p;
      if(!strncmp(p->atnam,"OXT ",4))
         O2 = p;
      if(p->next != NULL && !strncmp(p->next->atnam,"OXT ",4))
         O2prev = p;
   }

   /* If necessary, move O2 to the end of the residue                   */
   if(O2->next != end)
   {
      /* Unlink O2                                                      */
      if(O2prev == NULL) return(FALSE);
      O2prev->next = O2->next;
      
      /* Link O2 to end of residue                                      */
         for(p=start; p->next!=end; NEXT(p)) ;
      O2->next = p->next;
         p->next  = O2;
   }

   if(style == CTER_STYLE_GROMOS)
   {
      /* Correct the names                                              */
      if(O1 != NULL)
      {
         strcpy(O1->atnam,"O1  ");
         strcpy(O1->atnam_raw," O1 ");
         O1->altpos = ' ';            /* 03.06.05                       */
      }
      if(O2 != NULL)
      {
         strcpy(O2->atnam,"O2  ");
         strcpy(O2->atnam_raw," O2 ");
         O2->altpos = ' ';            /* 03.06.05                       */
      }
      
   }
   else if(style == CTER_STYLE_CHARMM)
   {
      /* Correct the names                                              */
      if(O1 != NULL)
      {
         strcpy(O1->atnam,"OT1 ");
         strcpy(O1->atnam_raw," OT1");
         O1->altpos = ' ';            /* 03.06.05                       */
      }
      if(O2 != NULL)
      {
         strcpy(O2->atnam,"OT2 ");
         strcpy(O2->atnam_raw," OT2");
         O2->altpos = ' ';            /* 03.06.05                       */
      }

      /* Change the name and number for O2                              */
      strcpy(O2->resnam,"CTER");
      strcpy(O2->insert," ");
      O2->resnum = start->resnum + 1;
   }

   return(TRUE);
}

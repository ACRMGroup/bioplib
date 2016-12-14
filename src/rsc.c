/************************************************************************/
/**

   \file       rsc.c
   
   \version    V1.16
   \date       14.12.16
   \brief      Modify sequence of a PDB linked list
   
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

   Takes a PDB linked list and replaces sidechains

**************************************************************************

   Usage:
   ======
   The main entry point is RepSChain() which takes a PDB linked list,
   the sequence, and names of the equivalent chi table and the reference
   coordinate file. These files are only opened on the first call.

**************************************************************************

   Revision History:
   =================
-  V1.0  12.05.92 Original
-  V1.1  09.07.93 Simplified and improved memory allocation checking
-  V1.2  11.03.94 RepSChain() now returns BOOL
-  V1.3  14.03.94 RepSChain() now returns TRUE for success and FALSE for
                  failure. THIS IS INCOMPATIBLE WITH PREVIOUS VERSIONS
                  which returned 0 for success.
                  Error messages are now placed in gRSCError rather than
                  being printed.
-  V1.4  19.05.94 Fixed bug resulting from change in RotatePDB(); now
                  calls ApplyMatrixPDB() instead
-  V1.5  05.10.94 Modified for KillSidechain() returning BOOL
-  V1.6  09.11.94 Modified to use OpenFile() fo open the data files
-  V1.7  12.08.96 Added RepOneSChain() and EndRepSChain() and a few
                  internal changes to RepSChain()
-  V1.8  15.08.96 Removed unused variables from RepOneSChain()
-  V1.9  30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.10 09.02.05 Fixed to handle atnam_raw and sensible defaults for
                  occ/bval
-  V1.11 03.06.05 Added altpos
-  V1.12 07.07.14 Use bl prefix for functions By: CTP
-  V1.13 15.08.14 Updated ReadRefCoords() to use CLEAR_PDB() By: CTP
-  V1.14 23.02.15 Modified for new blRenumAtomsPDB() which takes an 
                  offset   By: ACRM
-  V1.15 25.02.15 Sets the element type for new atoms
-  V1.16 14.12.16 FixTorsions() checks return from blCalcChi()

*************************************************************************/
/* Defines required for includes
*/
#define RSC_MAIN

/************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Modifying the structure
   #FUNCTION  blRepSChain()
   Replace sidechains. 

   #FUNCTION  blRepOneSChain()
   Replace a single sidechain. 

   #FUNCTION  blRepOneSChainForce()
   Replace a single sidechain - even if the sidechain was correct 
   already. 

   #FUNCTION  blEndRepSChain()
   Cleans up open files and memory used by the sidechain replacement
   routines.
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "SysDefs.h"
#include "MathType.h"
#include "seq.h"
#include "pdb.h"           /* Defines gRSCError                         */
#include "macros.h"
#include "array.h"
#include "angle.h"
#include "general.h"


/************************************************************************/
/* Defines and macros
*/
#define NUMAAKNOWN 20
#define MAXBUFF    80     /* Used by ReadRefCoords()                    */


/************************************************************************/
/* Globals
*/
static int  sFirstCall = TRUE;      /* Flag for things done first time  */
static FILE *sFp_RefCoords;         /* The reference coordinates file   */
static int  **sChiTab = NULL;       /* Equivalent torsions              */

/************************************************************************/
/* Prototypes
*/
static PDB *DoReplace(PDB *ResStart, PDB *NextRes, char seq, int **chitab,
                      FILE *fp_RefCoords);
static int ReplaceWithGly(PDB *ResStart, PDB *NextRes);
static int ReplaceWithAla(PDB *ResStart, PDB *NextRes);
static int ReplaceGly(PDB *ResStart, PDB *NextRes, char seq, 
                      FILE *fp_RefCoords);
static int FitByFragment(PDB *destination, PDB *fragment, PDB *mobile);
static int InsertSC(PDB *insert, PDB *ResStart, PDB *NextRes, BOOL doCB);
static int Replace(PDB *ResStart, PDB *NextRes, char seq, int **chitab,
                   FILE *fp_RefCoords);
static PDB *ReadRefCoords(FILE *fp, char seq);
static void ReadChiTable(FILE *fp, int **chitab);
static int FindChiIndex(char *resnam);
static PDB *FixTorsions(PDB *pdb, PDB *ResStart, PDB *NextRes, 
                        int **chitab);
static BOOL doRepOneSChain(PDB *pdb, char *ResSpec, char aa, 
                           char *ChiTable, char *RefCoords, BOOL force);

/************************************************************************/
/*>BOOL blRepSChain(PDB *pdb, char *sequence, char *ChiTable,
                  char *RefCoords)
   ----------------------------------------------------------
*//**

   \param[in,out] *pdb        PDB linked list to modify
   \param[in]     *sequence   The 1-letter code required for the structure
   \param[in]     *ChiTable   The equivalent Chi table
   \param[in]     *RefCoords  The reference coordinates file
   \return                      Success?

   Replace sidechains. Takes a PDB linked list and a 1-letter code 
   sequence and replaces the sidechains. Also requires filenames of the 
   two datafiles. DEL residues in the pdb linked list will be skipped 
   as will -'s in the sequence
   
-  12.05.92 Original
-  14.05.92 Corrected handling of matching DEL and -
-  21.06.93 Changed to allocate chitab using Array2D
-  14.03.94 Changed logic of return value. Now places error messages in
            gRSCError.
-  09.11.94 Uses OpenFile() to look in $(DATADIR) is the file wasn't
            found as specified.
-  12.08.96 Made static variables external to this routine as they are
            shared by RepOneSChain(). sChiTab was being allocated on every
            call instead of just the first one.
-  07.07.14 Use bl prefix for functions By: CTP
-  23.02.15 Modified for new blRenumAtomsPDB() which takes an offset

*/
BOOL blRepSChain(PDB  *pdb,         /* PDB linked list                  */
                 char *sequence,    /* Sequence 1-letter code           */
                 char *ChiTable,    /* Equivalent torsion table filename*/
                 char *RefCoords)   /* Reference coordinate filename    */
{
   int   resnum;                    /* Current residue number           */
   char  insert[8],                 /* Current insert code              */
         *seq;                      /* Pointer to aa in sequence        */
   PDB   *p,                        /* General PDB pointer              */
         *ResStart,                 /* Start of residue                 */
         *NextRes;                  /* Start of next residue            */
   BOOL   noenv = FALSE;            /* Flag for no env. var. found      */
   
   if(sFirstCall)
   {
      FILE  *fp_ChiTable;

      /* Allocate 2D array for equivalent torsions                      */
      if((sChiTab = (int **)blArray2D(sizeof(int), NUMAAKNOWN, NUMAAKNOWN)) 
         == NULL)
         return(FALSE);
   
      /* Open files                                                     */
      if((fp_ChiTable = blOpenFile(ChiTable,"DATADIR","r",&noenv)) == NULL)
      {
         if(noenv)
         {
            sprintf(gRSCError,"DATADIR environment variable not set\n");
         }
         else
         {
            sprintf(gRSCError,"Unable to open chi link table: %s\n",
                    ChiTable);
         }
         return(FALSE);
      }

      if((sFp_RefCoords = blOpenFile(RefCoords,"DATADIR","r",&noenv)) == 
         NULL)
      {
         if(noenv)
         {
            sprintf(gRSCError,"DATADIR environment variable not set\n");
         }
         else
         {
            sprintf(gRSCError,"Unable to open reference coordinates: \
%s\n", RefCoords);
         }
         return(FALSE);
      }

      /* Read the equivalent chi table and the atom table               */
      ReadChiTable(fp_ChiTable, sChiTab);
      
      /* Close the table files                                          */
      fclose(fp_ChiTable);
      
      sFirstCall = FALSE;
   }
   
   /* Step along the PDB linked list isolating a residue at a time and
      replacing it if necessary. The loop also steps along the sequence
      and will end if the sequence string ends.
   */
   for(p=pdb, seq = sequence; p && *seq != '\0'; NEXT(p))
   {
      ResStart = p;           /* Start of residue                       */
      resnum = p->resnum;
      strcpy(insert,p->insert);
      
      /* Move to end of residue and store pointer to next residue       */
      while(p->next && p->next->resnum == resnum 
                    && !strcmp(p->next->insert, insert))
         NEXT(p);
      NextRes = p->next;
      
      /* If its a DEL, step over any - in required sequence             */
      if(!strncmp(ResStart->resnam,"DEL",3))
      {
         while(*seq == '-') seq++;
      }
      else     /* Not a DEL, see if it needs replacing                  */
      {
         /* If there is a sequence mismatch, replace the residue        */
         if(*seq != blThrone(p->resnam))
            p = DoReplace(ResStart,NextRes,*seq,sChiTab,sFp_RefCoords);
         if(p == NULL) return(FALSE);

         /* Step to the next sequence item which isn't a -              */
         while(*(++seq) == '-') ;
      }
   }
   
   blRenumAtomsPDB(pdb, 1);
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL blRepOneSChain(PDB *pdb, char *ResSpec, char aa, 
                       char *ChiTable, char *RefCoords)
   -----------------------------------------------------------------------
*//**

   \param[in,out] *pdb        PDB linked list to modify
   \param[in]     *ResSpec    Residue spec for residue to replace in the
                              format [c]nnn[i]
   \param[in]     aa          The 1-letter code for the new sidechain
   \param[in]     *ChiTable   The equivalent Chi table
   \param[in]     *RefCoords  The reference coordinates file
   \return                    Success?

   Replace a single sidechain. Takes a PDB linked list, a residues
   specfication (in the form [c]nnn[i] where [c] is an optional chain
   name, nnn is a residue number and [i] is an optional insert code)
   and a 1-letter code of the required sidechain and does a simple
   maximum overlap replacement of the sidechain. Also requires filenames 
   of the two datafiles.
   
-  12.08.96 Original based on RepSChain()
-  15.08.96 Removed unused variables
-  07.07.14 Use bl prefix for functions By: CTP
-  23.02.15 Modified for new blRenumAtomsPDB() which takes an offset
-  29.09.15 Now just a wrapper to doRepOneSChain()
*/
BOOL blRepOneSChain(PDB *pdb, char *ResSpec, char aa, char *ChiTable,
                    char *RefCoords)
{
   return(doRepOneSChain(pdb, ResSpec, aa, ChiTable, RefCoords, FALSE));
}

/************************************************************************/
/*>BOOL blRepOneSChainForce(PDB *pdb, char *ResSpec, char aa, 
                            char *ChiTable, char *RefCoords)
   -----------------------------------------------------------------------
*//**

   \param[in,out] *pdb        PDB linked list to modify
   \param[in]     *ResSpec    Residue spec for residue to replace in the
                              format [c]nnn[i]
   \param[in]     aa          The 1-letter code for the new sidechain
   \param[in]     *ChiTable   The equivalent Chi table
   \param[in]     *RefCoords  The reference coordinates file
   \return                    Success?

   Replace a single sidechain. Takes a PDB linked list, a residues
   specfication (in the form [c]nnn[i] where [c] is an optional chain
   name, nnn is a residue number and [i] is an optional insert code)
   and a 1-letter code of the required sidechain and does a simple
   maximum overlap replacement of the sidechain. Also requires filenames 
   of the two datafiles.

   Replaces even if the sidechain was coorect already
   
-  29.09.15 Original - wrapper to doRepOneSChain()
*/
BOOL blRepOneSChainForce(PDB *pdb, char *ResSpec, char aa, 
                         char *ChiTable, char *RefCoords)
{
   return(doRepOneSChain(pdb, ResSpec, aa, ChiTable, RefCoords, TRUE));
}



/************************************************************************/
/*>static BOOL doRepOneSChain(PDB *pdb, char *ResSpec, char aa, 
                              char *ChiTable, char *RefCoords, BOOL force)
   -----------------------------------------------------------------------
*//**

   \param[in,out] *pdb        PDB linked list to modify
   \param[in]     *ResSpec    Residue spec for residue to replace in the
                              format [c]nnn[i]
   \param[in]     aa          The 1-letter code for the new sidechain
   \param[in]     *ChiTable   The equivalent Chi table
   \param[in]     *RefCoords  The reference coordinates file
   \param[in]     force       Replace even if it's the right AA already
   \return                    Success?

   Replace a single sidechain. Takes a PDB linked list, a residues
   specfication (in the form [c]nnn[i] where [c] is an optional chain
   name, nnn is a residue number and [i] is an optional insert code)
   and a 1-letter code of the required sidechain and does a simple
   maximum overlap replacement of the sidechain. Also requires filenames 
   of the two datafiles.
   
-  12.08.96 Original based on RepSChain()
-  15.08.96 Removed unused variables
-  07.07.14 Use bl prefix for functions By: CTP
-  23.02.15 Modified for new blRenumAtomsPDB() which takes an offset
-  29.09.15 Renamed as doRepOneSChain() and wrappers written as
            blRepOneSChain() and blRepOneSChainForce()
*/
static BOOL doRepOneSChain(PDB *pdb, char *ResSpec, char aa, 
                           char *ChiTable, char *RefCoords, BOOL force)
{
   PDB   *ResStart,                 /* Start of residue                 */
         *NextRes;                  /* Start of next residue            */
   BOOL  noenv = FALSE;             /* Flag for no env. var. found      */
   
   if(sFirstCall)
   {
      FILE  *fp_ChiTable;

      /* Allocate 2D array for equivalent torsions                      */
      if((sChiTab = (int **)blArray2D(sizeof(int), NUMAAKNOWN, NUMAAKNOWN)) 
         == NULL)
         return(FALSE);
   
      /* Open files                                                     */
      if((fp_ChiTable = blOpenFile(ChiTable,"DATADIR","r",&noenv)) == NULL)
      {
         if(noenv)
         {
            sprintf(gRSCError,"DATADIR environment variable not set\n");
         }
         else
         {
            sprintf(gRSCError,"Unable to open chi link table: %s\n",
                    ChiTable);
         }
         return(FALSE);
      }

      if((sFp_RefCoords = blOpenFile(RefCoords,"DATADIR","r",&noenv)) == 
         NULL)
      {
         if(noenv)
         {
            sprintf(gRSCError,"DATADIR environment variable not set\n");
         }
         else
         {
            sprintf(gRSCError,"Unable to open reference coordinates: \
%s\n", RefCoords);
         }
         return(FALSE);
      }

      /* Read the equivalent chi table and the atom table               */
      ReadChiTable(fp_ChiTable, sChiTab);
      
      /* Close the table files                                          */
      fclose(fp_ChiTable);
      
      sFirstCall = FALSE;
   }
   
   /* Find the specified residue and the one following                  */
   if((ResStart = blFindResidueSpec(pdb, ResSpec))==NULL)
   {
      sprintf(gRSCError,"Residue specification not in PDB list\n");
      return(FALSE);
   }
   NextRes = blFindNextResidue(ResStart);
   
   /* If there is a sequence mismatch, or force is set, 
      replace the residue              
   */
   if((aa != blThrone(ResStart->resnam)) || force)
   {
      if(DoReplace(ResStart, NextRes, aa, sChiTab, sFp_RefCoords) == 
         NULL) 
      {
         return(FALSE);
      }
   }
   
   blRenumAtomsPDB(pdb, 1);
   
   return(TRUE);
}


/************************************************************************/
/*>void blEndRepSChain(void)
   -------------------------
*//**

   Cleans up open files and memory used by the sidechain replacement
   routines.

-  12.08.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blEndRepSChain(void)
{
   if(!sFirstCall)
   {
      sFirstCall = TRUE;

      if(sFp_RefCoords != NULL)
      {
         fclose(sFp_RefCoords);
         sFp_RefCoords = NULL;
      }

      if(sChiTab != NULL)
      {
         blFreeArray2D((char **)sChiTab, NUMAAKNOWN, NUMAAKNOWN);
         sChiTab = NULL;
      }
   }
}


/************************************************************************/
/*>static PDB *DoReplace(PDB  *ResStart, PDB  *NextRes, char seq, 
                         int  **chitab, FILE *fp_RefCoords)
   --------------------------------------------------------------
*//**

   \param[in,out] *ResStart        Pointer to start of residue
   \param[in]     *NextRes         Pointer to start of next res
   \param[in]     seq              1-letter code for this aa
   \param[in]     **chitab         Equivalent chis
   \param[in]     *fp_RefCoords    Reference coordinate file
   \return                           If OK, Pointer to end of replaced 
                                  residue; NULL if error.
   
   Main dispatch routine for replacing one residue with another.

-  12.05.92 Original
-  21.06.93 Changed to use Array2D allocated chitab 
-  29.09.15 Removed check that we need to do the replacement
*/
static PDB *DoReplace(PDB  *ResStart,  /* Pointer to start of residue   */
                      PDB  *NextRes,   /* Pointer to start of next res  */
                      char seq,        /* 1-letter code for this aa     */
                      int  **chitab,   /* Equivalent chis               */
                      FILE *fp_RefCoords) /* Reference coordinate file  */
{
   int   retval = 0;
   PDB   *p;
   
   /* Check we need to do a replacement                                 */
/*   if(seq == blThrone(ResStart->resnam)) return(NULL); */
   
   if(!strncmp(ResStart->resnam,"GLY",3)) /* Replace Gly with X         */
      retval = ReplaceGly(ResStart,NextRes,seq,fp_RefCoords);
   else if(seq == 'G')                    /* Replace X with Gly         */
      retval = ReplaceWithGly(ResStart, NextRes);
   else if(seq == 'A')                    /* Replace X with Ala         */
      retval = ReplaceWithAla(ResStart, NextRes);
   else                                   /* Replace X with Y           */
      retval = Replace(ResStart,NextRes,seq,chitab,fp_RefCoords);
   
   if(retval) return(NULL);               /* Problem                    */
   
   /* Step through from ResStart to find the last replaced atom         */
   for(p=ResStart; p && p->next != NextRes; NEXT(p));
   
   return(p);
}


/************************************************************************/
/*>static int ReplaceWithGly(PDB *ResStart, PDB *NextRes)
   ------------------------------------------------------
*//**

   \param[in,out] *ResStart   Start of residue to be modified
   \param[in]     *NextRes    Pointer to start of next residue
   \return                      0: OK, 1: error

   Replace residue X with Gly. This is simply a case of trimming the 
   sidechain and changing the residue name.
   
-  12.05.92 Original
-  05.10.94 Changed for BOOL return from KillSidechain()
-  07.07.14 Use bl prefix for functions By: CTP
*/
static int ReplaceWithGly(PDB *ResStart,  /* Start of res               */
                          PDB *NextRes)   /* Start of next res          */
{
   BOOL  retval;
   
   /* Kill sidechain. Third parameter is a flag to kill CB              */
   retval = blKillSidechain(ResStart,NextRes,TRUE);
   
   /* If OK, set the residue name                                       */
   if(retval) blSetResnam(ResStart,NextRes,"GLY ",ResStart->resnum,
                          ResStart->insert,ResStart->chain);
   
   return(retval?0:1);
}


/************************************************************************/
/*>static int ReplaceWithAla(PDB *ResStart, PDB *NextRes)
   ------------------------------------------------------
*//**

   \param[in,out] *ResStart   Start of residue to be modified
   \param[in]     *NextRes    Pointer to start of next residue
   \return                      0: OK, 1: error

   Replace residue X with Ala. This is simply a case of trimming the 
   sidechain and changing the residue name.
   
-  12.05.92 Original
-  05.10.94 Changed for BOOL return from KillSidechain()
-  07.07.14 Use bl prefix for functions By: CTP
*/
static int ReplaceWithAla(PDB *ResStart,  /* Start of residue           */
                          PDB *NextRes)   /* Start of next residue      */
{
   BOOL retval;
   
   /* Kill sidechain. Third parameter is a flag not to kill CB          */
   retval = blKillSidechain(ResStart,NextRes,FALSE);
   
   /* If OK, set the residue name                                       */
   if(retval) blSetResnam(ResStart,NextRes,"ALA ",ResStart->resnum,
                          ResStart->insert,ResStart->chain);
   
   return(retval?0:1);
}


/************************************************************************/
/*>static int ReplaceGly(PDB *ResStart, PDB *NextRes, char seq,
                         FILE *fp_RefCoords)
   ------------------------------------------------------------
*//**

   \param[in,out] *ResStart     Start of residue to be modified
   \param[in]     *NextRes      Pointer to start of next residue
   \param[in]     seq           1-letter code for replacement residue
   \param[in]     *fp_RefCoords Reference coordinate file pointer
   \return                        0: OK, 1: error

   Replace a Gly with another residue type.
   
-  12.05.92 Original
-  21.06.93 Changed for new version of onethr(). Removed unused param.
-  09.07.93 Simplified allocation checking
-  07.07.14 Use bl prefix for functions By: CTP
*/
static int ReplaceGly(PDB  *ResStart,
                      PDB  *NextRes,
                      char seq,
                      FILE *fp_RefCoords)
{
   int   retval = 0,                /* Assume everything OK             */
         natoms;
   PDB   *p,                        /* General PDB pointer              */
         *q,                        /* General PDB pointer              */
         *reference = NULL,         /* Reference PDB linked list        */
         *ref_mc    = NULL,         /* Mainchain from reference PDB list*/
         *parent_mc = NULL;         /* Parent mainchain                 */
   char  *three;                    /* Three letter code                */
   
   three = blOnethr(seq);
   
   /* Read the required residue type out of the reference coordinate file.
      Returns NULL if unable to allocate memory, or residue not found.
   */
   reference = ReadRefCoords(fp_RefCoords, seq);
   if(reference == NULL)
   {
      retval = 1;
      goto Cleanup;
   }
   
   /* Build a PDB linked list from the parent N, CA, C                  */
   natoms = 0;
   for(p=ResStart; p && p!=NextRes; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4) ||
         !strncmp(p->atnam,"CA  ",4) ||
         !strncmp(p->atnam,"C   ",4))
      {
         natoms++;
         
         if(parent_mc == NULL)             /* Initialise start of list  */
         {
            INIT(parent_mc,PDB);
            q = parent_mc;
         }
         else                              /* Next item in list         */
         {
            ALLOCNEXT(q,PDB);
         }
         
         /* Check allocation                                            */
         if(q == NULL)
         {
            retval = 1;
            goto Cleanup;
         }

         
         blCopyPDB(q,p);                   /* Copy PDB item             */
      }
   }
   if(natoms != 3)                         /* Atoms missing             */
   {
      retval = 1;
      goto Cleanup;
   }
   
   
   /* Build a PDB linked list from the reference N, CA, C               */
   natoms = 0;
   for(p=reference; p; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4) ||
         !strncmp(p->atnam,"CA  ",4) ||
         !strncmp(p->atnam,"C   ",4))
      {
         natoms++;
         
         if(ref_mc == NULL)                /* Initialise start of list  */
         {
            INIT(ref_mc,PDB);
            q = ref_mc;
         }
         else                              /* Next item in list         */
         {
            ALLOCNEXT(q,PDB);
         }
         
         /* Check allocation                                            */
         if(q == NULL)
         {
            retval = 1;
            goto Cleanup;
         }
         
         blCopyPDB(q,p);                   /* Copy PDB item             */
      }
   }
   if(natoms != 3)                         /* Atoms missing             */
   {
      retval = 1;
      goto Cleanup;
   }
   
   /* Move the reference aa by fitting on the atoms stored in _mc       */
   if(FitByFragment(parent_mc, ref_mc, reference))
   {
      retval = 1;
      goto Cleanup;
   }
   
   /* Now insert the s/c into the main linked list, last parameter is a
      flag to insert the CB
   */
   if(InsertSC(reference, ResStart, NextRes, TRUE))
      retval = 1;
   
   if(retval == 0) blSetResnam(ResStart,NextRes,three,ResStart->resnum,
                               ResStart->insert,ResStart->chain);

Cleanup:
   if(reference != NULL)   FREELIST(reference, PDB);
   if(ref_mc    != NULL)   FREELIST(ref_mc,    PDB);
   if(parent_mc != NULL)   FREELIST(parent_mc, PDB);

   return(retval);
}


/************************************************************************/
/*>static int FitByFragment(PDB *destination, PDB *fragment, PDB *mobile)
   ----------------------------------------------------------------------
*//**

   \param[in]     *destination       Fragment we're fitting to
   \param[in]     *fragment          Part of mobile to fit
   \param[in,out] *mobile            Whole section to move
   Returns  int                    0:OK, 1:Error

   Calculate the best fit for moving fragment onto destination. Then move
   mobile in this manner. This is used to fit one mainchain (fragment) 
   onto another (destination) then move the whole amino acid (mobile) to 
   overlap the mainchains.
   
-  12.05.92 Original
-  19.05.94 Calls ApplyMatrixPDB() rather than RotatePDB()
-  07.07.14 Use bl prefix for functions By: CTP
*/
static int FitByFragment(PDB *destination, /* Fragment we're fitting to */
                         PDB *fragment,    /* Part of mobile to fit     */
                         PDB *mobile)      /* Whole section to move     */
{
   VEC3F cg_fragment,
         cg_destination;
   REAL  rm[3][3];                     /* Rotation matrix from fitting  */

   /* Correct the atom order for the PDB lists                          */
   destination = blShuffleBB(destination);
   fragment    = blShuffleBB(fragment);
   mobile      = blShuffleBB(mobile);

   /* Get CofG of fragment and destination                              */
   blGetCofGPDB(fragment,    &cg_fragment);
   blGetCofGPDB(destination, &cg_destination);
   
   /* Fit the fragment to the destination                               */
   if(!blFitCaCbPDB(destination, fragment, rm))
      return(1);
   
   /* Negate the vector for the fragment CofG so we can translate
      mobile to the origin.
   */
   cg_fragment.x = -cg_fragment.x;
   cg_fragment.y = -cg_fragment.y;
   cg_fragment.z = -cg_fragment.z;
   
   /* Move mobile so the CofG of fragment is at the origin              */
   blTranslatePDB(mobile, cg_fragment);
   
   /* RotatePDB mobile onto the destination                             */
   blApplyMatrixPDB(mobile, rm);
   
   /* TranslatePDB mobile back to coincide fragments and 
      destination CofG 
   */
   blTranslatePDB(mobile, cg_destination);
   
   return(0);
}

   
/************************************************************************/
/*>static int InsertSC(PDB *insert, PDB *ResStart, PDB *NextRes, 
                       BOOL doCB)
   -------------------------------------------------------------
*//**

   \param[in]     *insert    PDB linked list to insert
   \param[in]     *ResStart  Start of residue
   \param[in]     *NextRes   Start of next residue
   \param[in]     doCB       TRUE: CB will be inserted
                             FALSE: CB not inserted
   \return                      0 if OK, 1 if error.

   Inserts the sidechain from insert into the PDB linked list just before
   NextRes. Fresh allocations are made for the copy.
   
-  12.05.92 Original
-  19.06.92 Added num char param to strncmp()!
-  11.03.94 Changed doCB to BOOL
*/
static int InsertSC(PDB  *insert, 
                    PDB  *ResStart,
                    PDB  *NextRes,
                    BOOL doCB)
{
   PDB   *p,
         *start,
         *prev  = NULL;
   
   if(doCB)
   {
      /* If doCB is set, kill the CB between ResStart and NextRes 
         if found 
      */
      for(p=ResStart; p && p != NextRes; NEXT(p))
      {
         if(!strncmp(p->atnam,"CB  ",4))
         {
            if(prev == NULL) return(1);   /* No b/b atom before s/c     */
            blKillPDB(p, prev);
            break;
         }
         prev = p;
      }
   }

   /* Find the position to insert. This will be start                   */
   for(start=ResStart; start && start->next != NextRes; NEXT(start));
   
   /* Step through insert, copying non-backbone atoms into the linked list
      after start
   */
   for(p=insert; p; NEXT(p))
   {
      if(strcmp(p->atnam, "N   ") &&
         strcmp(p->atnam, "CA  ") &&
         strcmp(p->atnam, "C   ") &&
         strcmp(p->atnam, "O   ") &&
         strcmp(p->atnam, "CB  "))
      {
         ALLOCNEXT(start, PDB);
         if(start == NULL) return(1);
         blCopyPDB(start, p);
      }
      
      /* Insert the CB if required                                      */
      if(doCB && !strcmp(p->atnam,"CB  "))
      {
         ALLOCNEXT(start, PDB);
         if(start == NULL) return(1);
         blCopyPDB(start, p);
      }
   }
   
   /* Now link the end of this new section to NextRes                   */
   start->next = NextRes;
   
   return(0);
}


/************************************************************************/
/*>static int Replace(PDB *ResStart, PDB *NextRes, char seq,
                      int **chitab, FILE *fp_RefCoords)
   ---------------------------------------------------------
*//**

   \param[in,out] *ResStart     Start of residue to be modified
   \param[in]     *NextRes      Pointer to start of next residue
   \param[in]     seq           1-letter code for replacement residue
   \param[in]     **chitab      Equivalent chi table
   \param[in]     *fp_RefCoords Reference coordinate file pointer
   \return                        0: OK, 1: error

   Replace a non-Gly with another residue type.
   
-  12.05.92 Original
-  21.06.93 Changed for new version of onethr()
-  21.06.93 Changed to use Array2D allocated chitab 
-  09.07.93 Simplified allocation checking
-  05.10.94 Changed for BOOL return from KillSidechain()
-  07.07.14 Use bl prefix for functions By: CTP
*/
static int Replace(PDB  *ResStart,
                   PDB  *NextRes,
                   char seq,
                   int  **chitab,
                   FILE *fp_RefCoords)
{
   int   retval = 0,                /* Assume everything OK             */
         natoms;
   PDB   *p,                        /* General PDB pointer              */
         *q,                        /* General PDB pointer              */
         *reference = NULL,         /* Reference PDB linked list        */
         *ref_mc    = NULL,         /* Mainchain from reference PDB list*/
         *parent_mc = NULL;         /* Parent mainchain                 */
   char  *three;                    /* Three letter code                */
   
   three = blOnethr(seq);
   
   /* Read the required residue type out of the reference coordinate file.
      Returns NULL if unable to allocate memory, or residue not found.
   */
   reference = ReadRefCoords(fp_RefCoords, seq);
   if(reference == NULL)
   {
      retval = 1;
      goto Cleanup;
   }
   
   /* Build a PDB linked list from the parent N, CA, C, CB              */
   natoms = 0;
   for(p=ResStart; p && p!=NextRes; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4) ||
         !strncmp(p->atnam,"CA  ",4) ||
         !strncmp(p->atnam,"CB  ",4) ||
         !strncmp(p->atnam,"C   ",4))
      {
         natoms++;
         
         if(parent_mc == NULL)               /* Initialise start of list*/
         {
            INIT(parent_mc,PDB);
            q = parent_mc;
         }
         else                                /* Next item in list       */
         {
            ALLOCNEXT(q,PDB);
         }
         
         /* Check allocation                                            */
         if(q == NULL)
         {
            retval = 1;
            goto Cleanup;
         }
         
         blCopyPDB(q,p);                     /* Copy PDB item           */
      }
   }
   if(natoms != 4)                           /* Atoms missing           */
   {
      retval = 1;
      goto Cleanup;
   }

#ifdef DEBUG
   printf("Parent mainchain:\n");
   blWritePDB(stdout, parent_mc);
#endif   
   
   /* Build a PDB linked list from the reference N, CA, C, CB           */
   natoms = 0;
   for(p=reference; p; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4) ||
         !strncmp(p->atnam,"CA  ",4) ||
         !strncmp(p->atnam,"CB  ",4) ||
         !strncmp(p->atnam,"C   ",4))
      {
         natoms++;
         
         if(ref_mc == NULL)                  /* Initialise start of list*/
         {
            INIT(ref_mc,PDB);
            q = ref_mc;
         }
         else                                /* Next item in list       */
         {
            ALLOCNEXT(q,PDB);
         }
         
         /* Check allocation                                            */
         if(q == NULL)
         {
            retval = 1;
            goto Cleanup;
         }
         
         blCopyPDB(q,p);                       /* Copy PDB item           */
      }
   }
   if(natoms != 4)                           /* Atoms missing           */
   {
      retval = 1;
      goto Cleanup;
   }

#ifdef DEBUG
   printf("Reference mainchain:\n");
   blWritePDB(stdout, ref_mc);
#endif   
   
   /* Move the reference aa by fitting on the atoms stored in _mc       */
   if(FitByFragment(parent_mc, ref_mc, reference))
   {
      retval = 1;
      goto Cleanup;
   }

#ifdef DEBUG
   printf("Reference structure after fitting to parent:\n");
   blWritePDB(stdout, reference);
#endif   
   
   /* Fix the s/c torsion angles to match those in the parent. This may
      shuffle the atom order, so ResStart gets reset.
   */
   ResStart = FixTorsions(reference, ResStart, NextRes, chitab);
   
#ifdef DEBUG
   printf("Reference structure matching torsions:\n");
   blWritePDB(stdout, reference);
#endif   
   
   /* Kill sidechain. Third parameter is a flag kill CB                 */
   retval = ((blKillSidechain(ResStart,NextRes,TRUE)) ? 0 : 1);

   /* Now insert the s/c into the main linked list, last parameter is a
      flag to insert the CB
   */
   if(InsertSC(reference, ResStart, NextRes, TRUE))
      retval = 1;
   
   if(retval == 0) blSetResnam(ResStart,NextRes,three,ResStart->resnum,
                               ResStart->insert,ResStart->chain);
   
Cleanup:
   if(reference != NULL)   FREELIST(reference, PDB);
   if(ref_mc    != NULL)   FREELIST(ref_mc,    PDB);
   if(parent_mc != NULL)   FREELIST(parent_mc, PDB);

   return(retval);
}


/************************************************************************/
/*>static PDB *ReadRefCoords(FILE *fp, char seq)
   ---------------------------------------------
*//**

   \param[in]     *fp     Reference PDB file pointer
   \param[in]     seq     Residue for which to search
   \return                   PDB linked list containing residue information

   Reads the sidechain reference file (fp) to find residue type seq, then
   reads in the data for that sidechain into a PDB linked list which is 
   then returned.
   Returns NULL if there is a problem;

-  12.05.92 Original
-  21.06.93 Changed for new version of onethr()
-  09.07.93 Corrected check on allocations
-  04.01.94 Corrected string assignments of NULL to '\0'
-  09.02.05 Sets atnam_raw
            Sets default occ/bval to 1.0 and 20.0
-  03.06.05 Sets altpos
-  05.08.14 Use CLEAR_PDB() to set default values. By: CTP
-  25.02.15 Now sets the element type   By: ACRM
*/
static PDB *ReadRefCoords(FILE *fp,
                          char seq)
{
   PDB   *p    = NULL,
         *pdb  = NULL;
   char  buffer[MAXBUFF],
         *ptr,
         *three;
   int   done  = FALSE;
   
   rewind(fp);
   
   three = blOnethr(seq);

   /* Get lines from the file                                           */
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      ptr = buffer;

      if(!strncmp(buffer,"ATOM  ",6) && !strncmp(buffer+17, three, 3))
      {
         if(pdb == NULL)                     /* Initialise PDB list     */
         {
            INIT(pdb, PDB);
            p = pdb;
         }
         else                                /* Allocate next record    */
         {
            ALLOCNEXT(p, PDB);
         }
         
         /* Check allocation                                            */
         if(p==NULL)
         {
            if(pdb!=NULL) FREELIST(pdb,PDB);
            return(NULL);
         }
         
         /* Clear PDB                                                   */
         CLEAR_PDB(p);
         
         /* Copy the first 6 charcters into RECORD_TYPE                 */
         strncpy(p->record_type,ptr,6);
         p->record_type[6] = '\0';

         ptr += 6;

         /* Read atnum from here                                        */
         sscanf(ptr,"%d",&(p->atnum));

         ptr += 7; /* 2 spaces                                          */

         /* Copy the next 4 characters into ATNAM                       */
         strncpy(p->atnam,ptr,4);
         p->atnam[4] = '\0';

         /* 09.02.05 ...and into atnam_raw                              */
         p->atnam_raw[0] = ' ';
         strncpy(p->atnam_raw+1,ptr,3);
         p->atnam_raw[4] = '\0';

         /* 03.06.05 set alternate indicator to a blank                 */
         p->altpos = ' ';

         ptr += 4;

         /* Copy the next 4 characters into RESNAM                      */
         strncpy(p->resnam,ptr,4);
         p->resnam[4] = '\0';

         ptr += 4;

         /* Copy the next 1 character into CHAIN                        */
         strncpy(p->chain,ptr,1);
         p->chain[1] = '\0';

         ptr += 1;

         /* Read resnum from here                                       */
         sscanf(ptr,"%d",&(p->resnum));

         ptr += 4;


         /* Copy the next character into INSERT                         */
         strncpy(p->insert,ptr,1);
         p->insert[1] = '\0';

         ptr += 4;

         /* Read x from here                                            */
         sscanf(ptr,"%lf",&(p->x));

         ptr += 8;

         /* Read y from here                                            */
         sscanf(ptr,"%lf",&(p->y));

         ptr += 8;

         /* Read z from here                                            */
         sscanf(ptr,"%lf",&(p->z));

         /* We don't care about occ and BVal                            */
         p->occ  = 1.0;
         p->bval = 20.0;

         /* 25.02.15 Set the element type from atnam_raw                */
         blSetElementSymbolFromAtomName(p->element, p->atnam_raw);
         
         done = TRUE;
      }
      else
      {
         /* If we've got some data we exit as soon as another residue is
            found.
         */
         if(done) break;
      }
   }

   return(pdb);
}


/************************************************************************/
/*>static void ReadChiTable(FILE *fp, int **chitab)
   ------------------------------------------------
*//**

   \param[in]     *fp        Equivalent Chi table file pointer
   \param[out]    **chitab   Equivalent chi table.

   Read the chi table and creates the chitab array.
   
-  13.05.92 Original
-  21.06.93 Changed to use Array2D allocated chitab 
*/
static void ReadChiTable(FILE *fp,
                         int  **chitab)
{
   char  buffer[160],
         *ptr;
   int   StarFound = FALSE,
         row       = 0,
         col       = 0;
   
   while(fgets(buffer,160,fp))
   {
      if(StarFound)
      {
         ptr = buffer+4;
         for(col = 0; col<NUMAAKNOWN; col++, ptr += 4)
            sscanf(ptr,"%d",&(chitab[row][col]));
         if(row++ >= NUMAAKNOWN) break;
      }
      if(buffer[0] == '*') StarFound = TRUE;
   }

#ifdef DEBUG
   for(row=0; row<NUMAAKNOWN; row++)
   {
      for(col=0; col<NUMAAKNOWN; col++) printf("%3d",chitab[row][col]);
      printf("\n");
   }
#endif
}


/************************************************************************/
/*>static int FindChiIndex(char *resnam)
   -------------------------------------
*//**

   \param[in]     *resnam     Residue name
   \return                       Pointer into Chi table

   Find the index into the Chi table for a residue name
-  13.05.92 Original
*/
static int FindChiIndex(char *resnam)
{
   int ret = -1;
   
   if(!strncmp(resnam,"ALA",3)) ret = 0;
   if(!strncmp(resnam,"ARG",3)) ret = 1;
   if(!strncmp(resnam,"ASN",3)) ret = 2;
   if(!strncmp(resnam,"ASP",3)) ret = 3;
   if(!strncmp(resnam,"CYS",3)) ret = 4;
   if(!strncmp(resnam,"GLN",3)) ret = 5;
   if(!strncmp(resnam,"GLU",3)) ret = 6;
   if(!strncmp(resnam,"GLY",3)) ret = 7;
   if(!strncmp(resnam,"HIS",3)) ret = 8;
   if(!strncmp(resnam,"ILE",3)) ret = 9;
   if(!strncmp(resnam,"LEU",3)) ret = 10;
   if(!strncmp(resnam,"LYS",3)) ret = 11;
   if(!strncmp(resnam,"MET",3)) ret = 12;
   if(!strncmp(resnam,"PHE",3)) ret = 13;
   if(!strncmp(resnam,"PRO",3)) ret = 14;
   if(!strncmp(resnam,"SER",3)) ret = 15;
   if(!strncmp(resnam,"THR",3)) ret = 16;
   if(!strncmp(resnam,"TRP",3)) ret = 17;
   if(!strncmp(resnam,"TYR",3)) ret = 18;
   if(!strncmp(resnam,"VAL",3)) ret = 19;
   
   return(ret);
}


/************************************************************************/
/*>static PDB *FixTorsions(PDB *pdb, PDB *ResStart, PDB *NextRes, 
                           int **chitab)
   --------------------------------------------------------------
*//**

   \param[in,out] *pdb         Linked list to fix torsions
   \param[in]     *ResStart    Beginning of reference frag
   \param[in]     *NextRes     Start of next residue
   \param[in]     **chitab     Equivalent chis table
   \return                       Start of residue. Normally same as input
                              ResStart, but may differ as backbone
                              order will be corrected.
   
   Take the equivalent torsion angle between ResStart and NextRes and 
   apply them to pdb.
   
-  13.05.92 Original
-  21.06.93 Changed to use Array2D allocated chitab 
-  09.02.05 Chain name was getting set to last one in pdb
-  14.12.16 Added check on return from blCalcChi()
*/
static PDB *FixTorsions(PDB *pdb,      /* Linked list to fix torsions   */
                        PDB *ResStart, /* Beginning of reference frag   */
                        PDB *NextRes,  /* End of reference fragment     */
                        int **chitab)  /* Equivalent chis               */
{
   int   nchi,
         i,
         j;
   REAL  ParentChi;
   char  chain[8];
   PDB   *p;

   strcpy(chain, ResStart->chain);
   
   /* Correct the atom order of pdb and ResStart                        */
   ResStart = blShuffleBB(ResStart);

   /* 09.02.05 Fix the chain name                                       */
   for(p=ResStart; p!=NextRes; NEXT(p))
   {
      strcpy(p->chain, chain);
   }

   i = FindChiIndex(pdb->resnam);
   j = FindChiIndex(ResStart->resnam);

   /* If we haven't found one, we just return                           */
   if(i<0 || j<0) return(ResStart);
   
   /* look in the table for the number of chis                          */
   nchi = chitab[i][j];
   
   /* For each of the chis                                              */
   for(j=0;j<nchi;j++)
   {
      if((ParentChi = blCalcChi(ResStart, j)) < 9998.0)
      {
         blSetChi(pdb, NULL, ParentChi, j);
      }
   }
   
   return(ResStart);
}



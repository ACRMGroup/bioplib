/************************************************************************/
/**

   \file       deprecated.c
   
   \version    V1.0
   \date       31.07.14
   \brief      Source code for all deprecated functions.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2014
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

   Allows code utilising BiopLib to compile without modification to use 
   function names with the 'bl' prefix.(The 'bl' prefix was introduced in 
   July 2014.) 

   Allows use of deprecated functions that took a single character for the
   PDB chain identifier or insert value.

   Prototypes for the deprecated functions are in their original header 
   files:

-  pdb.h
-  general.h
-  BuffInp.h
-  ErrStack.h
-  MathUtil.h
-  WindIO.h
-  aalist.h
-  angle.h
-  array.h
-  cssr.h
-  fit.h
-  hbond.h
-  help.h
-  hpgl.h
-  matrix.h
-  parse.h
-  plotting.h
-  ps.h
-  safemem.h
-  seq.h



**************************************************************************

   Usage:
   ======

   Code using deprecated functions can be compiled without modification, 
   however, all deprecated functions use the DEPRECATED macro which prints
   a warning message to stderr when a deprecated function is called.

   The warning message can be suppressed as either an option at the 
   compilation of the library or by setting an environment variable at 
   runtime.

**************************************************************************

   Revision History:
   =================

-  V1.0  31.07.14 Original By: CTP

*************************************************************************/
/* Includes
*/
#include "deprecated.h"
#include "macros.h"

#include "pdb.h"
#include "general.h"
#include "BuffInp.h"
#include "ErrStack.h"
#include "MathUtil.h"
#include "WindIO.h"
#include "aalist.h"
#include "angle.h"
#include "array.h"
#include "cssr.h"
#include "fit.h"
#include "hbond.h"
#include "help.h"
#include "hpgl.h"
#include "matrix.h"
#include "parse.h"
#include "plotting.h"
#include "ps.h"
#include "safemem.h"
#include "seq.h"

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


/************************************************************************/
/*>PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
   ---------------------------------------------------------------------
   Finds a pointer to the start of a residue in a PDB linked list, but
   requires the residue is a HETATM record.
   Uses char for chain and insert.

-  26.10.11 Original   By: ACRM
-  24.02.14 Converted into wrapper for BiopFindHetatmResidue(). By: CTP
-  25.02.14 Added error message. By: CTP
-  07.05.14 Converted to macro. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/
PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
{
   char chain_a[2]  = " ",
        insert_a[2] = " ";

   DEPRECATED("FindHetatmResidue()","blFindHetatmResidue()");

   chain_a[0]  = chain;
   insert_a[0] = insert;

   return(blFindHetatmResidue(pdb, chain_a, resnum, insert_a));

}


/************************************************************************/
/*>PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
   ---------------------------------------------------------------
   Finds a pointer to the start of a residue in a PDB linked list.
   Uses char for string and insert.

-  06.02.96 Original   By: ACRM
-  24.02.14 Converted into wrapper for BiopFindResidue(). By: CTP
-  25.02.14 Added error message. By: CTP
-  07.05.14 Converted to macro. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/

PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
{
   char chain_a[2]  = " ",
        insert_a[2] = " ";

   DEPRECATED("FindResidue()","blFindResidue()");

   chain_a[0]  = chain;
   insert_a[0] = insert;

   return(blFindResidue(pdb, chain_a, resnum, insert_a));
}


/************************************************************************/
/*>BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, 
                    char stopinsert, char chain, int mode, 
                    PDB **pdb_start, PDB **pdb_stop)
   -------------------------------------------------------------
*//**

   \param[in]     *pdb        PDB linked list
   \param[in]     start       Resnum of start of zone
   \param[in]     startinsert Insert code for start of zone
   \param[in]     stop        Resnum of end of zone
   \param[in]     stopinsert  Insert code for end of zone
   \param[in]     chain       Chain name
   \param[in]     mode        ZONE_MODE_RESNUM:     Use PDB residue 
                                                    numbers/chain
                              ZONE_MODE_SEQUENTIAL: Use sequential 
                                                    numbering
   \param[out]    **pdb_start Start of zone
   \param[out]    **pdb_stop  End of zone
   \return                      OK?

   Finds pointers to the start and end of a zone in a PDB linked list. The
   end is the atom *after* the specified zone

-  30.09.92 Original
-  17.07.95 Chain name was being ignored in specs like L* (for whole
            of light chain)
-  18.08.95 Now handles inserts
-  31.07.95 Fixed bug when zone end==chain end
-  20.02.01 Changed to -999/-999 for beginning/end of chain rather than -1/-1
-  20.03.14 Function deprecated. Converted to wrapper for blFindZonePDB() 
            By: CTP
-  07.05.14 Converted to macro. By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/
BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, 
                 char stopinsert, char chain, int mode, PDB **pdb_start,
                 PDB **pdb_stop)
{
   char startinsert_a[2] = " ",
        stopinsert_a[2]  = " ",
        chain_a[2]       = " ";
        
   DEPRECATED("FindZonePDB()","blFindZonePDB()");
   
   strncpy(startinsert_a,&startinsert,1);
   strncpy(stopinsert_a ,&stopinsert ,1);
   strncpy(chain_a      ,&chain      ,1);
   
   return(blFindZonePDB(pdb,
                        start,
                        startinsert_a,
                        stop,
                        stopinsert_a,
                        chain_a,
                        mode,
                        pdb_start,
                        pdb_stop));
}


/************************************************************************/
/*>BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1,
                  int resnum2, char insert2)
   ----------------------------------------------------------
*//**

   \param[in]     *p         Pointer to a PDB record
   \param[in]     chain      Chain name
   \param[in]     resnum1    First residue
   \param[in]     insert1    First insert code
   \param[in]     resnum2    Second residue
   \param[in]     insert2    Second insert code
   \return                      Is p in the range specified?

   Checks that atom stored in PDB pointer p is within the specified 
   residue range.

   N.B. This assumes ASCII coding.

-  29.03.95 Original    By: ACRM
-  08.02.96 Insert residues inside a zone were not handled correctly!
-  18.06.96 Added to bioplib from QTree (was called InZone())
-  24.02.14 Converted into wrapper for BiopInPDBZone() By: CTP
-  25.02.14 Added error message. By: CTP
-  07.05.14 Converted to macro. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/
BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1,
                  int resnum2, char insert2)
{
   char chain_a[2]   = " ",
        insert1_a[2] = " ",
        insert2_a[2] = " ";
        
   DEPRECATED("InPDBZone()","blInPDBZone()");
   
   chain_a[0]   = chain;
   insert1_a[0] = insert1;
   insert2_a[0] = insert2;
   
   return(blInPDBZone(p,chain_a,resnum1,insert1_a,resnum2, insert2_a));
}

/************************************************************************/
/*>void WritePDB(FILE *fp, PDB *pdb)
   ---------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list by calls to WritePDBRecord()

-  08.03.89 Original
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 Uses NEXT macro; void type
-  08.07.93 Added insertion of TER cards
-  22.02.94 And a TER card at the end of the file
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  21.06.14 Function deprecated. Converted to wrapper for blWriteAsPDB().
            By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/
void WritePDB(FILE *fp, PDB *pdb)
{
   DEPRECATED("WritePDB","blWritePDB()");
   blWriteAsPDB(fp, pdb);
}

/************************************************************************/
/*>void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
   --------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes a PDB file including header and trailer information

-  30.05.02  Original   By: ACRM
-  21.06.14  Deprecated function. By CTP.
-  07.07.14  Use bl prefix for functions. Moved to deprecated.h By: CTP
-  31.07.14 Moved to deprecated.c and converted from macro to function. 
            By: CTP
*/
void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
{
   DEPRECATED("WriteWholePDB","blWriteWholePDB()");
   
   blWriteWholePDBHeader(fp, wpdb);
   blWritePDB(fp, wpdb->pdb);
   blWriteWholePDBTrailer(fp, wpdb);
}


/************************************************************************/
/*>void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
   --------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes the header of a PDB file 

-  30.05.02  Original   By: ACRM
-  21.06.14  Deprecated By: CTP
-  07.07.14  Moved to deprecated.h By: CTP
-  31.07.14  Moved to deprecated.c and converted from macro to function. 
             By: CTP
*/
void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
{
   DEPRECATED("WriteWholePDBHeader()","blWriteWholePDBHeader()");
   blWriteWholePDBHeader(fp, wpdb);
}


/************************************************************************/
/*>void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
   ---------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes the trailer of a PDB file 

-  30.05.02  Original   By: ACRM
-  21.06.14  Deprecated By: CTP
-  07.07.14  Moved to deprecated.h By: CTP
-  31.07.14  Moved to deprecated.c and converted from macro to function. 
             By: CTP
*/
void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
{
   DEPRECATED("WriteWholePDBTrailer()","blWriteWholePDBTrailer()");
   blWriteWholePDBTrailer(fp, wpdb);
}


/************************************************************************/
/*>char *GetPDBChainLabels(PDB *pdb)
   ---------------------------------
*//**

   \param[in]     *pdb      PDB linked list
   \return                     Allocated string containing chain labels
                             NULL if unable to allocate memory

   Scans a PDB linked list for chain names. Allocates memory for a 
   string containing these labels which is returned.

   N.B. You must free the allocated memory when you've finished with it!

-  25.07.95 Original    By: ACRM
-  25.03.14 Added deprecated message. By: CTP
-  07.05.14 Use DEPRECATED() macro. By: CTP
-  31.07.14 Moved to deprecated.c By: CTP
*/
char *GetPDBChainLabels(PDB *pdb)
{
   char *chains;
   int  nchains   = 0,
        maxchains = 16;
   PDB  *p;

   DEPRECATED("GetPDBChainLabels()","blGetPDBChainLabels()");
   
   /* Just return if linked list is NULL                                */
   if(pdb==NULL)
      return(NULL);

   /* Allocate a chunk for storing the chains                           */
   if((chains = (char *)malloc(maxchains * sizeof(char)))==NULL)
      return(NULL);

   /* Set up first chain label                                          */
   chains[nchains] = pdb->chain[0];

   /* Run through the linked list                                       */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(p->chain[0] != chains[nchains])
      {
         /* Increment chain count and reallocate memory if needed       */
         if(++nchains == maxchains)
         {
            maxchains += 16;
            if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
               return(NULL);
         }
         /* Store this new chain label                                  */
         chains[nchains] = p->chain[0];
      }
   }

   /* Increment chain count and reallocate memory if needed             */
   if(++nchains == maxchains)
   {
      maxchains += 16;
      if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
         return(NULL);
   }

   /* Terminate the chain list with a NUL character                     */
   chains[nchains] = '\0';

   return(chains);
}



/************************************************************************/
/* Renamed functions: pdb.h                                             */

PDB *ReadPDB(FILE *fp, int *natom)
{
   DEPRECATED("ReadPDB()","blReadPDB()");
   return(blReadPDB(fp, natom));
}

PDB *ReadPDBAll(FILE *fp, int *natom)
{
   DEPRECATED("ReadPDBAll()","blReadPDBAll()");
   return(blReadPDBAll(fp, natom));
}

PDB *ReadPDBAtoms(FILE *fp, int *natom)
{
   DEPRECATED("ReadPDBAtoms()","blReadPDBAtoms()");
   return(blReadPDBAtoms(fp, natom));
}

PDB *ReadPDBOccRank(FILE *fp, int *natom, int OccRank)
{
   DEPRECATED("ReadPDBOccRank()","blReadPDBOccRank()");
   return(blReadPDBOccRank(fp, natom, OccRank));
}

PDB *ReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank)
{
   DEPRECATED("ReadPDBAtomsOccRank()","blReadPDBAtomsOccRank()");
   return(blReadPDBAtomsOccRank(fp, natom, OccRank));
}

PDB *doReadPDB(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, int ModelNum)
{
   DEPRECATED("doReadPDB()","blDoReadPDB()");
   return(blDoReadPDB(fp, natom, AllAtoms, OccRank, ModelNum));
}

PDB *doReadPDBML(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, int ModelNum)
{
   DEPRECATED("doReadPDBML()","blDoReadPDBML()");
   return(blDoReadPDBML(fp, natom, AllAtoms, OccRank, ModelNum));
}
                
BOOL CheckFileFormatPDBML(FILE *fp)
{
   DEPRECATED("CheckFileFormatPDBML()","blCheckFileFormatPDBML()");
   return(blCheckFileFormatPDBML(fp));
}

void WriteAsPDB(FILE *fp, PDB  *pdb)
{
   DEPRECATED("WriteAsPDB()","blWriteAsPDB()");
   blWriteAsPDB(fp, pdb);
}

void WriteAsPDBML(FILE *fp, PDB  *pdb)
{
   DEPRECATED("WriteAsPDBML()","blWriteAsPDBML()");
   blWriteAsPDBML(fp, pdb);
}

BOOL FormatCheckWritePDB(PDB *pdb)
{
   DEPRECATED("FormatCheckWritePDB()","blFormatCheckWritePDB()");
   return(blFormatCheckWritePDB(pdb));
}

void WritePDBRecord(FILE *fp, PDB *pdb)
{
   DEPRECATED("WritePDBRecord()","blWritePDBRecord()");
   blWritePDBRecord(fp, pdb);
}

void WritePDBRecordAtnam(FILE *fp, PDB  *pdb)
{
   DEPRECATED("WritePDBRecordAtnam()","blWritePDBRecordAtnam()");
   blWritePDBRecordAtnam(fp, pdb);
}

void WriteGromosPDB(FILE *fp, PDB *pdb)
{
   DEPRECATED("WriteGromosPDB()","blWriteGromosPDB()");
   blWriteGromosPDB(fp, pdb);
}

void WriteGromosPDBRecord(FILE *fp, PDB *pdb)
{
   DEPRECATED("WriteGromosPDBRecord()","blWriteGromosPDBRecord()");
   blWriteGromosPDBRecord(fp, pdb);
}

void GetCofGPDB(PDB   *pdb, VEC3F *cg)
{
   DEPRECATED("GetCofGPDB()","blGetCofGPDB()");
   blGetCofGPDB(pdb, cg);
}

void GetCofGPDBRange(PDB *start, PDB *stop, VEC3F *cg)
{
   DEPRECATED("GetCofGPDBRange()","blGetCofGPDBRange()");
   blGetCofGPDBRange(start, stop, cg);
}

void GetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg)
{
   DEPRECATED("GetCofGPDBSCRange()","blGetCofGPDBSCRange()");
   blGetCofGPDBSCRange(start, stop, cg);
}

void OriginPDB(PDB *pdb)
{
   DEPRECATED("OriginPDB()","blOriginPDB()");
   blOriginPDB(pdb);
}

void RotatePDB(PDB  *pdb, REAL rm[3][3])
{
   DEPRECATED("RotatePDB()","blRotatePDB()");
   blRotatePDB(pdb, rm);
}

void TranslatePDB(PDB   *pdb, VEC3F tvect)
{
   DEPRECATED("TranslatePDB()","blTranslatePDB()");
   blTranslatePDB(pdb, tvect);
}

BOOL FitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   DEPRECATED("FitPDB()","blFitPDB()");
   return(blFitPDB(ref_pdb, fit_pdb, rm));
}

BOOL FitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   DEPRECATED("FitCaPDB()","blFitCaPDB()");
   return(blFitCaPDB(ref_pdb, fit_pdb, rm));
}

BOOL FitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   DEPRECATED("FitNCaCPDB()","blFitNCaCPDB()");
   return(blFitNCaCPDB(ref_pdb, fit_pdb, rm));
}

BOOL FitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   DEPRECATED("FitCaCbPDB()","blFitCaCbPDB()");
   return(blFitCaCbPDB(ref_pdb, fit_pdb, rm));
}

REAL CalcRMSPDB(PDB *pdb1, PDB *pdb2)
{
   DEPRECATED("CalcRMSPDB()","blCalcRMSPDB()");
   return(blCalcRMSPDB(pdb1, pdb2));
}

int GetPDBCoor(PDB *pdb, COOR **coor)
{
   DEPRECATED("GetPDBCoor()","blGetPDBCoor()");
   return(blGetPDBCoor(pdb, coor));
}
      
int HAddPDB(FILE *fp, PDB *pdb)
{
   DEPRECATED("HAddPDB()","blHAddPDB()");
   return(blHAddPDB(fp, pdb));
}

int ReadPGP(FILE *fp)
{
   DEPRECATED("ReadPGP()","blReadPGP()");
   return(blReadPGP(fp));
}

FILE *OpenPGPFile(char *pgpfile, BOOL AllHyd)
{
   DEPRECATED("OpenPGPFile()","blOpenPGPFile()");
   return(blOpenPGPFile(pgpfile, AllHyd));
}

PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom)
{
   DEPRECATED("SelectAtomsPDB()","blSelectAtomsPDB()");
   return(blSelectAtomsPDB(pdbin, nsel, sel, natom));
}

PDB *StripHPDB(PDB *pdbin, int *natom)
{
   DEPRECATED("StripHPDB()","blStripHPDB()");
   return(blStripHPDB(pdbin, natom));
}

SECSTRUC *ReadSecPDB(FILE *fp, int *nsec)
{
   DEPRECATED("ReadSecPDB()","blReadSecPDB()");
   return(blReadSecPDB(fp, nsec));
}

void RenumAtomsPDB(PDB *pdb)
{
   DEPRECATED("RenumAtomsPDB()","blRenumAtomsPDB()");
   blRenumAtomsPDB(pdb);
}

PDB *FindEndPDB(PDB *start)
{
   DEPRECATED("FindEndPDB()","blFindEndPDB()");
   return(blFindEndPDB(start));
}

PDB *FixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum)
{
   DEPRECATED("FixOrderPDB()","blFixOrderPDB()");
   return(blFixOrderPDB(pdb, Pad, Renum));
}

PDB *ShuffleResPDB(PDB *start, PDB *end, BOOL Pad)
{
   DEPRECATED("ShuffleResPDB()","blShuffleResPDB()");
   return(blShuffleResPDB(start, end, Pad));
}

BOOL GetAtomTypes(char *resnam, char **AtomTypes)
{
   DEPRECATED("GetAtomTypes()","blGetAtomTypes()");
   return(blGetAtomTypes(resnam, AtomTypes));
}

PDB *KillPDB(PDB *pdb, PDB *prev)
{
   DEPRECATED("KillPDB()","blKillPDB()");
   return(blKillPDB(pdb, prev));
}

void CopyPDB(PDB *out, PDB *in)
{
   DEPRECATED("CopyPDB()","blCopyPDB()");
   blCopyPDB(out, in);
}

BOOL MovePDB(PDB *move, PDB **from, PDB **to)
{
   DEPRECATED("MovePDB()","blMovePDB()");
   return(blMovePDB(move, from, to));
}

PDB *AppendPDB(PDB *first, PDB *second)
{
   DEPRECATED("AppendPDB()","blAppendPDB()");
   return(blAppendPDB(first, second));
}

PDB *ShuffleBB(PDB *pdb)
{
   DEPRECATED("ShuffleBB()","blShuffleBB()");
   return(blShuffleBB(pdb));
}

REAL CalcChi(PDB *pdb, int type)
{
   DEPRECATED("CalcChi()","blCalcChi()");
   return(blCalcChi(pdb, type));
}

PDB *GetPDBByN(PDB *pdb, int n)
{
   DEPRECATED("GetPDBByN()","blGetPDBByN()");
   return(blGetPDBByN(pdb, n));
}

void SetChi(PDB *pdb, PDB *next, REAL chi, int type)
{
   DEPRECATED("SetChi()","blSetChi()");
   blSetChi(pdb, next, chi, type);
}

BOOL KillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB)
{
   DEPRECATED("KillSidechain()","blKillSidechain()");
   return(blKillSidechain(ResStart, NextRes, doCB));
}

void SetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum, char *insert, char *chain)
{
   DEPRECATED("SetResnam()","blSetResnam()");
   blSetResnam(ResStart, NextRes, resnam, resnum, insert, chain);
}

void ApplyMatrixPDB(PDB *pdb, REAL matrix[3][3])
{
   DEPRECATED("ApplyMatrixPDB()","blApplyMatrixPDB()");
   blApplyMatrixPDB(pdb, matrix);
}

BOOL GetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor, int *StrucType)
{
   DEPRECATED("GetResolPDB()","blGetResolPDB()");
   return(blGetResolPDB(fp, resolution, RFactor, StrucType));
}

BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR, int *StrucType)
{
   DEPRECATED("GetExptl()","blGetExptl()");
   return(blGetExptl(fp, resolution, RFactor, FreeR, StrucType));
}

BOOL GetExptlOld(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR, int *StrucType)
{
   DEPRECATED("GetExptlOld()","blGetExptlOld()");
   return(blGetExptlOld(fp, resolution, RFactor, FreeR, StrucType));
}

char *ReportStructureType(int type)
{
   DEPRECATED("ReportStructureType()","blReportStructureType()");
   return(blReportStructureType(type));
}

PDB **IndexPDB(PDB *pdb, int *natom)
{
   DEPRECATED("IndexPDB()","blIndexPDB()");
   return(blIndexPDB(pdb, natom));
}

DISULPHIDE *ReadDisulphidesPDB(FILE *fp, BOOL *error)
{
   DEPRECATED("ReadDisulphidesPDB()","blReadDisulphidesPDB()");
   return(blReadDisulphidesPDB(fp, error));
}

BOOL ParseResSpec(char *spec, char *chain, int *resnum, char *insert)
{
   DEPRECATED("ParseResSpec()","blParseResSpec()");
   return(blParseResSpec(spec, chain, resnum, insert));
}

BOOL ParseResSpecNoUpper(char *spec, char *chain, int *resnum, char *insert)
{
   DEPRECATED("ParseResSpecNoUpper()","blParseResSpecNoUpper()");
   return(blParseResSpecNoUpper(spec, chain, resnum, insert));
}

BOOL DoParseResSpec(char *spec, char *chain, int *resnum, char *insert, BOOL uppercaseresspec)
{
   DEPRECATED("DoParseResSpec()","blDoParseResSpec()");
   return(blDoParseResSpec(spec, chain, resnum, insert, uppercaseresspec));
}

BOOL RepSChain(PDB *pdb, char *sequence, char *ChiTable, char *RefCoords)
{
   DEPRECATED("RepSChain()","blRepSChain()");
   return(blRepSChain(pdb, sequence, ChiTable, RefCoords));
}

PDB *FindNextChainPDB(PDB *pdb)
{
   DEPRECATED("FindNextChainPDB()","blFindNextChainPDB()");
   return(blFindNextChainPDB(pdb));
}

BOOL FixCterPDB(PDB *pdb, int style)
{
   DEPRECATED("FixCterPDB()","blFixCterPDB()");
   return(blFixCterPDB(pdb, style));
}

BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p)
{
   DEPRECATED("CalcCterCoords()","blCalcCterCoords()");
   return(blCalcCterCoords(p, ca_p, c_p, o_p));
}

int CalcTetraHCoords(PDB *nter, COOR *coor)
{
   DEPRECATED("CalcTetraHCoords()","blCalcTetraHCoords()");
   return(blCalcTetraHCoords(nter, coor));
}

int AddNTerHs(PDB **ppdb, BOOL Charmm)
{
   DEPRECATED("AddNTerHs()","blAddNTerHs()");
   return(blAddNTerHs(ppdb, Charmm));
}

char *FNam2PDB(char *filename)
{
   DEPRECATED("FNam2PDB()","blFNam2PDB()");
   return(blFNam2PDB(filename));
}

PDB *TermPDB(PDB *pdb, int length)
{
   DEPRECATED("TermPDB()","blTermPDB()");
   return(blTermPDB(pdb, length));
}

PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec)
{
   DEPRECATED("FindHetatmResidueSpec()","blFindHetatmResidueSpec()");
   return(blFindHetatmResidueSpec(pdb, resspec));
}

PDB *FindResidueSpec(PDB *pdb, char *resspec)
{
   DEPRECATED("FindResidueSpec()","blFindResidueSpec()");
   return(blFindResidueSpec(pdb, resspec));
}

PDB *FindNextResidue(PDB *pdb)
{
   DEPRECATED("FindNextResidue()","blFindNextResidue()");
   return(blFindNextResidue(pdb));
}

PDB *DupePDB(PDB *in)
{
   DEPRECATED("DupePDB()","blDupePDB()");
   return(blDupePDB(in));
}

BOOL CopyPDBCoords(PDB *out, PDB *in)
{
   DEPRECATED("CopyPDBCoords()","blCopyPDBCoords()");
   return(blCopyPDBCoords(out, in));
}

void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans)
{
   DEPRECATED("CalcCellTrans()","blCalcCellTrans()");
   blCalcCellTrans(UnitCell, CellAngles, xtrans, ytrans, ztrans);
}

int GetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles, char *spacegroup, REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
{
   DEPRECATED("GetCrystPDB()","blGetCrystPDB()");
   return(blGetCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix));
}

void WriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles, char *spacegroup, REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
{
   DEPRECATED("WriteCrystPDB()","blWriteCrystPDB()");
   blWriteCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix);
}

PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1, char *chain2, int resnum2, char *insert2)
{
   DEPRECATED("ExtractZonePDB()","blExtractZonePDB()");
   return(blExtractZonePDB(inpdb, chain1, resnum1, insert1, chain2, resnum2, insert2));
}

PDB *FindAtomInRes(PDB *pdb, char *atnam)
{
   DEPRECATED("FindAtomInRes()","blFindAtomInRes()");
   return(blFindAtomInRes(pdb, atnam));
}

BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2)
{
   DEPRECATED("InPDBZoneSpec()","blInPDBZoneSpec()");
   return(blInPDBZoneSpec(p, resspec1, resspec2));
}

BOOL AtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   DEPRECATED("AtomNameMatch()","blAtomNameMatch()");
   return(blAtomNameMatch(atnam, spec, ErrorWarn));
}

BOOL AtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn)
{
   DEPRECATED("AtomNameRawMatch()","blAtomNameRawMatch()");
   return(blAtomNameRawMatch(atnam, spec, ErrorWarn));
}

BOOL LegalAtomSpec(char *spec)
{
   DEPRECATED("LegalAtomSpec()","blLegalAtomSpec()");
   return(blLegalAtomSpec(spec));
}

BOOL RepOneSChain(PDB *pdb, char *ResSpec, char aa, char *ChiTable, char *RefCoords)
{
   DEPRECATED("RepOneSChain()","blRepOneSChain()");
   return(blRepOneSChain(pdb, ResSpec, aa, ChiTable, RefCoords));
}

void EndRepSChain(void)
{
   DEPRECATED("EndRepSChain()","blEndRepSChain()");
   blEndRepSChain();
}

char **ReadSeqresPDB(FILE *fp, int *nchains)
{
   DEPRECATED("ReadSeqresPDB()","blReadSeqresPDB()");
   return(blReadSeqresPDB(fp, nchains));
}

PDB *SelectCaPDB(PDB *pdb)
{
   DEPRECATED("SelectCaPDB()","blSelectCaPDB()");
   return(blSelectCaPDB(pdb));
}

char *FixAtomName(char *name, REAL occup)
{
   DEPRECATED("FixAtomName()","blFixAtomName()");
   return(blFixAtomName(name, occup));
}

void FreeWholePDB(WHOLEPDB *wpdb)
{
   DEPRECATED("FreeWholePDB()","blFreeWholePDB()");
   blFreeWholePDB(wpdb);
}

WHOLEPDB *ReadWholePDB(FILE *fpin)
{
   DEPRECATED("ReadWholePDB()","blReadWholePDB()");
   return(blReadWholePDB(fpin));
}

WHOLEPDB *ReadWholePDBAtoms(FILE *fpin)
{
   DEPRECATED("ReadWholePDBAtoms()","blReadWholePDBAtoms()");
   return(blReadWholePDBAtoms(fpin));
}

BOOL AddCBtoGly(PDB *pdb)
{
   DEPRECATED("AddCBtoGly()","blAddCBtoGly()");
   return(blAddCBtoGly(pdb));
}

BOOL AddCBtoAllGly(PDB *pdb)
{
   DEPRECATED("AddCBtoAllGly()","blAddCBtoAllGly()");
   return(blAddCBtoAllGly(pdb));
}

PDB *StripGlyCB(PDB *pdb)
{
   DEPRECATED("StripGlyCB()","blStripGlyCB()");
   return(blStripGlyCB(pdb));
}

PDB *RemoveAlternates(PDB *pdb)
{
   DEPRECATED("RemoveAlternates()","blRemoveAlternates()");
   return(blRemoveAlternates(pdb));
}

PDB *BuildAtomNeighbourPDBList(PDB *pdb, PDB *pRes, REAL NeighbDist)
{
   DEPRECATED("BuildAtomNeighbourPDBList()","blBuildAtomNeighbourPDBList()");
   return(blBuildAtomNeighbourPDBList(pdb, pRes, NeighbDist));
}

PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern)
{
   DEPRECATED("FindAtomWildcardInRes()","blFindAtomWildcardInRes()");
   return(blFindAtomWildcardInRes(pdb, pattern));
}

PDB *DupeResiduePDB(PDB *in)
{
   DEPRECATED("DupeResiduePDB()","blDupeResiduePDB()");
   return(blDupeResiduePDB(in));
}

PDB *StripWatersPDB(PDB *pdbin, int *natom)
{
   DEPRECATED("StripWatersPDB()","blStripWatersPDB()");
   return(blStripWatersPDB(pdbin, natom));
}

PDBSTRUCT *AllocPDBStructure(PDB *pdb)
{
   DEPRECATED("AllocPDBStructure()","blAllocPDBStructure()");
   return(blAllocPDBStructure(pdb));
}

PDB *FindNextChain(PDB *pdb)
{
   DEPRECATED("FindNextChain()","blFindNextChain()");
   return(blFindNextChain(pdb));
}

void FreePDBStructure(PDBSTRUCT *pdbstruct)
{
   DEPRECATED("FreePDBStructure()","blFreePDBStructure()");
   blFreePDBStructure(pdbstruct);
}

void SetElementSymbolFromAtomName(char *element, char *atom_name)
{
   DEPRECATED("SetElementSymbolFromAtomName()","blSetElementSymbolFromAtomName()");
   blSetElementSymbolFromAtomName(element, atom_name);
}



/************************************************************************/
/* Renamed functions: general.h                                         */


void StringToLower(char *string1, char *string2)
{
   DEPRECATED("StringToLower()","blStringToLower()");
   blStringToLower(string1, string2);
}

void StringToUpper(char *string1, char *string2)
{
   DEPRECATED("StringToUpper()","blStringToUpper()");
   blStringToUpper(string1, string2);
}

char *KillLeadSpaces(char *string)
{
   DEPRECATED("KillLeadSpaces()","blKillLeadSpaces()");
   return(blKillLeadSpaces(string));
}

void KillLine(FILE *fp)
{
   DEPRECATED("KillLine()","blKillLine()");
   blKillLine(fp);
}

void SetExtn(char *File, char *Ext)
{
   DEPRECATED("SetExtn()","blSetExtn()");
   blSetExtn(File, Ext);
}

int chindex(char *string, char ch)
{
   DEPRECATED("chindex()","blChindex()");
   return(blChindex(string, ch));
}

void Word(char *string1, char *string2)
{
   DEPRECATED("Word()","blWord()");
   blWord(string1, string2);
}

void WordN(char *string1, char *string2, int  MaxChar)
{
   DEPRECATED("WordN()","blWordN()");
   blWordN(string1, string2, MaxChar);
}

void padterm(char *string, int length)
{
   DEPRECATED("padterm()","blPadterm()");
   blPadterm(string, length);
}

void padchar(char *string, int length, char ch)
{
   DEPRECATED("padchar()","blPadchar()");
   blPadchar(string, length, ch);
}

BOOL CheckExtn(char *string, char *ext)
{
   DEPRECATED("CheckExtn()","blCheckExtn()");
   return(blCheckExtn(string, ext));
}

char *ftostr(char *str, int maxlen, REAL x, int precision)
{
   DEPRECATED("ftostr()","blFtostr()");
   return(blFtostr(str, maxlen, x, precision));
}

void GetFilestem(char *filename, char *stem)
{
   DEPRECATED("GetFilestem()","blGetFilestem()");
   blGetFilestem(filename, stem);
}

int upstrcmp(char *word1, char *word2)
{
   DEPRECATED("upstrcmp()","blUpstrcmp()");
   return(blUpstrcmp(word1, word2));
}

int upstrncmp(char *word1, char *word2, int ncomp)
{
   DEPRECATED("upstrncmp()","blUpstrncmp()");
   return(blUpstrncmp(word1, word2, ncomp));
}

char *GetWord(char *buffer, char *word, int maxsize)
{
   DEPRECATED("GetWord()","blGetWord()");
   return(blGetWord(buffer, word, maxsize));
}

BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
{
   DEPRECATED("OpenStdFiles()","blOpenStdFiles()");
   return(blOpenStdFiles(infile, outfile, in, out));
}

FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
{
   DEPRECATED("OpenFile()","blOpenFile()");
   return(blOpenFile(filename, envvar, mode, noenv));
}

int countchar(char *string, char ch)
{
   DEPRECATED("countchar()","blCountchar()");
   return(blCountchar(string, ch));
}

char *fgetsany(FILE *fp)
{
   DEPRECATED("fgetsany()","blFgetsany()");
   return(blFgetsany(fp));
}

char *strcatalloc(char *instr, char *catstr)
{
   DEPRECATED("strcatalloc()","blStrcatalloc()");
   return(blStrcatalloc(instr, catstr));
}

STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
{
   DEPRECATED("StoreString()","blStoreString()");
   return(blStoreString(StringList, string));
}

BOOL InStringList(STRINGLIST *StringList, char *string)
{
   DEPRECATED("InStringList()","blInStringList()");
   return(blInStringList(StringList, string));
}

void FreeStringList(STRINGLIST *StringList)
{
   DEPRECATED("FreeStringList()","blFreeStringList()");
   blFreeStringList(StringList);
}

char *QueryStrStr(char *string, char *substring)
{
   DEPRECATED("QueryStrStr()","blQueryStrStr()");
   return(blQueryStrStr(string, substring));
}

void IndexReal(REAL *arrin, int *indx, int n)
{
   DEPRECATED("IndexReal()","blIndexReal()");
   blIndexReal(arrin, indx, n);
}

FILE *OpenOrPipe(char *filename)
{
   DEPRECATED("OpenOrPipe()","blOpenOrPipe()");
   return(blOpenOrPipe(filename));
}

int CloseOrPipe(FILE *fp)
{
   DEPRECATED("CloseOrPipe()","blCloseOrPipe()");
   return(blCloseOrPipe(fp));
}

BOOL WrapString(char *in, char *out, int maxlen)
{
   DEPRECATED("WrapString()","blWrapString()");
   return(blWrapString(in, out, maxlen));
}

BOOL WrapPrint(FILE *out, char *string)
{
   DEPRECATED("WrapPrint()","blWrapPrint()");
   return(blWrapPrint(out, string));
}

void RightJustify(char *string)
{
   DEPRECATED("RightJustify()","blRightJustify()");
   blRightJustify(string);
}

char *GetWordNC(char *buffer, char *word, int maxlen)
{
   DEPRECATED("GetWordNC()","blGetWordNC()");
   return(blGetWordNC(buffer, word, maxlen));
}

void getfield(char *buffer, int start, int width, char *str)
{
   DEPRECATED("getfield()","blGetfield()");
   blGetfield(buffer, start, width, str);
}



/************************************************************************/
/* Renamed functions: BuffInp.h                                         */

INBUFFER *OpenBufferedFile(char *filename, int maxstr)
{
   DEPRECATED("OpenBufferedFile()","blOpenBufferedFile()");
   return(blOpenBufferedFile(filename, maxstr));
}

BOOL ReadBufferedFile(INBUFFER *bfp, char *string, int length)
{
   DEPRECATED("ReadBufferedFile()","blReadBufferedFile()");
   return(blReadBufferedFile(bfp, string, length));
}

BOOL ProbeBufferedFile(INBUFFER *bfp, char *string, int length)
{
   DEPRECATED("ProbeBufferedFile()","blProbeBufferedFile()");
   return(blProbeBufferedFile(bfp, string, length));
}



/************************************************************************/
/* Renamed functions: ErrStack.h                                        */

void StoreError(char *routine, char *error)
{
   DEPRECATED("StoreError()","blStoreError()");
   blStoreError(routine, error);
}

void ShowErrors(void *PrintRoutine(char *), BOOL Trace)
{
   DEPRECATED("ShowErrors()","blShowErrors()");
   blShowErrors(PrintRoutine, Trace);                     /* CHECK THIS */
}



/************************************************************************/
/* Renamed functions: MathUtil.h                                        */

void CalcSD(REAL val, int action, REAL *mean, REAL *SD)
{
   DEPRECATED("CalcSD()","blCalcSD()");
   blCalcSD(val, action, mean, SD);
}

void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, int *NValues, REAL *mean, REAL *SD)
{
   DEPRECATED("CalcExtSD()","blCalcExtSD()");
   blCalcExtSD(val, action, Sx, SxSq, NValues, mean, SD);
}

REAL pearson(REAL *x, REAL *y, int NItem)
{
   DEPRECATED("pearson()","blPearson()");
   return(blPearson(x, y, NItem));
}

REAL pearson1(REAL *x, REAL *y, int NItem)
{
   DEPRECATED("pearson1()","blPearson1()");
   return(blPearson1(x, y, NItem));
}


void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("CrossProd3()","blCrossProd3()");
   blCrossProd3(Out, In1, In2);
}

void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("VecSub3()","blVecSub3()");
   blVecSub3(Out, In1, In2);
}

void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("VecAdd3()","blVecAdd3()");
   blVecAdd3(Out, In1, In2);
}

REAL VecLen3(VEC3F Vec)
{
   DEPRECATED("VecLen3()","blVecLen3()");
   return(blVecLen3(Vec));
}

REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2)
{
   DEPRECATED("DistPtVect()","blDistPtVect()");
   return(blDistPtVect(Point, End1, End2));
}


REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac)
{
   DEPRECATED("PointLineDistance()","blPointLineDistance()");
   return(blPointLineDistance(Px,  Py,  Pz,
                              P1x, P1y, P1z,
                              P2x, P2y, P2z,
                              Rx,  Ry,  Rz,  frac));
}

ULONG factorial(int n)
{
   DEPRECATED("factorial()","blFactorial()");
   return(blFactorial(n));
}

ULONG factdiv(int n1, int n2)
{
   DEPRECATED("factdiv()","blFactdiv()");
   return(blFactdiv(n1, n2));
}

ULONG NPerm(int n, int r)
{
   DEPRECATED("NPerm()","blNPerm()");
   return(blNPerm(n, r));
}

ULONG NComb(int n, int r)
{
   DEPRECATED("NComb()","blNComb()");
   return(blNComb(n, r));
}



/************************************************************************/
/* Renamed functions: WindIO.h                                        */

void screen(char *string)
{
   DEPRECATED("Screen()","blScreen()");
   blScreen(string);
}

void prompt(char *string)
{
   DEPRECATED("Prompt()","blPrompt()");
   blPrompt(string);
}

void RePrompt(void)
{
   DEPRECATED("RePrompt()","blRePrompt()");
   blRePrompt();
}

void GetKybdString(char *string, int maxlen)
{
   DEPRECATED("GetKybdString()","blGetKybdString()");
   blGetKybdString(string, maxlen);
}

void PagingOn(void)
{
   DEPRECATED("PagingOn()","blPagingOn()");
   blPagingOn();
}

void PagingOff(void)
{
   DEPRECATED("PagingOff()","blPagingOff()");
   blPagingOff();
}

void WindowMode(BOOL mode)
{
   DEPRECATED("WindowMode()","blWindowMode()");
   blWindowMode(mode);
}

void WindowInteractive(BOOL mode)
{
   DEPRECATED("WindowInteractive()","blWindowInteractive()");
   blWindowInteractive(mode);
}

int YorN(char deflt)
{
   DEPRECATED("YorN()","blYorN()");
   return(blYorN(deflt));
}


/************************************************************************/
/* Renamed functions: aalist.h                                        */

AA *InsertNextResiduesInAAList(AA *a, char res, int nres)
{
   DEPRECATED("InsertNextResiduesInAAList()","blInsertNextResiduesInAAList()");
   return(blInsertNextResiduesInAAList(a, res, nres));
}

AA *InsertNextResidueInAAList(AA *a, char res)
{
   DEPRECATED("InsertNextResidueInAAList()","blInsertNextResidueInAAList()");
   return(blInsertNextResidueInAAList(a, res));
}

char *BuildSeqFromAAList(AA *aa)
{
   DEPRECATED("BuildSeqFromAAList()","blBuildSeqFromAAList()");
   return(blBuildSeqFromAAList(aa));
}

AA *InsertResidueInAAListAt(AA *aa, char res, int pos)
{
   DEPRECATED("InsertResidueInAAListAt()","blInsertResidueInAAListAt()");
   return(blInsertResidueInAAListAt(aa, res, pos));
}

AA *InsertResiduesInAAListAt(AA *aa, char res, int nres, int pos)
{
   DEPRECATED("InsertResiduesInAAListAt()","blInsertResiduesInAAListAt()");
   return(blInsertResiduesInAAListAt(aa, res, nres, pos));
}

AA *BuildAAList(char *seq)
{
   DEPRECATED("BuildAAList()","blBuildAAList()");
   return(blBuildAAList(seq));
}

int FindAAListOffsetByResnum(AA *aa, int resnum)
{
   DEPRECATED("FindAAListOffsetByResnum()","blFindAAListOffsetByResnum()");
   return(blFindAAListOffsetByResnum(aa, resnum));
}

AA *FindAAListItemByResnum(AA *aa, int resnum)
{
   DEPRECATED("FindAAListItemByResnum()","blFindAAListItemByResnum()");
   return(blFindAAListItemByResnum(aa, resnum));
}

void SetAAListFlagByResnum(AA *aa, int resnum)
{
   DEPRECATED("SetAAListFlagByResnum()","blSetAAListFlagByResnum()");
   blSetAAListFlagByResnum(aa, resnum);
}

char *BuildFlagSeqFromAAList(AA *aa, char ch)
{
   DEPRECATED("BuildFlagSeqFromAAList()","blBuildFlagSeqFromAAList()");
   return(blBuildFlagSeqFromAAList(aa, ch));
}

int GetAAListLen(AA *aa)
{
   DEPRECATED("GetAAListLen()","blGetAAListLen()");
   return(blGetAAListLen(aa));
}




/************************************************************************/
/* Renamed functions: angle.h                                        */

REAL angle(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj, REAL xk, REAL yk, REAL zk)
{
   DEPRECATED("Angle()","blAngle()");
   return(blAngle(xi, yi, zi, xj, yj, zj, xk, yk, zk));
}

REAL phi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj, REAL xk, REAL yk, REAL zk, REAL xl, REAL yl, REAL zl)
{
   DEPRECATED("Phi()","blPhi()");
   return(blPhi(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl));
}

REAL simpleangle(REAL ang)
{
   DEPRECATED("Simpleangle()","blSimpleangle()");
   return(blSimpleangle(ang));
}

REAL TrueAngle(REAL opp, REAL adj)
{
   DEPRECATED("TrueAngle()","blTrueAngle()");
   return(blTrueAngle(opp, adj));
}

BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, REAL bond, REAL theta, REAL torsion, VEC3F *coords)
{
   DEPRECATED("TorToCoor()","blTorToCoor()");
   return(blTorToCoor(ant1, ant2, ant3, bond, theta, torsion, coords));
}




/************************************************************************/
/* Renamed functions: array.h                                        */

char **Array2D(int size, int dim1, int dim2)
{
   DEPRECATED("Array2D()","blArray2D()");
   return(blArray2D(size, dim1, dim2));
}

void FreeArray2D(char **array, int dim1, int dim2)
{
   DEPRECATED("FreeArray2D()","blFreeArray2D()");
   blFreeArray2D(array, dim1, dim2);
}


char ***Array3D(int size, int dim1, int dim2, int dim3)
{
   DEPRECATED("Array3D()","blArray3D()");
   return(blArray3D(size, dim1, dim2, dim3));
}

void FreeArray3D(char ***array, int dim1, int dim2, int dim3)
{
   DEPRECATED("FreeArray3D()","blFreeArray3D()");
   blFreeArray3D(array, dim1, dim2, dim3);
}




/************************************************************************/
/* Renamed functions: cssr.h                                        */

CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title)
{
   DEPRECATED("ReadCSSR()","blReadCSSR()");
   return(blReadCSSR(fp, natom, name, title));
}

PDB *ReadCSSRasPDB(FILE *fp, int *natom)
{
   DEPRECATED("ReadCSSRasPDB()","blReadCSSRasPDB()");
   return(blReadCSSRasPDB(fp, natom));
}

void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta, REAL gamma)
{
   DEPRECATED("NormaliseCSSR()","blNormaliseCSSR()");
   blNormaliseCSSR(cssr, cell, alpha, beta, gamma);
}

void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta, REAL gamma)
{
   DEPRECATED("NormalisePDB()","blNormalisePDB()");
   blNormalisePDB(pdb, cell, alpha, beta, gamma);
}

void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma, REAL amatrx[3][3], int isw, int ncode)
{
   DEPRECATED("Ortho()","blOrtho()");
   blOrtho(cell, alpha, beta, gamma, amatrx, isw, ncode);
}

void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title)
{
   DEPRECATED("WriteCSSR()","blWriteCSSR()");
   blWriteCSSR(fp, cssr, name, title);
}




/************************************************************************/
/* Renamed functions: fit.h                                        */

BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, BOOL column)
{
   DEPRECATED("matfit()","blMatfit()");
   return(blMatfit(x1, x2, rm, n, wt1, column));
}



/************************************************************************/
/* Renamed functions: hbond.h                                        */

int  IsHBonded(PDB *res1, PDB *res2, int type);
BOOL ValidHBond(PDB *AtomH, PDB *AtomD, PDB *AtomA, PDB *AtomP)
{
   DEPRECATED("ValidHBond()","blValidHBond()");
   return(blValidHBond(AtomH, AtomD, AtomA, AtomP));
}

int IsMCDonorHBonded(PDB *res1, PDB *res2, int type)
{
   DEPRECATED("IsMCDonorHBonded()","blIsMCDonorHBonded()");
   return(blIsMCDonorHBonded(res1, res2, type));
}

int IsMCAcceptorHBonded(PDB *res1, PDB *res2, int type)
{
   DEPRECATED("IsMCAcceptorHBonded()","blIsMCAcceptorHBonded()");
   return(blIsMCAcceptorHBonded(res1, res2, type));
}



/************************************************************************/
/* Renamed functions: help.h                                        */

void Help(char *string, char *HelpFile)
{
   DEPRECATED("Help()","blHelp()");
   blHelp(string, HelpFile);
}

void DoHelp(char *string, char *HelpFile)
{
   DEPRECATED("DoHelp()","blDoHelp()");
   blDoHelp(string, HelpFile);
}



/************************************************************************/
/* Renamed functions: hpgl.h                                        */

BOOL HPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin)
{
   DEPRECATED("HPGLInit()","blHPGLInit()");
   return(blHPGLInit(filename, AltFont, xmargin, ymargin));
}

void HPGLPen(int num)
{
   DEPRECATED("HPGLPen()","blHPGLPen()");
   blHPGLPen(num);
}

void HPGLMove(REAL x, REAL y)
{
   DEPRECATED("HPGLMove()","blHPGLMove()");
   blHPGLMove(x, y);
}

void HPGLDraw(REAL x, REAL y)
{
   DEPRECATED("HPGLDraw()","blHPGLDraw()");
   blHPGLDraw(x, y);
}

void HPGLSetDash(int style)
{
   DEPRECATED("HPGLSetDash()","blHPGLSetDash()");
   blHPGLSetDash(style);
}

void HPGLFont(int font, REAL size)
{
   DEPRECATED("HPGLFont()","blHPGLFont()");
   blHPGLFont(font, size);
}

void HPGLLText(REAL x, REAL y, char *string)
{
   DEPRECATED("HPGLLText()","blHPGLLText()");
   blHPGLLText(x, y, string);
}

void HPGLCBText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLCBText()","blHPGLCBText()");
   blHPGLCBText(x, y, offset, text);
}

void HPGLROffText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLROffText()","blHPGLROffText()");
   blHPGLROffText(x, y, offset, text);
}

void HPGLLCText(REAL x, REAL y, char *text)
{
   DEPRECATED("HPGLLCText()","blHPGLLCText()");
   blHPGLLCText(x, y, text);
}

void HPGLCTText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLCTText()","blHPGLCTText()");
   blHPGLCTText(x, y, offset, text);
}

void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, REAL TitleSize, char *label, int LabelFont, REAL LabelSize)
{
   DEPRECATED("HPGLVText()","blHPGLVText()");
   blHPGLVText(x, y, xoff, text, TitleFont, TitleSize, label, LabelFont, LabelSize);
}

void HPGLEnd(void)
{
   DEPRECATED("HPGLEnd()","blHPGLEnd()");
   blHPGLEnd();
}

void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase)
{
   DEPRECATED("HPGLShowText()","blHPGLShowText()");
   blHPGLShowText(text, orientation, XBase, YBase);
}



/************************************************************************/
/* Renamed functions: matrix.h                                        */

void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout)
{
   DEPRECATED("MatMult3_33()","blMatMult3_33()");
   blMatMult3_33(vecin, matin, vecout);
}

void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
{
   DEPRECATED("MatMult33_33()","blMatMult33_33()");
   blMatMult33_33(a, b, out);
}

void invert33(REAL s[3][3], REAL ss[3][3])
{
   DEPRECATED("invert33()","blInvert33()");
   blInvert33(s, ss);
}

void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
{
   DEPRECATED("CreateRotMat()","blCreateRotMat()");
   blCreateRotMat(direction, angle, matrix);
}

REAL VecDist(REAL *a, REAL *b, int len)
{
   DEPRECATED("VecDist()","blVecDist()");
   return(blVecDist(a, b, len));
}



/************************************************************************/
/* Renamed functions: parse.h                                        */

int parse(char *comline, int nkeys, KeyWd *keywords, REAL *REALparam,
          char **strparam)
{
   DEPRECATED("parse()","blParse()");
   return(blParse(comline, nkeys, keywords, REALparam, strparam));
}

int mparse(char *comline, int nkeys, MKeyWd *keywords, REAL *REALparam,
           char **strparam, int *nparams)
{
   DEPRECATED("mparse()","blMparse()");
   return(blMparse(comline, nkeys, keywords, REALparam, strparam, 
                   nparams));
}
           
int match(char *comstring, char *string2, int *nletters)
{
   DEPRECATED("match()","blMatch()");
   return(blMatch(comstring, string2, nletters));
}

int GetString(char *command, char *strparam)
{
   DEPRECATED("GetString()","blGetString()");
   return(blGetString(command, strparam));
}

int GetParam(char  *command, REAL *value, int *nletters)
{
   DEPRECATED("GetParam()","blGetParam()");
   return(blGetParam(command, value, nletters));
}



/************************************************************************/
/* Renamed functions: plotting.h                                        */

BOOL AMInitPlot(char *filename, char *title, int dest, REAL OutXSize, 
                REAL OutYSize, REAL OutXOff, REAL OutYOff,
                char *AltFont, REAL xmargin, REAL ymargin,
                REAL DataXMin, REAL DataYMin, REAL DataXMax,
                REAL DataYMax)
{
   DEPRECATED("AMInitPlot()","blAMInitPlot()");
   return(blAMInitPlot(filename, title, dest, OutXSize, OutYSize, OutXOff,
                       OutYOff, AltFont, xmargin, ymargin, DataXMin, 
                       DataYMin, DataXMax,DataYMax));
}

void AMSetPen(int dest, int pen)
{
   DEPRECATED("AMSetPen()","blAMSetPen()");
   blAMSetPen(dest, pen);
}

void AMMove(int dest, REAL x, REAL y)
{
   DEPRECATED("AMMove()","blAMMove()");
   blAMMove(dest, x, y);
}

void AMDraw(int dest, REAL x, REAL y)
{
   DEPRECATED("AMDraw()","blAMDraw()");
   blAMDraw(dest, x, y);
}

void AMSetLineStyle(int dest, int style)
{
   DEPRECATED("AMSetLineStyle()","blAMSetLineStyle()");
   blAMSetLineStyle(dest, style);
}

void AMEndLine(int dest)
{
   DEPRECATED("AMEndLine()","blAMEndLine()");
   blAMEndLine(dest);
}

void AMSetFont(int dest, char *PSFontName, REAL FontSize)
{
   DEPRECATED("AMSetFont()","blAMSetFont()");
   blAMSetFont(dest, PSFontName, FontSize);
}

void AMText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMText()","blAMText()");
   blAMText(dest, x, y, text);
}

void AMCBText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMCBText()","blAMCBText()");
   blAMCBText(dest, x, y, text);
}

void AMRText(int dest, REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("AMRText()","blAMRText()");
   blAMRText(dest, x, y, offset, text);
}

void AMLCText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMLCText()","blAMLCText()");
   blAMLCText(dest, x, y, text);
}

void AMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text)
{
   DEPRECATED("AMCTText()","blAMCTText()");
   blAMCTText(dest, x, y, CTOffset, text);
}

void AMEndPlot(int dest)
{
   DEPRECATED("AMEndPlot()","blAMEndPlot()");
   blAMEndPlot(dest);
}

int  PS2HPGLFont(char *font);
char *SimplifyText(char *string)
{
   DEPRECATED("SimplifyText()","blSimplifyText()");
   return(blSimplifyText(string));
}



/************************************************************************/
/* Renamed functions: ps.h                                        */

BOOL PSInit(char *FName, char *creator, char *AltFont)
{
   DEPRECATED("PSInit()","blPSInit()");
   return(blPSInit(FName, creator, AltFont));
}

void PSThick(REAL thickness)
{
   DEPRECATED("PSThick()","blPSThick()");
   blPSThick(thickness);
}

void PSMove(REAL X, REAL Y)
{
   DEPRECATED("PSMove()","blPSMove()");
   blPSMove(X, Y);
}

void PSDraw(REAL X, REAL Y)
{
   DEPRECATED("PSDraw()","blPSDraw()");
   blPSDraw(X, Y);
}

void PSSetDash(char *linepatt)
{
   DEPRECATED("PSSetDash()","blPSSetDash()");
   blPSSetDash(linepatt);
}

void PSClearDash(void)
{
   DEPRECATED("PSClearDash()","blPSClearDash()");
   blPSClearDash();
}

void PSStroke(void)
{
   DEPRECATED("PSStroke()","blPSStroke()");
   blPSStroke();
}

void PSFont(char *fontname, REAL size)
{
   DEPRECATED("PSFont()","blPSFont()");
   blPSFont(fontname, size);
}

void PSLText(REAL X, REAL Y, char *label)
{
   DEPRECATED("PSLText()","blPSLText()");
   blPSLText(X, Y, label);
}

void PSCBText(REAL X, REAL Y, REAL Offset, char *label)
{
   DEPRECATED("PSCBText()","blPSCBText()");
   blPSCBText(X, Y, Offset, label);
}

void PSROffText(REAL X, REAL Y, REAL offset, char *label)
{
   DEPRECATED("PSROffText()","blPSROffText()");
   blPSROffText(X, Y, offset, label);
}

void PSLCText(REAL X, REAL Y, char *label)
{
   DEPRECATED("PSLCText()","blPSLCText()");
   blPSLCText(X, Y, label);
}

void PSCTText(REAL X, REAL Y, REAL Offset, char *label)
{
   DEPRECATED("PSCTText()","blPSCTText()");
   blPSCTText(X, Y, Offset, label);
}

void PSVText(REAL x, REAL y, REAL xoff, char *text, char *font, REAL size,
             char *label, char *lfont, REAL lsize)
{
   DEPRECATED("PSVText()","blPSVText()");
   blPSVText(x, y, xoff, text, font, size, label, lfont, lsize);
}

void PSShowText(char *text)
{
   DEPRECATED("PSShowText()","blPSShowText()");
   blPSShowText(text);
}

void PSEnd(void)
{
   DEPRECATED("PSEnd()","blPSEnd()");
   blPSEnd();
}

char *PSCorrectCase(char *font)
{
   DEPRECATED("PSCorrectCase()","blPSCorrectCase()");
   return(blPSCorrectCase(font));
}



/************************************************************************/
/* Renamed functions: safemem.h                                        */

void *safemalloc(int nbytes)
{
   DEPRECATED("Safemalloc()","blSafemalloc()");
   return(blSafemalloc(nbytes));
}

BOOL safefree(void *ptr)
{
   DEPRECATED("Safefree()","blSafefree()");
   return(blSafefree(ptr));
}

void safeleaks(void)
{
   DEPRECATED("Safeleaks()","blSafeleaks()");
   blSafeleaks();
}



/************************************************************************/
/* Renamed functions: seq.h                                             */

char throne(char *three)
{
   DEPRECATED("throne()","blThrone()");
   return(blThrone(three));
}

char thronex(char *three)
{
   DEPRECATED("thronex()","blThronex()");
   return(blThronex(three));
}

char *onethr(char one)
{
   DEPRECATED("onethr()","blOnethr()");
   return(blOnethr(one));
}

char *DoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX)
{
   DEPRECATED("DoPDB2Seq()","blDoPDB2Seq()");
   return(blDoPDB2Seq(pdb, DoAsxGlx, ProtOnly, NoX));
}

int SplitSeq(char *LinearSeq, char **seqs)
{
   DEPRECATED("SplitSeq()","blSplitSeq()");
   return(blSplitSeq(LinearSeq, seqs));
}

int ReadSimplePIR(FILE *fp, int  maxres, char **seqs)
{
   DEPRECATED("ReadSimplePIR()","blReadSimplePIR()");
   return(blReadSimplePIR(fp, maxres, seqs));
}


int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
            SEQINFO *seqinfo, BOOL *punct, BOOL *error)
{
   DEPRECATED("ReadPIR()","blReadPIR()");
   return(blReadPIR(fp, DoInsert, seqs, maxchain, seqinfo, punct, error));
}

int ReadRawPIR(FILE *fp, char **seqs, int maxchain, BOOL upcase,
               SEQINFO *seqinfo, BOOL *error)
{
   DEPRECATED("ReadRawPIR()","blReadRawPIR()");
   return(blReadRawPIR(fp, seqs, maxchain, upcase, seqinfo, error));
}

int align(char *seq1, int  length1, char *seq2, int  length2, 
          BOOL verbose, BOOL identity, int  penalty, 
          char *align1, char *align2, int  *align_len)
{
   DEPRECATED("align()","blAlign()");
   return(blAlign(seq1, length1, seq2, length2, 
                  verbose, identity, penalty, 
                  align1, align2, align_len));
}

int affinealign(char *seq1, int  length1, char *seq2, int  length2, 
                BOOL verbose, BOOL identity, int  penalty, int penext,
                char *align1, char *align2, int  *align_len)
{
   DEPRECATED("affinealign()","blAffinealign()");
   return(blAffinealign(seq1, length1, seq2, length2, 
                        verbose, identity, penalty, penext,
                        align1, align2, align_len));
}


int CalcMDMScore(char resa, char resb)
{
   DEPRECATED("CalcMDMScore()","blCalcMDMScore()");
   return(blCalcMDMScore(resa, resb));
}

int affinealignuc(char *seq1, int  length1, char *seq2, int  length2, 
                  BOOL verbose, BOOL identity, int  penalty, int penext,
                  char *align1, char *align2, int  *align_len)
{
   DEPRECATED("affinealignuc()","blAffinealignuc()");
   return(blAffinealignuc(seq1, length1, seq2, length2, 
                          verbose, identity, penalty, penext,
                          align1, align2, align_len));
}

int CalcMDMScoreUC(char resa, char resb)
{
   DEPRECATED("CalcMDMScoreUC()","blCalcMDMScoreUC()");
   return(blCalcMDMScoreUC(resa, resb));
}

BOOL ReadMDM(char *mdmfile)
{
   DEPRECATED("ReadMDM()","blReadMDM()");
   return(blReadMDM(mdmfile));
}

int ZeroMDM(void)
{
   DEPRECATED("ZeroMDM()","blZeroMDM()");
   return(blZeroMDM());
}

char DNAtoAA(char *dna)
{
   DEPRECATED("DNAtoAA()","blDNAtoAA()");
   return(blDNAtoAA(dna));
}

int TrueSeqLen(char *sequence)
{
   DEPRECATED("TrueSeqLen()","blTrueSeqLen()");
   return(blTrueSeqLen(sequence));
}

int KnownSeqLen(char *sequence)
{
   DEPRECATED("KnownSeqLen()","blKnownSeqLen()");
   return(blKnownSeqLen(sequence));
}

BOOL NumericReadMDM(char *mdmfile)
{
   DEPRECATED("NumericReadMDM()","blNumericReadMDM()");
   return(blNumericReadMDM(mdmfile));
}

int NumericCalcMDMScore(int resa, int resb)
{
   DEPRECATED("NumericCalcMDMScore()","blNumericCalcMDMScore()");
   return(blNumericCalcMDMScore(resa, resb));
}

int NumericAffineAlign(int *seq1, int length1, int *seq2, int length2, 
                       BOOL verbose, BOOL identity, int penalty,
                       int penext, int *align1, int *align2, 
                       int *align_len)
{
   DEPRECATED("NumericAffineAlign()","blNumericAffineAlign()");
   return(blNumericAffineAlign(seq1, length1, seq2, length2, 
                               verbose, identity, penalty,
                               penext, align1, align2, 
                               align_len));
}


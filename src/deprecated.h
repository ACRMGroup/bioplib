/************************************************************************/
/**

   \file       deprecated.h
   
   \version    v1.0
   \date       07.05.14
   \brief      Redirect calls to deprecated functions.
   
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


   Deprecated functions for BiopLib.

   Contains the DEPRECATED() macro which prints a warning message when a 
   call is made to a deprecated function.

   Contains macros with statement expressions that act as wrappers for 
   replacement functions. Calls to deprecated functions are redirected 
   their corresponding replacement functions.



**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   
-  V1.0  07.05.14 Original By: CTP
-  V1.1  07.07.14 Rename functions with 'bl' prefix. By: CTP

*************************************************************************/
#ifndef _DEPRECATED_H
#define _DEPRECATED_H

/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/************************************************************************/
/* Defines and macros
*/


/************************************************************************/
/*>DEPRECATED(s, t)
   ----------------
*//**

   The DEPRECATED macro gives a warning message if a function is 
   deprecated and indicates the replacement function.
   
   The default option is to give a warning message unless an environment 
   variable, BIOPLIB_DEPRECATED_QUIET, is set.
   
   Alternatively the compile options: -D BIOPLIB_DEPRECATED_CHECK or 
   -D BIOPLIB_DEPRECATED_QUIET will set the DEPRECATED macro to ignore the
   BIOPLIB_DEPRECATED_QUIET environment variable. 
   
   -D BIOPLIB_DEPRECATED_CHECK will display the warning message. 
   -D BIOPLIB_DEPRECATED_QUIET will always silence the warning message.
   
-  29.04.14 Original    By: CTP
*/
#if(defined BIOPLIB_DEPRECATED_CHECK || defined BIOPLIB_DEPRECATED_QUIET)
#  ifndef BIOPLIB_DEPRECATED_QUIET
#     define DEPRECATED(s, t)                                            \
               fprintf(stderr,                                           \
                       "This code uses %s which is now deprecated!\n"    \
                       "   Use %s instead\n", (s), (t))
#  else
#      define DEPRECATED(s, t)
#  endif
#else
#  define DEPRECATED(s, t)                                               \
   {                                                                     \
      if(!getenv("BIOPLIB_DEPRECATED_QUIET"))                            \
         fprintf(stderr,                                                 \
                 "This code uses %s which is now deprecated!\n"          \
                 "   Use %s instead\n", (s), (t));                       \
   }
#endif


/************************************************************************/
/*>PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list, but
  requires the residue is a HETATM record.
  Uses char for chain and insert.

  26.10.11 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindHetatmResidue(). By: CTP
  25.02.14 Added error message. By: CTP
  07.05.14 Converted to macro. By: CTP
  07.07.14 Use bl prefix for functions By: CTP
*/
#define FindHetatmResidue(pdb, chain, resnum, insert)                    \
({                                                                       \
   char chain_a[2]  = " ",                                               \
        insert_a[2] = " ";                                               \
                                                                         \
   DEPRECATED("FindHetatmResidue()","blFindHetatmResidue()");            \
                                                                         \
   chain_a[0]  = chain;                                                  \
   insert_a[0] = insert;                                                 \
                                                                         \
   blFindHetatmResidue(pdb, chain_a, resnum, insert_a);                  \
})


/************************************************************************/
/*>PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list.
  Uses char for string and insert.

  06.02.96 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindResidue(). By: CTP
  25.02.14 Added error message. By: CTP
  07.05.14 Converted to macro. By: CTP
  07.07.14 Use bl prefix for functions By: CTP
*/

#define FindResidue(pdb, chain, resnum, insert)                          \
({                                                                       \
   char chain_a[2]  = " ",                                               \
        insert_a[2] = " ";                                               \
                                                                         \
   DEPRECATED("FindResidue()","blFindResidue()");                        \
                                                                         \
   chain_a[0]  = chain;                                                  \
   insert_a[0] = insert;                                                 \
                                                                         \
   blFindResidue(pdb, chain_a, resnum, insert_a);                        \
})


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
*/
#define FindZonePDB(pdb, start, startinsert, stop, stopinsert, chain,    \
                    mode, pdb_start, pdb_stop)                           \
({                                                                       \
   char startinsert_a[2] = " ",                                          \
        stopinsert_a[2]  = " ",                                          \
        chain_a[2]       = " ";                                          \
                                                                         \
   DEPRECATED("FindZonePDB()","blFindZonePDB()");                        \
                                                                         \
   strncpy(startinsert_a,&startinsert,1);                                \
   strncpy(stopinsert_a ,&stopinsert ,1);                                \
   strncpy(chain_a      ,&chain      ,1);                                \
                                                                         \
   blFindZonePDB(pdb,                                                    \
                 start,                                                  \
                 startinsert_a,                                          \
                 stop,                                                   \
                 stopinsert_a,                                           \
                 chain_a,                                                \
                 mode,                                                   \
                 pdb_start,                                              \
                 pdb_stop);                                              \
})


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
*/
#define InPDBZone(p, chain, resnum1, insert1, resnum2, insert2)          \
({                                                                       \
   char chain_a[2]   = " ",                                              \
        insert1_a[2] = " ",                                              \
        insert2_a[2] = " ";                                              \
                                                                         \
   DEPRECATED("InPDBZone()","blInPDBZone()");                            \
                                                                         \
   chain_a[0]   = chain;                                                 \
   insert1_a[0] = insert1;                                               \
   insert2_a[0] = insert2;                                               \
                                                                         \
   blInPDBZone(p,chain_a,resnum1,insert1_a,resnum2, insert2_a);          \
})

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
*/
#define WritePDB(fp, pdb)                                                \
({                                                                       \
   DEPRECATED("WritePDB","blWritePDB()");                                \
   blWriteAsPDB(fp, pdb);                                                \
})

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
*/
#define WriteWholePDB(fp, wpdb)                                          \
({                                                                       \
   DEPRECATED("WriteWholePDB","blWriteWholePDB()");                      \
                                                                         \
   blWriteWholePDBHeader(fp, wpdb);                                      \
   blWritePDB(fp, wpdb->pdb);                                            \
   blWriteWholePDBTrailer(fp, wpdb);                                     \
})


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
*/
#define WriteWholePDBHeader(fp, wpdb)                                    \
({                                                                       \
   DEPRECATED("WriteWholePDBHeader()","blWriteWholePDBHeader()");        \
   blWriteWholePDBHeader(fp, wpdb);                                      \
})


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
*/
#define WriteWholePDBTrailer(fp, wpdb)                                   \
({                                                                       \
   DEPRECATED("WriteWholePDBTrailer()","blWriteWholePDBTrailer()");      \
   blWriteWholePDBTrailer(fp, wpdb);                                     \
})

/************************************************************************/
/* Renamed Functions 
*/
/*
#define XXX(a,b,c)                                                       \
({                                                                       \
   DEPRECATED("XXX","blXXX");                                            \
   blXXX(a,b,c);                                                         \
})
*/

/************************************************************************/
/* Renamed functions: ReadPDB.c                                         */

#define ReadPDB(fp, natom)                                               \
({                                                                       \
   DEPRECATED("ReadPDB()","blReadPDB()");                                \
   blReadPDB(fp, natom);                                                 \
})

#define ReadPDBAll(fp, natom)                                            \
({                                                                       \
   DEPRECATED("ReadPDBAll()","blReadPDBAll()");                          \
   blReadPDBAll(fp, natom);                                              \
})

#define ReadPDBAtoms(fp, natom)                                          \
({                                                                       \
   DEPRECATED("ReadPDBAtoms()","blReadPDBAtoms()");                      \
   blReadPDBAtoms(fp, natom);                                            \
})

#define ReadPDBOccRank(fp, natom, OccRank)                               \
({                                                                       \
   DEPRECATED("ReadPDBOccRank()","blReadPDBOccRank()");                  \
   blReadPDBOccRank(fp, natom, OccRank);                                 \
})

#define ReadPDBAtomsOccRank(fp, natom, OccRank)                          \
({                                                                       \
   DEPRECATED("ReadPDBAtomsOccRank()","blReadPDBAtomsOccRank()");        \
   blReadPDBAtomsOccRank(fp, natom, OccRank);                            \
})

#define doReadPDB(fpin, natom, AllAtoms, OccRank, ModelNum)              \
({                                                                       \
   DEPRECATED("doReadPDB()","bldoReadPDB()");                            \
   bldoReadPDB(fpin, natom, AllAtoms, OccRank, ModelNum);                \
})

#define FixAtomName(name, occup)                                         \
({                                                                       \
   DEPRECATED("FixAtomName()","blFixAtomName()");                        \
   blFixAtomName(name, occup);                                           \
})

#define RemoveAlternates(pdb)                                            \
({                                                                       \
   DEPRECATED("RemoveAlternates","blRemoveAlternates");                  \
   blRemoveAlternates(pdb);                                              \
})

#define doReadPDBML(fp, natom, AllAtoms, OccRank, ModelNum)              \
({                                                                       \
   DEPRECATED("doReadPDBML()","bldoReadPDBML()");                        \
   bldoReadPDBML(fp, natom, AllAtoms, OccRank, ModelNum);                \
})

#define CheckFileFormatPDBML(fp)                                         \
({                                                                       \
   DEPRECATED("CheckFileFormatPDBML()","blCheckFileFormatPDBML()");      \
   blCheckFileFormatPDBML(fp);                                           \
})

/************************************************************************/
/* Renamed functions: WritePDB.c                                        */

#define WritePDBRecord(fp, pdb)                                          \
({                                                                       \
   DEPRECATED("WritePDBRecord()","blWritePDBRecord()");                  \
   blWritePDBRecord(fp, pdb);                                            \
})

#define WritePDBRecordAtnam(fp, pdb)                                     \
({                                                                       \
   DEPRECATED("WritePDBRecordAtnam()","blWritePDBRecordAtnam()");        \
   blWritePDBRecordAtnam(fp, pdb);                                       \
})

/************************************************************************/
/* Renamed functions: WholePDB.c                                        */

#define FreeWholePDB(wpdb)                                               \
({                                                                       \
   DEPRECATED("FreeWholePDB()","blFreeWholePDB()");                      \
   blFreeWholePDB(wpdb);                                                 \
})

#define ReadWholePDB(fpin)                                               \
({                                                                       \
   DEPRECATED("ReadWholePDB()","blReadWholePDB()");                      \
   blReadWholePDB(fpin);                                                 \
})

#define ReadWholePDBAtoms(fpin)                                          \
({                                                                       \
   DEPRECATED("ReadWholePDBAtoms()","blReadWholePDBAtoms()");            \
   blReadWholePDBAtoms(fpin);                                            \
})

/************************************************************************/
/* Renamed functions: pdb.h                                             */

#define WriteGromosPDB(fp, pdb)                                          \
({                                                                       \
   DEPRECATED("WriteGromosPDB()","blWriteGromosPDB()");                  \
   blWriteGromosPDB(fp, pdb);                                            \
})
#define WriteGromosPDBRecord(fp, pdb)                                    \
({                                                                       \
   DEPRECATED("WriteGromosPDBRecord()","blWriteGromosPDBRecord()");      \
   blWriteGromosPDBRecord(fp, pdb);                                      \
})

#define GetCofGPDB(pdb, cg)                                              \
({                                                                       \
   DEPRECATED("GetCofGPDB()","blGetCofGPDB()");                          \
   blGetCofGPDB(pdb, cg);                                                \
})

#define GetCofGPDBRange(start, stop, cg)                                 \
({                                                                       \
   DEPRECATED("GetCofGPDBRange()","blGetCofGPDBRange()");                \
   blGetCofGPDBRange(start, stop, cg);                                   \
})

#define GetCofGPDBSCRange(start, stop, cg)                               \
({                                                                       \
   DEPRECATED("GetCofGPDBSCRange()","blGetCofGPDBSCRange()");            \
   blGetCofGPDBSCRange(start, stop, cg);                                 \
})

#define OriginPDB(pdb)                                                   \
({                                                                       \
   DEPRECATED("OriginPDB()","blOriginPDB()");                            \
   blOriginPDB(pdb);                                                     \
})

#define RotatePDB(pdb, rm)                                               \
({                                                                       \
   DEPRECATED("RotatePDB()","blRotatePDB()");                            \
   blRotatePDB(pdb, rm);                                                 \
})

#define TranslatePDB(pdb, tvect)                                         \
({                                                                       \
   DEPRECATED("TranslatePDB()","blTranslatePDB()");                      \
   blTranslatePDB(pdb, tvect);                                           \
})


/*345678901234567890123456789012345678901234567890123456789012345678901234567890*/
#define FitPDB(ref_pdb, fit_pdb, rm)                                     \
({                                                                       \
   DEPRECATED("FitPDB()","blFitPDB()");                                  \
   blFitPDB(ref_pdb, fit_pdb, rm);                                       \
})

#define FitCaPDB(ref_pdb, fit_pdb, rm)                                   \
({                                                                       \
   DEPRECATED("FitCaPDB()","blFitCaPDB()");                              \
   blFitCaPDB(ref_pdb, fit_pdb, rm);                                     \
})

#define FitNCaCPDB(ref_pdb, fit_pdb, rm)                                 \
({                                                                       \
   DEPRECATED("FitNCaCPDB()","blFitNCaCPDB()");                          \
   blFitNCaCPDB(ref_pdb, fit_pdb, rm);                                   \
})

#define FitCaCbPDB(ref_pdb, fit_pdb, rm)                                 \
({                                                                       \
   DEPRECATED("FitCaCbPDB()","blFitCaCbPDB()");                          \
   blFitCaCbPDB(ref_pdb, fit_pdb, rm);                                   \
})

#define CalcRMSPDB(pdb1, pdb2)                                           \
({                                                                       \
   DEPRECATED("CalcRMSPDB()","blCalcRMSPDB()");                          \
   blCalcRMSPDB(pdb1, pdb2);                                             \
})

#define GetPDBCoor(pdb, coor)                                            \
({                                                                       \
   DEPRECATED("GetPDBCoor()","blGetPDBCoor()");                          \
   blGetPDBCoor(pdb, coor);                                              \
})

#define HAddPDB(fp, pdb)                                                 \
({                                                                       \
   DEPRECATED("HAddPDB()","blHAddPDB()");                                \
   blHAddPDB(fp, pdb);                                                   \
})

#define ReadPGP(fp)                                                      \
({                                                                       \
   DEPRECATED("ReadPGP()","blReadPGP()");                                \
   blReadPGP(fp);                                                        \
})

#define OpenPGPFile(pgpfile, AllHyd)                                     \
({                                                                       \
   DEPRECATED("OpenPGPFile()","blOpenPGPFile()");                        \
   blOpenPGPFile(pgpfile, AllHyd);                                       \
})

#define SelectAtomsPDB(pdbin, nsel, sel, natom)                          \
({                                                                       \
   DEPRECATED("SelectAtomsPDB()","blSelectAtomsPDB()");                  \
   blSelectAtomsPDB(pdbin, nsel, sel, natom);                            \
})

#define StripHPDB(pdbin, natom)                                          \
({                                                                       \
   DEPRECATED("StripHPDB()","blStripHPDB()");                            \
   blStripHPDB(pdbin, natom);                                            \
})

#define ReadSecPDB(fp, nsec)                                             \
({                                                                       \
   DEPRECATED("ReadSecPDB()","blReadSecPDB()");                          \
   blReadSecPDB(fp, nsec);                                               \
})

#define RenumAtomsPDB(pdb)                                               \
({                                                                       \
   DEPRECATED("RenumAtomsPDB()","blRenumAtomsPDB()");                    \
   blRenumAtomsPDB(pdb);                                                 \
})

#define FindEndPDB(start)                                                \
({                                                                       \
   DEPRECATED("FindEndPDB()","blFindEndPDB()");                          \
   blFindEndPDB(start);                                                  \
})

#define FixOrderPDB(pdb, Pad, Renum)                                     \
({                                                                       \
   DEPRECATED("FixOrderPDB()","blFixOrderPDB()");                        \
   blFixOrderPDB(pdb, Pad, Renum);                                       \
})

#define ShuffleResPDB(start, end, Pad)                                   \
({                                                                       \
   DEPRECATED("ShuffleResPDB()","blShuffleResPDB()");                    \
   blShuffleResPDB(start, end, Pad);                                     \
})

#define GetAtomTypes(resnam, AtomTypes)                                  \
({                                                                       \
   DEPRECATED("GetAtomTypes()","blGetAtomTypes()");                      \
   blGetAtomTypes(resnam, AtomTypes);                                    \
})

#define KillPDB(pdb, prev)                                               \
({                                                                       \
   DEPRECATED("KillPDB()","blKillPDB()");                                \
   blKillPDB(pdb, prev);                                                 \
})

#define CopyPDB(out, in)                                                 \
({                                                                       \
   DEPRECATED("CopyPDB()","blCopyPDB()");                                \
   blCopyPDB(out, in);                                                   \
})

#define MovePDB(move, from, to)                                          \
({                                                                       \
   DEPRECATED("MovePDB()","blMovePDB()");                                \
   blMovePDB(move, from, to);                                            \
})

#define AppendPDB(first, second)                                         \
({                                                                       \
   DEPRECATED("AppendPDB()","blAppendPDB()");                            \
   blAppendPDB(first, second);                                           \
})

#define ShuffleBB(pdb)                                                   \
({                                                                       \
   DEPRECATED("ShuffleBB()","blShuffleBB()");                            \
   blShuffleBB(pdb);                                                     \
})

#define CalcChi(pdb, type)                                               \
({                                                                       \
   DEPRECATED("CalcChi()","blCalcChi()");                                \
   blCalcChi(pdb, type);                                                 \
})

#define GetPDBByN(pdb, n)                                                \
({                                                                       \
   DEPRECATED("GetPDBByN()","blGetPDBByN()");                            \
   blGetPDBByN(pdb, n);                                                  \
})

#define SetChi(pdb, next, chi, type)                                     \
({                                                                       \
   DEPRECATED("SetChi()","blSetChi()");                                  \
   blSetChi(pdb, next, chi, type);                                       \
})

#define KillSidechain(ResStart, NextRes, doCB)                           \
({                                                                       \
   DEPRECATED("KillSidechain()","blKillSidechain()");                    \
   blKillSidechain(ResStart, NextRes, doCB);                             \
})

#define SetResnam(ResStart, NextRes, resnam, resnum, insert, chain)      \
({                                                                       \
   DEPRECATED("SetResnam()","blSetResnam()");                            \
   blSetResnam(ResStart, NextRes, resnam, resnum, insert, chain);        \
})

#define ApplyMatrixPDB(pdb, rm)                                          \
({                                                                       \
   DEPRECATED("ApplyMatrixPDB()","blApplyMatrixPDB()");                  \
   blApplyMatrixPDB(pdb, rm);                                            \
})


#define GetResolPDB(fp, resolution, RFactor, StrucType)                  \
({                                                                       \
   DEPRECATED("GetResolPDB()","blGetResolPDB()");                        \
   blGetResolPDB(fp, resolution, RFactor, StrucType);                    \
})

#define GetExptl(fp, resolution, RFactor, FreeR, StrucType)              \
({                                                                       \
   DEPRECATED("GetExptl()","blGetExptl()");                              \
   blGetExptl(fp, resolution, RFactor, FreeR, StrucType);                \
})

#define GetExptlOld(fp, resolution, RFactor, FreeR, StrucType)           \
({                                                                       \
   DEPRECATED("GetExptlOld()","blGetExptlOld()");                        \
   blGetExptlOld(fp, resolution, RFactor, FreeR, StrucType);             \
})

#define ReportStructureType(type)                                        \
({                                                                       \
   DEPRECATED("ReportStructureType()","blReportStructureType()");        \
   blReportStructureType(type);                                          \
})

#define IndexPDB(pdb, natom)                                             \
({                                                                       \
   DEPRECATED("IndexPDB()","blIndexPDB()");                              \
   blIndexPDB(pdb, natom);                                               \
})

#define ReadDisulphidesPDB(fp, error)                                    \
({                                                                       \
   DEPRECATED("ReadDisulphidesPDB()","blReadDisulphidesPDB()");          \
   blReadDisulphidesPDB(fp, error);                                      \
})

#define ParseResSpec(spec, chain, resnum, insert)                        \
({                                                                       \
   DEPRECATED("ParseResSpec()","blParseResSpec()");                      \
   blParseResSpec(spec, chain, resnum, insert);                          \
})

#define ParseResSpecNoUpper(spec, chain, resnum, insert)                 \
({                                                                       \
   DEPRECATED("ParseResSpecNoUpper()","blParseResSpecNoUpper()");        \
   blParseResSpecNoUpper(spec, chain, resnum, insert);                   \
})

#define DoParseResSpec(spec, chain, resnum, insert, uppercaseresspec)    \
({                                                                       \
   DEPRECATED("DoParseResSpec()","blDoParseResSpec()");                  \
   blDoParseResSpec(spec, chain, resnum, insert, uppercaseresspec);      \
})

#define RepSChain(pdb, sequence, ChiTable, RefCoords)                    \
({                                                                       \
   DEPRECATED("RepSChain()","blRepSChain()");                            \
   blRepSChain(pdb, sequence, ChiTable, RefCoords);                      \
})

#define FindNextChainPDB(pdb)                                            \
({                                                                       \
   DEPRECATED("FindNextChainPDB()","blFindNextChainPDB()");              \
   blFindNextChainPDB(pdb);                                              \
})

#define FixCterPDB(pdb, style)                                           \
({                                                                       \
   DEPRECATED("FixCterPDB()","blFixCterPDB()");                          \
   blFixCterPDB(pdb, style);                                             \
})

#define CalcCterCoords(p, ca_p, c_p, o_p)                                \
({                                                                       \
   DEPRECATED("CalcCterCoords()","blCalcCterCoords()");                  \
   blCalcCterCoords(p, ca_p, c_p, o_p);                                  \
})

#define CalcTetraHCoords(nter, coor)                                     \
({                                                                       \
   DEPRECATED("CalcTetraHCoords()","blCalcTetraHCoords()");              \
   blCalcTetraHCoords(nter, coor);                                       \
})

#define AddNTerHs(ppdb, Charmm)                                          \
({                                                                       \
   DEPRECATED("AddNTerHs()","blAddNTerHs()");                            \
   blAddNTerHs(ppdb, Charmm);                                            \
})

#define FNam2PDB(filename)                                               \
({                                                                       \
   DEPRECATED("FNam2PDB()","blFNam2PDB()");                              \
   blFNam2PDB(filename);                                                 \
})

#define TermPDB(pdb, length)                                             \
({                                                                       \
   DEPRECATED("TermPDB()","blTermPDB()");                                \
   blTermPDB(pdb, length);                                               \
})

#define FindHetatmResidueSpec(pdb, resspec)                              \
({                                                                       \
   DEPRECATED("FindHetatmResidueSpec()","blFindHetatmResidueSpec()");    \
   blFindHetatmResidueSpec(pdb, resspec);                                \
})

#define FindResidueSpec(pdb, resspec)                                    \
({                                                                       \
   DEPRECATED("FindResidueSpec()","blFindResidueSpec()");                \
   blFindResidueSpec(pdb, resspec);                                      \
})

#define FindNextResidue(pdb)                                             \
({                                                                       \
   DEPRECATED("FindNextResidue()","blFindNextResidue()");                \
   blFindNextResidue(pdb);                                               \
})

#define DupePDB(in)                                                      \
({                                                                       \
   DEPRECATED("DupePDB()","blDupePDB()");                                \
   blDupePDB(in);                                                        \
})

#define CopyPDBCoords(out, in)                                           \
({                                                                       \
   DEPRECATED("CopyPDBCoords()","blCopyPDBCoords()");                    \
   blCopyPDBCoords(out, in);                                             \
})

#define CalcCellTrans(UnitCell, CellAngles, xtrans, ytrans, ztrans)      \
({                                                                       \
   DEPRECATED("CalcCellTrans()","blCalcCellTrans()");                    \
   blCalcCellTrans(UnitCell, CellAngles, xtrans, ytrans, ztrans);        \
})

#define GetCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix) \
({                                                                       \
   DEPRECATED("GetCrystPDB()","blGetCrystPDB()");                        \
   blGetCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix); \
})

#define WriteCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix) \
({                                                                       \
   DEPRECATED("WriteCrystPDB()","blWriteCrystPDB()");                    \
   blWriteCrystPDB(fp, UnitCell, CellAngles, spacegroup, OrigMatrix, ScaleMatrix); \
})

#define ExtractZonePDB(inpdb, chain1, resnum1, insert1, chain2, resnum2, insert2) \
({                                                                       \
   DEPRECATED("ExtractZonePDB()","blExtractZonePDB()");                  \
   blExtractZonePDB(inpdb, chain1, resnum1, insert1, chain2, resnum2, insert2); \
})

#define FindAtomInRes(pdb, atnam)                                        \
({                                                                       \
   DEPRECATED("FindAtomInRes()","blFindAtomInRes()");                    \
   blFindAtomInRes(pdb, atnam);                                          \
})

#define InPDBZoneSpec(p, resspec1, resspec2)                             \
({                                                                       \
   DEPRECATED("InPDBZoneSpec()","blInPDBZoneSpec()");                    \
   blInPDBZoneSpec(p, resspec1, resspec2);                               \
})

#define AtomNameMatch(atnam, spec, ErrorWarn)                            \
({                                                                       \
   DEPRECATED("AtomNameMatch()","blAtomNameMatch()");                    \
   blAtomNameMatch(atnam, spec, ErrorWarn);                              \
})

#define AtomNameRawMatch(atnam, spec, ErrorWarn)                         \
({                                                                       \
   DEPRECATED("AtomNameRawMatch()","blAtomNameRawMatch()");              \
   blAtomNameRawMatch(atnam, spec, ErrorWarn);                           \
})

#define LegalAtomSpec(spec)                                              \
({                                                                       \
   DEPRECATED("LegalAtomSpec()","blLegalAtomSpec()");                    \
   blLegalAtomSpec(spec);                                                \
})

#define RepOneSChain(pdb, ResSpec, aa, ChiTable, RefCoords)              \
({                                                                       \
   DEPRECATED("RepOneSChain()","blRepOneSChain()");                      \
   blRepOneSChain(pdb, ResSpec, aa, ChiTable, RefCoords);                \
})

/* Check this OK */
#define EndRepSChain()                                                   \
({                                                                       \
   DEPRECATED("EndRepSChain()","blEndRepSChain()");                      \
   blEndRepSChain();                                                     \
})

#define ReadSeqresPDB(fp, nchains)                                       \
({                                                                       \
   DEPRECATED("ReadSeqresPDB()","blReadSeqresPDB()");                    \
   blReadSeqresPDB(fp, nchains);                                         \
})

#define SelectCaPDB(pdb)                                                 \
({                                                                       \
   DEPRECATED("SelectCaPDB()","blSelectCaPDB()");                        \
   blSelectCaPDB(pdb);                                                   \
})

#define AddCBtoGly(pdb)                                                  \
({                                                                       \
   DEPRECATED("AddCBtoGly()","blAddCBtoGly()");                          \
   blAddCBtoGly(pdb);                                                    \
})

#define AddCBtoAllGly(pdb)                                               \
({                                                                       \
   DEPRECATED("AddCBtoAllGly()","blAddCBtoAllGly()");                    \
   blAddCBtoAllGly(pdb);                                                 \
})

#define StripGlyCB(pdb)                                                  \
({                                                                       \
   DEPRECATED("StripGlyCB()","blStripGlyCB()");                          \
   blStripGlyCB(pdb);                                                    \
})

#define BuildAtomNeighbourPDBList(pdb, pRes, NeighbDist)                 \
({                                                                       \
   DEPRECATED("BuildAtomNeighbourPDBList()","blBuildAtomNeighbourPDBList()"); \
   blBuildAtomNeighbourPDBList(pdb, pRes, NeighbDist);                   \
})

#define FindAtomWildcardInRes(pdb, pattern)                              \
({                                                                       \
   DEPRECATED("FindAtomWildcardInRes()","blFindAtomWildcardInRes()");    \
   blFindAtomWildcardInRes(pdb, pattern);                                \
})

#define DupeResiduePDB(in)                                               \
({                                                                       \
   DEPRECATED("DupeResiduePDB()","blDupeResiduePDB()");                  \
   blDupeResiduePDB(in);                                                 \
})

#define StripWatersPDB(pdbin, natom)                                     \
({                                                                       \
   DEPRECATED("StripWatersPDB()","blStripWatersPDB()");                  \
   blStripWatersPDB(pdbin, natom);                                       \
})

#define AllocPDBStructure(pdb)                                           \
({                                                                       \
   DEPRECATED("AllocPDBStructure()","blAllocPDBStructure()");            \
   blAllocPDBStructure(pdb);                                             \
})

#define FindNextChain(pdb)                                               \
({                                                                       \
   DEPRECATED("FindNextChain()","blFindNextChain()");                    \
   blFindNextChain(pdb);                                                 \
})

#define FreePDBStructure(pdbstruct)                                      \
({                                                                       \
   DEPRECATED("FreePDBStructure()","blFreePDBStructure()");              \
   blFreePDBStructure(pdbstruct);                                        \
})


/************************************************************************/
/* Renamed functions: hbond.h                                           */

#define IsHBonded(res1, res2, type)                                      \
({                                                                       \
   DEPRECATED("IsHBonded()","blIsHBonded()");                            \
   blIsHBonded(res1, res2, type);                                        \
})

#define ValidHBond(AtomH, AtomD, AtomA, AtomP)                           \
({                                                                       \
   DEPRECATED("ValidHBond()","blValidHBond()");                          \
   blValidHBond(AtomH, AtomD, AtomA, AtomP);                             \
})

#define IsMCDonorHBonded(res1, res2, type)                               \
({                                                                       \
   DEPRECATED("IsMCDonorHBonded()","blIsMCDonorHBonded()");              \
   blIsMCDonorHBonded(res1, res2, type);                                 \
})

#define IsMCAcceptorHBonded(res1, res2, type)                            \
({                                                                       \
   DEPRECATED("IsMCAcceptorHBonded()","blIsMCAcceptorHBonded()");        \
   blIsMCAcceptorHBonded(res1, res2, type);                              \
})

/************************************************************************/
/* Renamed functions: BuffInp.h                                         */


#define OpenBufferedFile(filename, maxstr)                               \
({                                                                       \
   DEPRECATED("OpenBufferedFile()","blOpenBufferedFile()");              \
   blOpenBufferedFile(filename, maxstr);                                 \
})

#define ReadBufferedFile(bfp, string, length)                            \
({                                                                       \
   DEPRECATED("ReadBufferedFile()","blReadBufferedFile()");              \
   blReadBufferedFile(bfp, string, length);                              \
})

#define ProbeBufferedFile(bfp, string, length)                           \
({                                                                       \
   DEPRECATED("ProbeBufferedFile()","blProbeBufferedFile()");            \
   blProbeBufferedFile(bfp, string, length);                             \
})



/************************************************************************/
/* Renamed functions: ErrStack.h                                        */

#define StoreError(routine, error)                                       \
({                                                                       \
   DEPRECATED("StoreError()","blStoreError()");                          \
   blStoreError(routine, error);                                         \
})

#define ShowErrors(PrintRoutine, Trace)                                  \
({                                                                       \
   DEPRECATED("ShowErrors()","blShowErrors()");                          \
   blShowErrors(PrintRoutine, Trace);                                    \
})



/************************************************************************/
/* Renamed functions: MathUtil.h                                        */

#define CalcSD(val, action, mean, SD)                                    \
({                                                                       \
   DEPRECATED("CalcSD()","blCalcSD()");                                  \
   blCalcSD(val, action, mean, SD);                                      \
})
#define CalcExtSD(val, action, Sx, SxSq, NValues, mean, SD)              \
({                                                                       \
   DEPRECATED("CalcExtSD()","blCalcExtSD()");                            \
   blCalcExtSD(val, action, Sx, SxSq, NValues, mean, SD);                \
})
#define pearson(x, y, NItem)                                             \
({                                                                       \
   DEPRECATED("pearson()","blpearson()");                                \
   blpearson(x, y, NItem);                                               \
})
#define pearson1(x, y, NItem)                                            \
({                                                                       \
   DEPRECATED("pearson1()","blpearson1()");                              \
   blpearson1(x, y, NItem);                                              \
})

#define CrossProd3(Out, In1, In2)                                        \
({                                                                       \
   DEPRECATED("CrossProd3()","blCrossProd3()");                          \
   blCrossProd3(Out, In1, In2);                                          \
})

#define VecSub3(Out, In1, In2)                                           \
({                                                                       \
   DEPRECATED("VecSub3()","blVecSub3()");                                \
   blVecSub3(Out, In1, In2);                                             \
})

#define VecAdd3(Out, In1, In2)                                           \
({                                                                       \
   DEPRECATED("VecAdd3()","blVecAdd3()");                                \
   blVecAdd3(Out, In1, In2);                                             \
})

#define VecLen3(Vec)                                                     \
({                                                                       \
   DEPRECATED("VecLen3()","blVecLen3()");                                \
   blVecLen3(Vec);                                                       \
})

#define DistPtVect(Point, End1, End2)                                    \
({                                                                       \
   DEPRECATED("DistPtVect()","blDistPtVect()");                          \
   blDistPtVect(Point, End1, End2);                                      \
})

#define PointLineDistance(Px, Py, Pz, P1x, P1y, P1z, P2x, P2y, P2z, Rx, Ry, Rz, frac) \
({                                                                       \
   DEPRECATED("PointLineDistance()","blPointLineDistance()");            \
   blPointLineDistance(Px, Py, Pz, P1x, P1y, P1z, P2x, P2y, P2z, Rx, Ry, Rz, frac)  \
})

#define factorial(n)                                                     \
({                                                                       \
   DEPRECATED("factorial()","blfactorial()");                            \
   blfactorial(n);                                                       \
})

#define factdiv(n1, n2)                                                  \
({                                                                       \
   DEPRECATED("factdiv()","blfactdiv()");                                \
   blfactdiv(n1, n2);                                                    \
})

#define NPerm(n, r)                                                      \
({                                                                       \
   DEPRECATED("NPerm()","blNPerm()");                                    \
   blNPerm(n, r);                                                        \
})

#define NComb(n, r)                                                      \
({                                                                       \
   DEPRECATED("NComb()","blNComb()");                                    \
   blNComb(n, r);                                                        \
})

/************************************************************************/
/* Renamed functions: WindIO.h                                          */

#define screen(string)                                                   \
({                                                                       \
   DEPRECATED("screen()","blscreen()");                                  \
   blscreen(string);                                                     \
})

#define prompt(string)                                                   \
({                                                                       \
   DEPRECATED("prompt()","blprompt()");                                  \
   blprompt(string);                                                     \
})

#define RePrompt()                                                       \
({                                                                       \
   DEPRECATED("RePrompt()","blRePrompt()");                              \
   blRePrompt();                                                         \
})

#define GetKybdString(string, maxlen)                                    \
({                                                                       \
   DEPRECATED("GetKybdString()","blGetKybdString()");                    \
   blGetKybdString(string, maxlen);                                      \
})

#define PagingOn()                                                       \
({                                                                       \
   DEPRECATED("PagingOn()","blPagingOn()");                              \
   blPagingOn();                                                         \
})

#define PagingOff()                                                      \
({                                                                       \
   DEPRECATED("PagingOff()","blPagingOff()");                            \
   blPagingOff();                                                        \
})

#define WindowMode(mode)                                                 \
({                                                                       \
   DEPRECATED("WindowMode()","blWindowMode()");                          \
   blWindowMode(mode);                                                   \
})

#define WindowInteractive(mode)                                          \
({                                                                       \
   DEPRECATED("WindowInteractive()","blWindowInteractive()");            \
   blWindowInteractive(mode);                                            \
})

#define YorN(deflt)                                                      \
({                                                                       \
   DEPRECATED("YorN()","blYorN()");                                      \
   blYorN(deflt);                                                        \
})






/************************************************************************/
/* Renamed functions: aalist.h                                          */

#define InsertNextResiduesInAAList(a, res, nres)                         \
({                                                                       \
   DEPRECATED("InsertNextResiduesInAAList()","blInsertNextResiduesInAAList()");\
   blInsertNextResiduesInAAList(a, res, nres);                           \
})

#define InsertNextResidueInAAList(a, res)                                \
({                                                                       \
   DEPRECATED("InsertNextResidueInAAList()","blInsertNextResidueInAAList()");\
   blInsertNextResidueInAAList(a, res);                                  \
})

#define BuildSeqFromAAList(aa)                                           \
({                                                                       \
   DEPRECATED("BuildSeqFromAAList()","blBuildSeqFromAAList()");          \
   blBuildSeqFromAAList(aa);                                             \
})

#define InsertResidueInAAListAt(aa, res, pos)                            \
({                                                                       \
   DEPRECATED("InsertResidueInAAListAt()","blInsertResidueInAAListAt()");\
   blInsertResidueInAAListAt(aa, res, pos);                              \
})

#define InsertResiduesInAAListAt(aa, res, nres, pos)                     \
({                                                                       \
   DEPRECATED("InsertResiduesInAAListAt()","blInsertResiduesInAAListAt()");\
   blInsertResiduesInAAListAt(aa, res, nres, pos);                       \
})

#define BuildAAList(seq)                                                 \
({                                                                       \
   DEPRECATED("BuildAAList()","blBuildAAList()");                        \
   blBuildAAList(seq);                                                   \
})

#define FindAAListOffsetByResnum(aa, resnum)                             \
({                                                                       \
   DEPRECATED("FindAAListOffsetByResnum()","blFindAAListOffsetByResnum()");\
   blFindAAListOffsetByResnum(aa, resnum);                               \
})

#define FindAAListItemByResnum(aa, resnum)                               \
({                                                                       \
   DEPRECATED("FindAAListItemByResnum()","blFindAAListItemByResnum()");  \
   blFindAAListItemByResnum(aa, resnum);                                 \
})

#define SetAAListFlagByResnum(aa, resnum)                                \
({                                                                       \
   DEPRECATED("SetAAListFlagByResnum()","blSetAAListFlagByResnum()");    \
   blSetAAListFlagByResnum(aa, resnum);                                  \
})

#define BuildFlagSeqFromAAList(aa, ch)                                   \
({                                                                       \
   DEPRECATED("BuildFlagSeqFromAAList()","blBuildFlagSeqFromAAList()");  \
   blBuildFlagSeqFromAAList(aa, ch);                                     \
})

#define GetAAListLen(aa)                                                 \
({                                                                       \
   DEPRECATED("GetAAListLen()","blGetAAListLen()");                      \
   blGetAAListLen(aa);                                                   \
})



/************************************************************************/
/* Renamed functions: angle.h                                           */

#define angle(xi, yi, zi, xj, yj, zj, xk, yk, zk)                        \
({                                                                       \
   DEPRECATED("angle()","blangle()");                                    \
   blangle(xi, yi, zi, xj, yj, zj, xk, yk, zk);                          \
})

#define phi(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl)              \
({                                                                       \
   DEPRECATED("phi()","blphi()");                                        \
   blphi(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl);                \
})

#define simpleangle(ang)                                                 \
({                                                                       \
   DEPRECATED("simpleangle()","blsimpleangle()");                        \
   blsimpleangle(ang);                                                   \
})

#define TrueAngle(opp, adj)                                              \
({                                                                       \
   DEPRECATED("TrueAngle()","blTrueAngle()");                            \
   blTrueAngle(opp, adj);                                                \
})

#define TorToCoor(ant1, ant2, ant3, bond, theta, torsion, coords)        \
({                                                                       \
   DEPRECATED("TorToCoor()","blTorToCoor()");                            \
   blTorToCoor(ant1, ant2, ant3, bond, theta, torsion, coords);          \
})



/************************************************************************/
/* Renamed functions: array.h                                           */

#define Array2D(size, dim1, dim2)                                        \
({                                                                       \
   DEPRECATED("Array2D()","blArray2D()");                                \
   blArray2D(size, dim1, dim2);                                          \
})

#define FreeArray2D(array, dim1, dim2)                                   \
({                                                                       \
   DEPRECATED("FreeArray2D()","blFreeArray2D()");                        \
   blFreeArray2D(array, dim1, dim2);                                     \
})


#define Array3D(size, dim1, dim2, dim3)                                  \
({                                                                       \
   DEPRECATED("Array3D()","blArray3D()");                                \
   blArray3D(size, dim1, dim2, dim3);                                    \
})

#define FreeArray3D(array, dim1, dim2, dim3)                             \
({                                                                       \
   DEPRECATED("FreeArray3D()","blFreeArray3D()");                        \
   blFreeArray3D(array, dim1, dim2, dim3);                               \
})



/************************************************************************/
/* Renamed functions: cssr.h                                            */
/*
#define ReadCSSR(fp, natom, name, title)                                 \
({                                                                       \
   DEPRECATED("ReadCSSR()","blReadCSSR()");                              \
   blReadCSSR(fp, natom, name, title);                                   \
})

#define ReadCSSRasPDB(fp, natom)                                         \
({                                                                       \
   DEPRECATED("ReadCSSRasPDB()","blReadCSSRasPDB()");                    \
   blReadCSSRasPDB(fp, natom);                                           \
})

#define NormaliseCSSR(cssr, cell, alpha, beta, gamma)                    \
({                                                                       \
   DEPRECATED("NormaliseCSSR()","blNormaliseCSSR()");                    \
   blNormaliseCSSR(cssr, cell, alpha, beta, gamma);                      \
})

#define NormalisePDB(pdb, cell, alpha, beta, gamma)                      \
({                                                                       \
   DEPRECATED("NormalisePDB()","blNormalisePDB()");                      \
   blNormalisePDB(pdb, cell, alpha, beta, gamma);                        \
})

#define ortho(cell, alpha, beta, gamma, amatrx, isw, ncode)              \
({                                                                       \
   DEPRECATED("ortho()","blortho()");                                    \
   blortho(cell, alpha, beta, gamma, amatrx, isw, ncode);                \
})

#define padterm(string, len)                                             \
({                                                                       \
   DEPRECATED("padterm()","blpadterm()");                                \
   blpadterm(string, len);                                               \
})

#define WriteCSSR(fp, cssr, name, title)                                 \
({                                                                       \
   DEPRECATED("WriteCSSR()","blWriteCSSR()");                            \
   blWriteCSSR(fp, cssr, name, title);                                   \
})

*/

/************************************************************************/
/* Renamed functions: fit.h                                             */
/*
#define matfit(x1, x2, rm, n, wt1, column)                               \
({                                                                       \
   DEPRECATED("matfit()","blmatfit()");                                  \
   blmatfit(x1, x2, rm, n, wt1, column);                                 \
})

*/

/************************************************************************/
/* Renamed functions: fsscanf.h                                         */
/* NOTE: Can this be converted to macro? */

/*
#define fsscanf(buffer, format, ...);
*/

/************************************************************************/
/* Renamed functions: general.h                                         */
/*
#define StringToLower(string1, string2)                                  \
({                                                                       \
   DEPRECATED("StringToLower()","blStringToLower()");                    \
   blStringToLower(string1, string2);                                    \
})

#define StringToUpper(string1, string2)                                  \
({                                                                       \
   DEPRECATED("StringToUpper()","blStringToUpper()");                    \
   blStringToUpper(string1, string2);                                    \
})

#define KillLeadSpaces(string)                                           \
({                                                                       \
   DEPRECATED("KillLeadSpaces()","blKillLeadSpaces()");                  \
   blKillLeadSpaces(string);                                             \
})

#define KillLine(fp)                                                     \
({                                                                       \
   DEPRECATED("KillLine()","blKillLine()");                              \
   blKillLine(fp);                                                       \
})

#define SetExtn(File, Ext)                                               \
({                                                                       \
   DEPRECATED("SetExtn()","blSetExtn()");                                \
   blSetExtn(File, Ext);                                                 \
})

#define chindex(string, ch)                                              \
({                                                                       \
   DEPRECATED("chindex()","blchindex()");                                \
   blchindex(string, ch);                                                \
})

#define Word(string1, string2)                                           \
({                                                                       \
   DEPRECATED("Word()","blWord()");                                      \
   blWord(string1, string2);                                             \
})

#define WordN(string1, string2, MaxChar)                                 \
({                                                                       \
   DEPRECATED("WordN()","blWordN()");                                    \
   blWordN(string1, string2, MaxChar);                                   \
})

#define padterm(string, length)                                          \
({                                                                       \
   DEPRECATED("padterm()","blpadterm()");                                \
   blpadterm(string, length);                                            \
})

#define padchar(string, length, ch)                                      \
({                                                                       \
   DEPRECATED("padchar()","blpadchar()");                                \
   blpadchar(string, length, ch);                                        \
})

#define CheckExtn(string, ext)                                           \
({                                                                       \
   DEPRECATED("CheckExtn()","blCheckExtn()");                            \
   blCheckExtn(string, ext);                                             \
})

#define ftostr(str, maxlen, x, precision)                                \
({                                                                       \
   DEPRECATED("ftostr()","blftostr()");                                  \
   blftostr(str, maxlen, x, precision);                                  \
})


#define GetFilestem(filename, stem)                                      \
({                                                                       \
   DEPRECATED("GetFilestem()","blGetFilestem()");                        \
   blGetFilestem(filename, stem);                                        \
})

#define upstrcmp(word1, word2)                                           \
({                                                                       \
   DEPRECATED("upstrcmp()","blupstrcmp()");                              \
   blupstrcmp(word1, word2);                                             \
})

#define upstrncmp(word1, word2, ncomp)                                   \
({                                                                       \
   DEPRECATED("upstrncmp()","blupstrncmp()");                            \
   blupstrncmp(word1, word2, ncomp);                                     \
})

#define GetWord(buffer, word, maxsize)                                   \
({                                                                       \
   DEPRECATED("GetWord()","blGetWord()");                                \
   blGetWord(buffer, word, maxsize);                                     \
})

#define OpenStdFiles(infile, outfile, in, out)                           \
({                                                                       \
   DEPRECATED("OpenStdFiles()","blOpenStdFiles()");                      \
   blOpenStdFiles(infile, outfile, in, out);                             \
})

#define OpenFile(filename, envvar, mode, noenv)                          \
({                                                                       \
   DEPRECATED("OpenFile()","blOpenFile()");                              \
   blOpenFile(filename, envvar, mode, noenv);                            \
})

#define countchar(string, ch)                                            \
({                                                                       \
   DEPRECATED("countchar()","blcountchar()");                            \
   blcountchar(string, ch);                                              \
})

#define fgetsany(fp)                                                     \
({                                                                       \
   DEPRECATED("fgetsany()","blfgetsany()");                              \
   blfgetsany(fp);                                                       \
})

#define strcatalloc(instr, catstr)                                       \
({                                                                       \
   DEPRECATED("strcatalloc()","blstrcatalloc()");                        \
   blstrcatalloc(instr, catstr);                                         \
})


#define StoreString(StringList, string)                                  \
({                                                                       \
   DEPRECATED("StoreString()","blStoreString()");                        \
   blStoreString(StringList, string);                                    \
})

#define InStringList(StringList, string)                                 \
({                                                                       \
   DEPRECATED("InStringList()","blInStringList()");                      \
   blInStringList(StringList, string);                                   \
})

#define FreeStringList(StringList)                                       \
({                                                                       \
   DEPRECATED("FreeStringList()","blFreeStringList()");                  \
   blFreeStringList(StringList);                                         \
})


#define QueryStrStr(string, substring)                                   \
({                                                                       \
   DEPRECATED("QueryStrStr()","blQueryStrStr()");                        \
   blQueryStrStr(string, substring);                                     \
})


#define IndexReal(arrin, indx, n)                                        \
({                                                                       \
   DEPRECATED("IndexReal()","blIndexReal()");                            \
   blIndexReal(arrin, indx, n);                                          \
})


#define OpenOrPipe(filename)                                             \
({                                                                       \
   DEPRECATED("OpenOrPipe()","blOpenOrPipe()");                          \
   blOpenOrPipe(filename);                                               \
})

#define CloseOrPipe(fp)                                                  \
({                                                                       \
   DEPRECATED("CloseOrPipe()","blCloseOrPipe()");                        \
   blCloseOrPipe(fp);                                                    \
})


#define WrapString(in, out, maxlen)                                      \
({                                                                       \
   DEPRECATED("WrapString()","blWrapString()");                          \
   blWrapString(in, out, maxlen);                                        \
})

#define WrapPrint(out, string)                                           \
({                                                                       \
   DEPRECATED("WrapPrint()","blWrapPrint()");                            \
   blWrapPrint(out, string);                                             \
})

#define RightJustify(string)                                             \
({                                                                       \
   DEPRECATED("RightJustify()","blRightJustify()");                      \
   blRightJustify(string);                                               \
})

#define GetWordNC(buffer, word, maxlen)                                  \
({                                                                       \
   DEPRECATED("GetWordNC()","blGetWordNC()");                            \
   blGetWordNC(buffer, word, maxlen);                                    \
})

#define getfield(buffer, start, width, str)                              \
({                                                                       \
   DEPRECATED("getfield()","blgetfield()");                              \
   blgetfield(buffer, start, width, str);                                \
})

*/

/************************************************************************/
/* Renamed functions: help.h                                            */
/*
#define Help(string, HelpFile)                                           \
({                                                                       \
   DEPRECATED("Help()","blHelp()");                                      \
   blHelp(string, HelpFile);                                             \
})

#define DoHelp(string, HelpFile)                                         \
({                                                                       \
   DEPRECATED("DoHelp()","blDoHelp()");                                  \
   blDoHelp(string, HelpFile);                                           \
})

*/

/************************************************************************/
/* Renamed functions: hpgl.h                                            */
/*
#define HPGLInit(filename, AltFont, xmargin, ymargin)                    \
({                                                                       \
   DEPRECATED("HPGLInit()","blHPGLInit()");                              \
   blHPGLInit(filename, AltFont, xmargin, ymargin);                      \
})

#define HPGLPen(num)                                                     \
({                                                                       \
   DEPRECATED("HPGLPen()","blHPGLPen()");                                \
   blHPGLPen(num);                                                       \
})

#define HPGLMove(x, y)                                                   \
({                                                                       \
   DEPRECATED("HPGLMove()","blHPGLMove()");                              \
   blHPGLMove(x, y);                                                     \
})

#define HPGLDraw(x, y)                                                   \
({                                                                       \
   DEPRECATED("HPGLDraw()","blHPGLDraw()");                              \
   blHPGLDraw(x, y);                                                     \
})

#define HPGLSetDash(style)                                               \
({                                                                       \
   DEPRECATED("HPGLSetDash()","blHPGLSetDash()");                        \
   blHPGLSetDash(style);                                                 \
})

#define HPGLFont(font, size)                                             \
({                                                                       \
   DEPRECATED("HPGLFont()","blHPGLFont()");                              \
   blHPGLFont(font, size);                                               \
})

#define HPGLLText(x, y, string)                                          \
({                                                                       \
   DEPRECATED("HPGLLText()","blHPGLLText()");                            \
   blHPGLLText(x, y, string);                                            \
})

#define HPGLCBText(x, y, offset, text)                                   \
({                                                                       \
   DEPRECATED("HPGLCBText()","blHPGLCBText()");                          \
   blHPGLCBText(x, y, offset, text);                                     \
})

#define HPGLROffText(x, y, offset, text)                                 \
({                                                                       \
   DEPRECATED("HPGLROffText()","blHPGLROffText()");                      \
   blHPGLROffText(x, y, offset, text);                                   \
})

#define HPGLLCText(x, y, text)                                           \
({                                                                       \
   DEPRECATED("HPGLLCText()","blHPGLLCText()");                          \
   blHPGLLCText(x, y, text);                                             \
})

#define HPGLCTText(x, y, offset, text)                                   \
({                                                                       \
   DEPRECATED("HPGLCTText()","blHPGLCTText()");                          \
   blHPGLCTText(x, y, offset, text);                                     \
})

#define HPGLVText(x, y, xoff, text, TitleFont, TitleSize, label, LabelFont, LabelSize)\
({                                                                       \
   DEPRECATED("HPGLVText()","blHPGLVText()");                            \
   blHPGLVText(x, y, xoff, text, TitleFont, TitleSize, label, LabelFont, LabelSize);\
})

#define HPGLEnd()                                                        \
({                                                                       \
   DEPRECATED("HPGLEnd()","blHPGLEnd()");                                \
   blHPGLEnd();                                                          \
})

#define HPGLShowText(text, orientation, XBase, YBase)                    \
({                                                                       \
   DEPRECATED("HPGLShowText()","blHPGLShowText()");                      \
   blHPGLShowText(text, orientation, XBase, YBase);                      \
})

*/

/************************************************************************/
/* Renamed functions: matrix.h                                          */
/*
#define MatMult3_33(vecin, matin, vecout)                                \
({                                                                       \
   DEPRECATED("MatMult3_33()","blMatMult3_33()");                        \
   blMatMult3_33(vecin, matin, vecout);                                  \
})

#define MatMult33_33(a, b, out)                                          \
({                                                                       \
   DEPRECATED("MatMult33_33()","blMatMult33_33()");                      \
   blMatMult33_33(a, b, out);                                            \
})

#define invert33(s, ss)                                                  \
({                                                                       \
   DEPRECATED("invert33()","blinvert33()");                              \
   blinvert33(s, ss);                                                    \
})

#define CreateRotMat(direction, angle, matrix)                           \
({                                                                       \
   DEPRECATED("CreateRotMat()","blCreateRotMat()");                      \
   blCreateRotMat(direction, angle, matrix);                             \
})

#define VecDist(a, b, len)                                               \
({                                                                       \
   DEPRECATED("VecDist()","blVecDist()");                                \
   blVecDist(a, b, len);                                                 \
})

*/

/************************************************************************/
/* Renamed functions: parse.h                                           */
/*
#define parse(comline, nkeys, keywords, REALparam, strparam)             \
({                                                                       \
   DEPRECATED("parse()","blparse()");                                    \
   blparse(comline, nkeys, keywords, REALparam, strparam);               \
})

#define mparse(comline, nkeys, keywords, REALparam, strparam, nparams)   \
({                                                                       \
   DEPRECATED("mparse()","blmparse()");                                  \
   blmparse(comline, nkeys, keywords, REALparam, strparam, nparams);     \
})

#define match(comstring, string2, nletters)                              \
({                                                                       \
   DEPRECATED("match()","blmatch()");                                    \
   blmatch(comstring, string2, nletters);                                \
})

#define GetString(command, strparam)                                     \
({                                                                       \
   DEPRECATED("GetString()","blGetString()");                            \
   blGetString(command, strparam);                                       \
})

#define GetParam( command, value, nletters)                              \
({                                                                       \
   DEPRECATED("GetParam()","blGetParam()");                              \
   blGetParam( command, value, nletters);                                \
})
*/

/************************************************************************/
/* Renamed functions: plotting.h                                        */
/*
#define AMInitPlot(filename, title, dest, OutXSize, OutYSize, OutXOff, OutYOff, AltFont, xmargin, ymargin, DataXMin, DataYMin, DataXMax, DataYMax)\
({                                                                       \
   DEPRECATED("AMInitPlot()","blAMInitPlot()");                          \
   blAMInitPlot(filename, title, dest, OutXSize, OutYSize, OutXOff, OutYOff, AltFont, xmargin, ymargin, DataXMin, DataYMin, DataXMax, DataYMax);\
})

#define AMSetPen(dest, pen)                                              \
({                                                                       \
   DEPRECATED("AMSetPen()","blAMSetPen()");                              \
   blAMSetPen(dest, pen);                                                \
})

#define AMMove(dest, x, y)                                               \
({                                                                       \
   DEPRECATED("AMMove()","blAMMove()");                                  \
   blAMMove(dest, x, y);                                                 \
})

#define AMDraw(dest, x, y)                                               \
({                                                                       \
   DEPRECATED("AMDraw()","blAMDraw()");                                  \
   blAMDraw(dest, x, y);                                                 \
})

#define AMSetLineStyle(dest, style)                                      \
({                                                                       \
   DEPRECATED("AMSetLineStyle()","blAMSetLineStyle()");                  \
   blAMSetLineStyle(dest, style);                                        \
})

#define AMEndLine(dest)                                                  \
({                                                                       \
   DEPRECATED("AMEndLine()","blAMEndLine()");                            \
   blAMEndLine(dest);                                                    \
})

#define AMSetFont(dest, PSFontName, FontSize)                            \
({                                                                       \
   DEPRECATED("AMSetFont()","blAMSetFont()");                            \
   blAMSetFont(dest, PSFontName, FontSize);                              \
})

#define AMText(dest, x, y, text)                                         \
({                                                                       \
   DEPRECATED("AMText()","blAMText()");                                  \
   blAMText(dest, x, y, text);                                           \
})

#define AMCBText(dest, x, y, text)                                       \
({                                                                       \
   DEPRECATED("AMCBText()","blAMCBText()");                              \
   blAMCBText(dest, x, y, text);                                         \
})

#define AMRText(dest, x, y, offset, text)                                \
({                                                                       \
   DEPRECATED("AMRText()","blAMRText()");                                \
   blAMRText(dest, x, y, offset, text);                                  \
})

#define AMLCText(dest, x, y, text)                                       \
({                                                                       \
   DEPRECATED("AMLCText()","blAMLCText()");                              \
   blAMLCText(dest, x, y, text);                                         \
})

#define AMCTText(dest, x, y, CTOffset, text)                             \
({                                                                       \
   DEPRECATED("AMCTText()","blAMCTText()");                              \
   blAMCTText(dest, x, y, CTOffset, text);                               \
})

#define AMEndPlot(dest)                                                  \
({                                                                       \
   DEPRECATED("AMEndPlot()","blAMEndPlot()");                            \
   blAMEndPlot(dest);                                                    \
})

#define PS2HPGLFont(font)                                                \
({                                                                       \
   DEPRECATED("PS2HPGLFont()","blPS2HPGLFont()");                        \
   blPS2HPGLFont(font);                                                  \
})

#define SimplifyText(string)                                             \
({                                                                       \
   DEPRECATED("SimplifyText()","blSimplifyText()");                      \
   blSimplifyText(string);                                               \
})

*/


/************************************************************************/
/* Renamed functions: ps.h                                              */
/*
#define PSInit(FName, creator, AltFont)                                  \
({                                                                       \
   DEPRECATED("PSInit()","blPSInit()");                                  \
   blPSInit(FName, creator, AltFont);                                    \
})

#define PSThick(thickness)                                               \
({                                                                       \
   DEPRECATED("PSThick()","blPSThick()");                                \
   blPSThick(thickness);                                                 \
})

#define PSMove(X, Y)                                                     \
({                                                                       \
   DEPRECATED("PSMove()","blPSMove()");                                  \
   blPSMove(X, Y);                                                       \
})

#define PSDraw(X, Y)                                                     \
({                                                                       \
   DEPRECATED("PSDraw()","blPSDraw()");                                  \
   blPSDraw(X, Y);                                                       \
})

#define PSSetDash(linepatt)                                              \
({                                                                       \
   DEPRECATED("PSSetDash()","blPSSetDash()");                            \
   blPSSetDash(linepatt);                                                \
})

#define PSClearDash()                                                    \
({                                                                       \
   DEPRECATED("PSClearDash()","blPSClearDash()");                        \
   blPSClearDash();                                                      \
})

#define PSStroke()                                                       \
({                                                                       \
   DEPRECATED("PSStroke()","blPSStroke()");                              \
   blPSStroke();                                                         \
})

#define PSFont(fontname, size)                                           \
({                                                                       \
   DEPRECATED("PSFont()","blPSFont()");                                  \
   blPSFont(fontname, size);                                             \
})

#define PSLText(X, Y, label)                                             \
({                                                                       \
   DEPRECATED("PSLText()","blPSLText()");                                \
   blPSLText(X, Y, label);                                               \
})

#define PSCBText(X, Y, Offset, label)                                    \
({                                                                       \
   DEPRECATED("PSCBText()","blPSCBText()");                              \
   blPSCBText(X, Y, Offset, label);                                      \
})

#define PSROffText(X, Y, offset, label)                                  \
({                                                                       \
   DEPRECATED("PSROffText()","blPSROffText()");                          \
   blPSROffText(X, Y, offset, label);                                    \
})

#define PSLCText(X, Y, label)                                            \
({                                                                       \
   DEPRECATED("PSLCText()","blPSLCText()");                              \
   blPSLCText(X, Y, label);                                              \
})

#define PSCTText(X, Y, Offset, label)                                    \
({                                                                       \
   DEPRECATED("PSCTText()","blPSCTText()");                              \
   blPSCTText(X, Y, Offset, label);                                      \
})

#define PSVText(x, y, xoff, text, font, size, label, lfont, lsize)       \
({                                                                       \
   DEPRECATED("PSVText()","blPSVText()");                                \
   blPSVText(x, y, xoff, text, font, size, label, lfont, lsize);         \
})

#define PSShowText(text)                                                 \
({                                                                       \
   DEPRECATED("PSShowText()","blPSShowText()");                          \
   blPSShowText(text);                                                   \
})

#define PSEnd()                                                          \
({                                                                       \
   DEPRECATED("PSEnd()","blPSEnd()");                                    \
   blPSEnd();                                                            \
})

#define PSCorrectCase(font)                                              \
({                                                                       \
   DEPRECATED("PSCorrectCase()","blPSCorrectCase()");                    \
   blPSCorrectCase(font);                                                \
})

*/


/************************************************************************/
/* Renamed functions: seq.h                                             */
/*
#define throne(three)                                                    \
({                                                                       \
   DEPRECATED("throne()","blthrone()");                                  \
   blthrone(three);                                                      \
})

#define thronex(three)                                                   \
({                                                                       \
   DEPRECATED("thronex()","blthronex()");                                \
   blthronex(three);                                                     \
})

#define onethr(one)                                                      \
({                                                                       \
   DEPRECATED("onethr()","blonethr()");                                  \
   blonethr(one);                                                        \
})

#define DoPDB2Seq(pdb, DoAsxGlx, ProtOnly, NoX)                          \
({                                                                       \
   DEPRECATED("DoPDB2Seq()","blDoPDB2Seq()");                            \
   blDoPDB2Seq(pdb, DoAsxGlx, ProtOnly, NoX);                            \
})

#define SplitSeq(LinearSeq, seqs)                                        \
({                                                                       \
   DEPRECATED("SplitSeq()","blSplitSeq()");                              \
   blSplitSeq(LinearSeq, seqs);                                          \
})

#define ReadSimplePIR(fp, maxres, seqs)                                  \
({                                                                       \
   DEPRECATED("ReadSimplePIR()","blReadSimplePIR()");                    \
   blReadSimplePIR(fp, maxres, seqs);                                    \
})

#define ReadPIR(fp, DoInsert, seqs, maxchain, seqinfo, punct, error)     \
({                                                                       \
   DEPRECATED("ReadPIR()","blReadPIR()");                                \
   blReadPIR(fp, DoInsert, seqs, maxchain, seqinfo, punct, error);       \
})

#define ReadRawPIR(fp, seqs, maxchain, upcase, seqinfo, error)           \
({                                                                       \
   DEPRECATED("ReadRawPIR()","blReadRawPIR()");                          \
   blReadRawPIR(fp, seqs, maxchain, upcase, seqinfo, error);             \
})

#define align(seq1, length1, seq2, length2, verbose, identity, penalty, align1, align2, align_len)\
({                                                                       \
   DEPRECATED("align()","blalign()");                                    \
   blalign(seq1, length1, seq2, length2, verbose, identity, penalty, align1, align2, align_len);\
})

#define affinealign(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len)\
({                                                                       \
   DEPRECATED("affinealign()","blaffinealign()");                        \
   blaffinealign(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len);\
})

#define CalcMDMScore(resa, resb)                                         \
({                                                                       \
   DEPRECATED("CalcMDMScore()","blCalcMDMScore()");                      \
   blCalcMDMScore(resa, resb);                                           \
})

#define affinealignuc(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len)\
({                                                                       \
   DEPRECATED("affinealignuc()","blaffinealignuc()");                    \
   blaffinealignuc(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len);\
})

#define CalcMDMScoreUC(resa, resb)                                       \
({                                                                       \
   DEPRECATED("CalcMDMScoreUC()","blCalcMDMScoreUC()");                  \
   blCalcMDMScoreUC(resa, resb);                                         \
})

#define ReadMDM(mdmfile)                                                 \
({                                                                       \
   DEPRECATED("ReadMDM()","blReadMDM()");                                \
   blReadMDM(mdmfile);                                                   \
})

#define ZeroMDM()                                                        \
({                                                                       \
   DEPRECATED("ZeroMDM()","blZeroMDM()");                                \
   blZeroMDM();                                                          \
})

#define DNAtoAA(dna)                                                     \
({                                                                       \
   DEPRECATED("DNAtoAA()","blDNAtoAA()");                                \
   blDNAtoAA(dna);                                                       \
})

#define TrueSeqLen(sequence)                                             \
({                                                                       \
   DEPRECATED("TrueSeqLen()","blTrueSeqLen()");                          \
   blTrueSeqLen(sequence);                                               \
})

#define KnownSeqLen(sequence)                                            \
({                                                                       \
   DEPRECATED("KnownSeqLen()","blKnownSeqLen()");                        \
   blKnownSeqLen(sequence);                                              \
})

#define NumericReadMDM(mdmfile)                                          \
({                                                                       \
   DEPRECATED("NumericReadMDM()","blNumericReadMDM()");                  \
   blNumericReadMDM(mdmfile);                                            \
})

#define NumericCalcMDMScore(resa, resb)                                  \
({                                                                       \
   DEPRECATED("NumericCalcMDMScore()","blNumericCalcMDMScore()");        \
   blNumericCalcMDMScore(resa, resb);                                    \
})

#define NumericAffineAlign(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len)\
({                                                                       \
   DEPRECATED("NumericAffineAlign()","blNumericAffineAlign()");          \
   blNumericAffineAlign(seq1, length1, seq2, length2, verbose, identity, penalty, penext, align1, align2, align_len);\
})

*/
 



/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#endif

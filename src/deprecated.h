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
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#endif

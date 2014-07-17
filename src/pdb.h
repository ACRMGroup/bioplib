/************************************************************************/
/**

   \file       pdb.h
   
   \version    V1.58
   \date       17.07.14
   \brief      Include file for pdb routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin, UCL, Reading 1993-2014
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
-  V1.0  04.11.88 Original
-  V1.1  22.03.90 Added Secondary structure routines
-  V1.2  28.03.90 Corrected field widths for V1.2 of ReadPDB
-  V1.3  04.05.90 Added clear_pdb()
-  V1.4  19.06.90 Changed SEC structure to correct chain and ins widths
-  V1.5  19.07.90 Added INITINDEX macro
-  V1.6  09.09.91 Added define so won't screw up if included twice
-  V1.7  22.09.91 Altered character sizes for alignment
-  V1.8  10.06.93 Changed to use REAL rather than float. Changed 
                  order within structure
-  V1.9  22.02.94 Added MAXSTDAA and MAXATINRES definitions
-  V1.10 01.03.94 Added stuff for ResolPDB. Removed INIT_INDEX().
                  Added DISULPHIDE definition.
                  Added HADDINFO definition.
-  V1.11 18.03.94 Added prototypes for ReadPDBOccRank() and
                  ReadPDBAtomsOccRank()
                  Added gPDBPartialOcc
-  V1.12 23.05.94 Added FindNextChainPDB() prototype
-  V1.13 24.08.94 Added OpenPGPFile() prototype. Added prototypes for
                  new version of FixPDB(). Added CTER styles
-  V1.14 03.10.94 Added FindCofGPDBRange(), FindCofGPDBSCRange(),
                  ReadPDBALL()
-  V1.15 05.10.94 Changed KillSidechain()
-  V1.16 11.01.94 Added StripHPDB()
-  V1.17 06.03.95 doReadPDB() is now defined here rather than static
-  V1.18 17.07.95 ParseResSpec() is now a BOOL
-  V1.19 24.07.95 Added FNam2PDB(), TermPDB()
-  V1.20 25.07.95 Added GetPDBChainLabels()
-  V1.21 08.08.95 Added FindResidueSpec() and FindNextResidue()
-  V1.22 12.10.95 Added DupePDB(), CopyPDBCoords(), CalcCellTrans(),
                  GetCrystPDB(), WriteCrystPDB()
-  V1.23 10.01.96 Added ExtractZonePDB()
-  V1.24 08.02.96 Added FindResidue()
-  V1.25 14.03.96 Added FitCaPDB(), FindAtomInRes()
-  V1.26 18.06.96 Added InPDBZone() and ZONE_MODE_*. Modified prototype
                  for FindZonePDB()
-  V1.27 23.07.96 Added AtomNameMatch() and LegalAtomSpec()
-  V1.28 12.08.96 Added RepOneSChain() and EndRepSChain()
-  V1.29 19.09.96 Added InPDBZoneSpec()
-  V1.30 14.10.96 Added ReadSeqresPDB();
-  V1.31 16.10.96 Added SelectCaPDB()
-  V1.32 18.08.98 Changed SEC to SECSTRUC 'cos of conflict in SunOS
                  Also defines SEC macro if not defined to warn you to
                  change your code!
-  V1.33 28.04.99 Added GetExptl()
-  V1.34 15.02.01 Added atnam_raw[] to PDB
                  Added WriteGromosPDB(), WriteGromosPDBRecord(),
                        AtomNameRawMatch()
-  V1.35 12.12.01 Added FitNCaCPDB()
-  V1.36 30.05.02 Changed PDB field from 'junk' to 'record_type'
                  Added the WholePDB routines and definition
-  V1.37 03.06.05 Added altpos to PDB.
                  Added altpos and atnam_raw to CLEAR_PDB
-  V1.38 22.09.05 Added WritePDBRecordAtnam()
-  V1.39 29.09.05 Added ParseResSpecNoUpper() and DoParseResSpec()  By: TL
-  V1.40 04.01.06 Added AddCBtiGly(), AddCBtoAllGly(), 
                  StripGlyCB()      By: ACRM
-  V1.41 25.01.06 Added RemoveAlternates()
-  V1.42 08.11.07 Added BuildAtomNeighbourPDBList()
                        FindAtomWildcardInRes()
                        DupeResiduePDB()
-  V1.43 30.04.08 Added StripWatersPDB() and ISWATER() macro
-  V1.44 01.06.09 Added extras field to PDB structure
-  V1.45 24.11.09 Added PDBSTRUCT, PDBCHAIN, PDBRESIDUE
                  AllocPDBStructure(), FindNextChain(),
                  FreePDBStructure()
-  V1.46 26.10.11 Added FindHetatmResidueSpec() and FindHetatmResidue()
-  V1.47 12.12.11 Added GetExptlOld()
                  Added ResportStructureType()
                  Added new STRUCTURE_TYPE_* defines
-  V1.48 04.02.14 Added CHAINMATCH macro. By: CTP
-  V1.49 24.02.14 Added BiopFindResidue(), BiopFindHetatmResidue() and 
                        BiopInPDBZone(). By: CTP
-  V1.50 20.03.14 Added blFindZonePDB(). By: CTP
-  V1.51 25.03.14 Added blGetPDBChainLabels(). By: CTP
-  V1.52 22.04.14 Added CheckFileFormatPDBML(FILE *fp). By: CTP
-  V1.53 07.05.14 Added deprecated.h and removed definitions for 
                  deprecated funtions: FindHetatmResidue(), FindResidue(),
                  InPDBZone() and FindZonePDB(). By: CTP
-  V1.54 02.06.14 Added WritePDBML() By: CTP
-  V1.55 09.06.14 Added gPDBXML flag By: CTP
-  V1.56 21.06.14 Added gPDBXMLForce flag and updated functions:
                  blWritePDB(), blWriteAsPDB(), blWriteAsPDBML(), 
                  blFormatCheckWritePDB(), blWriteWholePDB(), 
                  blWriteWholePDBHeader() and blWriteWholePDBTrailer().
                  By: CTP
-  V1.57 07.07.14 Rename functions with 'bl' prefix. By: CTP
-  V1.58 17.07.14 Added access and radius to PDB structure. Also
                  added CREATEPDBEXTRAS() and FREEPDBEXTRAS()  By: ACRM
-  V1.59 17.07.14 Added blSetElementSymbolFromAtomName() By: CTP

*************************************************************************/
#ifndef _PDB_H
#define _PDB_H

#include <stdio.h>
#include <string.h>

#include "MathType.h"
#include "SysDefs.h"
#include "general.h"
#include "deprecated.h"

#define MAXSTDAA    21       /* Number of standard amino acids (w/ PCA)*/
#define MAXATINRES  14       /* Max number of atoms in a standard aa   */

/* This is our main PDB structure used for the PDB linked lists.

   The 'extras' field is used for flags or to attach another
   structure or array to each PDB record. For example:

   typedef struct
   {
      REAL angle;
      BOOL flag;
   }  EXTRAS;
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->extras = (APTR)malloc(sizeof(EXTRAS)))==NULL)
         return(FALSE);
   }

   for(p=pdb; p!=NULL; NEXT(p))
   {
      ((EXTRAS *)p->extras)->flag = FALSE;
      ((EXTRAS *)p->extras)->angle = (REAL)0.0;
   }
*/
typedef struct pdb_entry
{
   REAL x,y,z,occ,bval,access,radius;
   APTR extras;
   struct pdb_entry *next;
   int  atnum;
   int  resnum;
   char record_type[8];
   char atnam[8];
   char atnam_raw[8];
   char resnam[8];
   char insert[8];
   char chain[8];
   char altpos;
}  PDB;

typedef struct pdbresidue
{
   struct pdbresidue *next,  *prev;
   PDB               *start, *stop;
   APTR              *extras;
   int               resnum;
   char              chain[8];
   char              insert[8];
   char              resnam[8];
   char              resid[8];
} PDBRESIDUE;

typedef struct pdbchain
{
   struct pdbchain *next,  *prev;
   PDB             *start, *stop;
   PDBRESIDUE      *residues;
   APTR            *extras;
   char            chain[8];
} PDBCHAIN;

typedef struct
{
   PDB      *pdb;
   PDBCHAIN *chains;
   APTR     *extras;
} PDBSTRUCT;


#define SELECT(x,w) (x) = (char *)malloc(5 * sizeof(char)); \
                    if((x) != NULL) strncpy((x),(w),5)

typedef struct sec_entry
{
   struct sec_entry *next;
   char chain1[8];
   char ins1[8];
   char chain2[8];
   char ins2[8];
   int  res1;
   int  res2;
   char type;
}  SECSTRUC;

typedef struct _wholepdb
{
   PDB        *pdb;
   STRINGLIST *header;
   STRINGLIST *trailer;
   int        natoms;
}  WHOLEPDB;

/* This is designed to cause an error message which prints this line
   It has been tested with gcc and Irix cc and does as required in
   both cases
*/
#ifndef SEC
#   define SEC (The_type_SEC_is_now_called_SECSTRUC_You_must_change_your_code *)
#endif

typedef struct _disulphide
{
   struct _disulphide *next;
   int                res1,
                      res2;
   char               chain1[8],
                      chain2[8],
                      insert1[8],
                      insert2[8];
}  DISULPHIDE;

typedef struct
{
   int   Total,      /* Total hydrogens                                 */
         T1,         /* Type 1 C-H's                                    */
         T2,         /* Type 2 C-H2's                                   */
         T3,         /* Type 3 C-H3's                                   */
         T4,         /* Type 4 sp2 C-H's,>N-H                           */
         T5;         /* Type 5 O-H's =N-H's                             */
}  HADDINFO;

#define CLEAR_PDB(p) strcpy(p->record_type,"      "); \
                     p->atnum=0; \
                     strcpy(p->atnam,"    "); \
                     strcpy(p->atnam_raw,"    "); \
                     strcpy(p->resnam,"    "); \
                     p->resnum=0; \
                     strcpy(p->insert," "); \
                     strcpy(p->chain," "); \
                     p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                     p->altpos = ' '; \
                     p->occ = 0.0; p->bval = 0.0; \
                     p->next = NULL

#define ISWATER(z)   (!strncmp((z)->resnam,"HOH",3) || \
                      !strncmp((z)->resnam,"OH2",3) || \
                      !strncmp((z)->resnam,"OHH",3) || \
                      !strncmp((z)->resnam,"DOD",3) || \
                      !strncmp((z)->resnam,"OD2",3) || \
                      !strncmp((z)->resnam,"ODD",3) || \
                      !strncmp((z)->resnam,"WAT",3))

#define CHAINMATCH(chain1,chain2) !strcmp(chain1,chain2)

/* Called as CREATEPDBEXTRAS(pdb, EXTRATYPE)                            */
#define CREATEPDBEXTRAS(x, y)                                   \
   {  PDB *_cpe_p;                                              \
      for(_cpe_p=(x); _cpe_p!=NULL; _cpe_p=_cpe_p->next){       \
         _cpe_p->extras = (APTR)malloc(sizeof(y));              \
      }                                                         \
   }

/* Called as FREEPDBEXTRAS(pdb)                                         */
#define FREEPDBEXTRAS(x)                                        \
   {  PDB *_fpe_p;                                              \
      for(_fpe_p=(x); _fpe_p!=NULL; _fpe_p=_fpe_p->next){       \
         if(_fpe_p->extras != NULL){                            \
            free(_fpe_p->extras);                               \
            _fpe_p->extras = NULL;                              \
         }                                                      \
      }                                                         \
   }
      


/* These are the types returned by ResolPDB()                          */
#define STRUCTURE_TYPE_UNKNOWN   0
#define STRUCTURE_TYPE_XTAL      1
#define STRUCTURE_TYPE_NMR       2
#define STRUCTURE_TYPE_MODEL     3
#define STRUCTURE_TYPE_ELECTDIFF 4
#define STRUCTURE_TYPE_FIBER     5
#define STRUCTURE_TYPE_SSNMR     6
#define STRUCTURE_TYPE_NEUTRON   7
#define STRUCTURE_TYPE_EM        8
#define STRUCTURE_TYPE_SOLSCAT   9
#define STRUCTURE_TYPE_IR       10
#define STRUCTURE_TYPE_POWDER   11
#define STRUCTURE_TYPE_FRET     12

/* These are the styles used by FixCterPDB()                            */
#define CTER_STYLE_STD         0
#define CTER_STYLE_GROMOS      1
#define CTER_STYLE_CHARMM      2

/* Return flags from GetCrystPDB()                                      */
#define XTAL_DATA_CRYST        0x0001
#define XTAL_DATA_ORIGX        0x0002
#define XTAL_DATA_SCALE        0x0004

/* Modes for FindZonePDB()                                              */
#define ZONE_MODE_RESNUM       0
#define ZONE_MODE_SEQUENTIAL   1

/* Modes and macros for gPDBXMLForce                                    */
#define FORCEXML_NOFORCE       0
#define FORCEXML_PDB           1
#define FORCEXML_XML           2

/* For forcing writing in PDB or XML format                             */
#define FORCEPDB gPDBXMLForce = FORCEXML_PDB
#define FORCEXML gPDBXMLForce = FORCEXML_XML

/************************************************************************/
/* Globals
*/
#ifdef RSC_MAIN
   char gRSCError[80];
#else
   extern char gRSCError[80];
#endif

#ifdef READPDB_MAIN
   BOOL gPDBPartialOcc;
   BOOL gPDBMultiNMR;
   BOOL gPDBXML = FALSE;
#else
   extern BOOL gPDBPartialOcc;
   extern BOOL gPDBMultiNMR;
   extern BOOL gPDBXML;
#endif

#ifdef WRITEPDB_MAIN
   int gPDBXMLForce = FORCEXML_NOFORCE;
#else
   extern int gPDBXMLForce;
#endif

/************************************************************************/
/* Prototypes
*/

PDB *blReadPDB(FILE *fp, int *natom);
PDB *blReadPDBAll(FILE *fp, int *natom);
PDB *blReadPDBAtoms(FILE *fp, int *natom);
PDB *blReadPDBOccRank(FILE *fp, int *natom, int OccRank);
PDB *blReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank);
PDB *blDoReadPDB(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
                 int ModelNum);
PDB *blDoReadPDBML(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
                   int ModelNum);
BOOL blCheckFileFormatPDBML(FILE *fp);

BOOL blWritePDB(FILE *fp, PDB  *pdb);
void blWriteAsPDB(FILE *fp, PDB  *pdb);
void blWriteAsPDBML(FILE *fp, PDB  *pdb);
BOOL blFormatCheckWritePDB(PDB *pdb);
BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb);
void blWriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb);
void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb);

void blWritePDBRecord(FILE *fp, PDB *pdb);
void blWritePDBRecordAtnam(FILE *fp, PDB  *pdb);
void blWriteGromosPDB(FILE *fp, PDB *pdb);
void blWriteGromosPDBRecord(FILE *fp, PDB *pdb);
void blGetCofGPDB(PDB   *pdb, VEC3F *cg);
void blGetCofGPDBRange(PDB *start, PDB *stop, VEC3F *cg);
void blGetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg);
void blOriginPDB(PDB *pdb);
void blRotatePDB(PDB  *pdb, REAL rm[3][3]);
void blTranslatePDB(PDB   *pdb, VEC3F tvect);
BOOL blFitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL blFitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL blFitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL blFitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
REAL blCalcRMSPDB(PDB *pdb1, PDB *pdb2);
int blGetPDBCoor(PDB *pdb, COOR **coor);
BOOL blFindZonePDB(PDB *pdb, int start, char *startinsert, int stop, 
                   char *stopinsert, char *chain, int mode, 
                   PDB **pdb_start, PDB **pdb_stop);
int blHAddPDB(FILE *fp, PDB *pdb);
int blReadPGP(FILE *fp);
FILE *blOpenPGPFile(char *pgpfile, BOOL AllHyd);
PDB *blSelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom);
PDB *blStripHPDB(PDB *pdbin, int *natom);
SECSTRUC *blReadSecPDB(FILE *fp, int *nsec);
void blRenumAtomsPDB(PDB *pdb);
PDB *blFindEndPDB(PDB *start);
PDB *blFixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum);
PDB *blShuffleResPDB(PDB *start, PDB *end, BOOL Pad);
BOOL blGetAtomTypes(char *resnam, char **AtomTypes);
PDB *blKillPDB(PDB *pdb, PDB *prev);
void blCopyPDB(PDB *out, PDB *in);
BOOL blMovePDB(PDB *move, PDB **from, PDB **to);
PDB *blAppendPDB(PDB *first, PDB *second);
PDB *blShuffleBB(PDB *pdb);
REAL blCalcChi(PDB *pdb, int type);
PDB *blGetPDBByN(PDB *pdb, int n);
void blSetChi(PDB *pdb, PDB *next, REAL chi, int type);
BOOL blKillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB);
void blSetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum,   
                 char *insert, char *chain);
void blApplyMatrixPDB(PDB *pdb, REAL matrix[3][3]);
BOOL blGetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor, 
                 int *StrucType);
BOOL blGetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
BOOL blGetExptlOld(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
char *blReportStructureType(int type);
PDB **blIndexPDB(PDB *pdb, int *natom);
DISULPHIDE *blReadDisulphidesPDB(FILE *fp, BOOL *error);
BOOL blParseResSpec(char *spec, char *chain, int *resnum, char *insert);
BOOL blParseResSpecNoUpper(char *spec, char *chain, int *resnum, char *insert);
BOOL blDoParseResSpec(char *spec, char *chain, int *resnum, char *insert, 
                      BOOL uppercaseresspec);
BOOL blRepSChain(PDB *pdb, char *sequence, char *ChiTable, char *RefCoords);
PDB *blFindNextChainPDB(PDB *pdb);
BOOL blFixCterPDB(PDB *pdb, int style);
BOOL blCalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p);
int blCalcTetraHCoords(PDB *nter, COOR *coor);
int blAddNTerHs(PDB **ppdb, BOOL Charmm);
char *blFNam2PDB(char *filename);
PDB *blTermPDB(PDB *pdb, int length);
char *GetPDBChainLabels(PDB *pdb); /* Deprecated function */
char **blGetPDBChainLabels(PDB *pdb, int *nchains);
PDB *blFindHetatmResidueSpec(PDB *pdb, char *resspec);
PDB *blFindResidueSpec(PDB *pdb, char *resspec);
PDB *blFindNextResidue(PDB *pdb);
PDB *blDupePDB(PDB *in);
BOOL blCopyPDBCoords(PDB *out, PDB *in);
void blCalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                     VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans);
int blGetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles,
                  char *spacegroup,
                  REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
void blWriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                     char *spacegroup,
                     REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
PDB *blExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                      char *chain2, int resnum2, char *insert2);
PDB *blFindResidue(PDB *pdb, char *chain, int resnum, char *insert);
PDB *blFindHetatmResidue(PDB *pdb, char *chain, int resnum, char *insert);
PDB *blFindAtomInRes(PDB *pdb, char *atnam);
BOOL blInPDBZone(PDB *p, char *chain, int resnum1, char *insert1, 
                 int resnum2, char *insert2);
BOOL blInPDBZoneSpec(PDB *p, char *resspec1, char *resspec2);
BOOL blAtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL blAtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL blLegalAtomSpec(char *spec);
BOOL blRepOneSChain(PDB *pdb, char *ResSpec, char aa, char *ChiTable,
                    char *RefCoords);
void blEndRepSChain(void);
char **blReadSeqresPDB(FILE *fp, int *nchains);
PDB *blSelectCaPDB(PDB *pdb);
char *blFixAtomName(char *name, REAL occup);

void blFreeWholePDB(WHOLEPDB *wpdb);
WHOLEPDB *blReadWholePDB(FILE *fpin);
WHOLEPDB *blReadWholePDBAtoms(FILE *fpin);
BOOL blAddCBtoGly(PDB *pdb);
BOOL blAddCBtoAllGly(PDB *pdb);
PDB *blStripGlyCB(PDB *pdb);
PDB *blRemoveAlternates(PDB *pdb);
PDB *blBuildAtomNeighbourPDBList(PDB *pdb, PDB *pRes, REAL NeighbDist);
PDB *blFindAtomWildcardInRes(PDB *pdb, char *pattern);
PDB *blDupeResiduePDB(PDB *in);
PDB *blStripWatersPDB(PDB *pdbin, int *natom);
PDBSTRUCT *blAllocPDBStructure(PDB *pdb);
PDB *blFindNextChain(PDB *pdb);
void blFreePDBStructure(PDBSTRUCT *pdbstruct);
void blSetElementSymbolFromAtomName(char *element, char * atom_name);
#endif

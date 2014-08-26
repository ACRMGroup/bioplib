/************************************************************************/
/**

   \file       deprecated.h
   
   \version    v1.3
   \date       14.08.14
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

   Allows code utilising BiopLib to compile without modification to use 
   function names with the 'bl' prefix.(The 'bl' prefix was introduced in 
   July 2014.) 

   Contains prototypes for all deprecated functions.

   Contains the DEPRECATED() macro which prints a warning message when a 
   call is made to a deprecated function.

   Allows use of deprecated functions that took a single character for the
   PDB chain identifier or insert value.

   Used by the following header files:

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


\par Documentation
   Doxygen is not set to document deprecated functions. To document 
   deprecated functions, add "deprecated" to the ENABLED SECTIONS tag in
   the Doxygen config file:
   bioplib/doc/doxygen/Doxyfile


**************************************************************************


   Usage:
   ======

   The deprecated.h header file is included at the end of a header file 
   with deprecated functions. This is to ensure that data structures are
   defined _before_ the prototypes that use them are declared in
   deprecated.h. For example, the PDB data structure must be defined 
   before declaring the prototypes for the functions returning a PDB 
   linked list.

   Prototypes for the deprecated functions for each header file are 
   enclosed within compiler directives to ensure that functions are not 
   declared using undefined data structures.


**************************************************************************

   Revision History:
   =================
   
-  V1.0  07.05.14 Original By: CTP
-  V1.1  07.07.14 Rename functions with 'bl' prefix. By: CTP
-  V1.2  31.07.14 Rewrite of deprecation system. Moved deprecated function
                  code to deprecated.c  By: CTP 
-  V1.3  14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP

*************************************************************************/
#ifndef _DEPRECATED_H
#define _DEPRECATED_H

/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"


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
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

# endif

/************************************************************************/
/** \cond deprecated                                                    */
/************************************************************************/
/* Deprecated functions: general.h                                      */

#ifdef _GENERAL_H_DEPRECATED
#undef _GENERAL_H_DEPRECATED

void StringToLower(char *string1, char *string2);
void StringToUpper(char *string1, char *string2);
char *KillLeadSpaces(char *string);
void KillLine(FILE *fp);
void SetExtn(char *File, char *Ext);
int chindex(char *string, char ch);
void Word(char *string1, char *string2);
void WordN(char *string1, char *string2, int  MaxChar);
void padterm(char *string, int length); /* defined in cssr.h */
void padchar(char *string, int length, char ch);
BOOL CheckExtn(char *string, char *ext);
char *ftostr(char *str, int maxlen, REAL x, int precision);

void GetFilestem(char *filename, char *stem);
int upstrcmp(char *word1, char *word2);
int upstrncmp(char *word1, char *word2, int ncomp);
char *GetWord(char *buffer, char *word, int maxsize);
BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out);
FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv);
int countchar(char *string, char ch);
char *fgetsany(FILE *fp);
char *strcatalloc(char *instr, char *catstr);

STRINGLIST *StoreString(STRINGLIST *StringList, char *string);
BOOL InStringList(STRINGLIST *StringList, char *string);
void FreeStringList(STRINGLIST *StringList);

char *QueryStrStr(char *string, char *substring);

void IndexReal(REAL *arrin, int *indx, int n);

FILE *OpenOrPipe(char *filename);
int CloseOrPipe(FILE *fp);

BOOL WrapString(char *in, char *out, int maxlen);
BOOL WrapPrint(FILE *out, char *string);
void RightJustify(char *string);
char *GetWordNC(char *buffer, char *word, int maxlen);
void getfield(char *buffer, int start, int width, char *str);

#endif

/************************************************************************/
/* Deprecated functions: pdb.h                                          */

#ifdef _PDB_H_DEPRECATED
#undef _PDB_H_DEPRECATED

char *GetPDBChainLabels(PDB *pdb);
PDB *ReadPDB(FILE *fp, int *natom);
PDB *ReadPDBAll(FILE *fp, int *natom);
PDB *ReadPDBAtoms(FILE *fp, int *natom);
PDB *ReadPDBOccRank(FILE *fp, int *natom, int OccRank);
PDB *ReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank);
PDB *doReadPDB(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
                 int ModelNum);
PDB *doReadPDBML(FILE *fp, int  *natom, BOOL AllAtoms, int OccRank, 
                   int ModelNum);
BOOL CheckFileFormatPDBML(FILE *fp);

void WritePDB(FILE *fp, PDB  *pdb);
void WriteAsPDB(FILE *fp, PDB  *pdb);
void WriteAsPDBML(FILE *fp, PDB  *pdb);
BOOL FormatCheckWritePDB(PDB *pdb);
void WriteWholePDB(FILE *fp, WHOLEPDB *wpdb);
void WriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb);
void WriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb);

void WritePDBRecord(FILE *fp, PDB *pdb);
void WritePDBRecordAtnam(FILE *fp, PDB  *pdb);
void WriteGromosPDB(FILE *fp, PDB *pdb);
void WriteGromosPDBRecord(FILE *fp, PDB *pdb);
void GetCofGPDB(PDB   *pdb, VEC3F *cg);
void GetCofGPDBRange(PDB *start, PDB *stop, VEC3F *cg);
void GetCofGPDBSCRange(PDB *start, PDB *stop, VEC3F *cg);
void OriginPDB(PDB *pdb);
void RotatePDB(PDB  *pdb, REAL rm[3][3]);
void TranslatePDB(PDB   *pdb, VEC3F tvect);
BOOL FitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitNCaCPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
BOOL FitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3]);
REAL CalcRMSPDB(PDB *pdb1, PDB *pdb2);
int GetPDBCoor(PDB *pdb, COOR **coor);
BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, 
                   char stopinsert, char chain, int mode, 
                   PDB **pdb_start, PDB **pdb_stop);
int HAddPDB(FILE *fp, PDB *pdb);
int ReadPGP(FILE *fp);
FILE *OpenPGPFile(char *pgpfile, BOOL AllHyd);
PDB *SelectAtomsPDB(PDB *pdbin, int nsel, char **sel, int *natom);
PDB *StripHPDB(PDB *pdbin, int *natom);
SECSTRUC *ReadSecPDB(FILE *fp, int *nsec);
void RenumAtomsPDB(PDB *pdb);
PDB *FindEndPDB(PDB *start);
PDB *FixOrderPDB(PDB *pdb, BOOL Pad, BOOL Renum);
PDB *ShuffleResPDB(PDB *start, PDB *end, BOOL Pad);
BOOL GetAtomTypes(char *resnam, char **AtomTypes);
PDB *KillPDB(PDB *pdb, PDB *prev);
void CopyPDB(PDB *out, PDB *in);
BOOL MovePDB(PDB *move, PDB **from, PDB **to);
PDB *AppendPDB(PDB *first, PDB *second);
PDB *ShuffleBB(PDB *pdb);
REAL CalcChi(PDB *pdb, int type);
PDB *GetPDBByN(PDB *pdb, int n);
void SetChi(PDB *pdb, PDB *next, REAL chi, int type);
BOOL KillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB);
void SetResnam(PDB *ResStart, PDB *NextRes, char *resnam, int resnum,   
                 char *insert, char *chain);
void ApplyMatrixPDB(PDB *pdb, REAL matrix[3][3]);
BOOL GetResolPDB(FILE *fp, REAL *resolution, REAL *RFactor, 
                 int *StrucType);
BOOL GetExptl(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
BOOL GetExptlOld(FILE *fp, REAL *resolution, REAL *RFactor, REAL *FreeR,
              int *StrucType);
char *ReportStructureType(int type);
PDB **IndexPDB(PDB *pdb, int *natom);
DISULPHIDE *ReadDisulphidesPDB(FILE *fp, BOOL *error);
BOOL ParseResSpec(char *spec, char *chain, int *resnum, char *insert);
BOOL ParseResSpecNoUpper(char *spec, char *chain, int *resnum, char *insert);
BOOL DoParseResSpec(char *spec, char *chain, int *resnum, char *insert, 
                      BOOL uppercaseresspec);
BOOL RepSChain(PDB *pdb, char *sequence, char *ChiTable, char *RefCoords);
PDB *FindNextChainPDB(PDB *pdb);
BOOL FixCterPDB(PDB *pdb, int style);
BOOL CalcCterCoords(PDB *p, PDB *ca_p, PDB *c_p, PDB *o_p);
int CalcTetraHCoords(PDB *nter, COOR *coor);
int AddNTerHs(PDB **ppdb, BOOL Charmm);
char *FNam2PDB(char *filename);
PDB *TermPDB(PDB *pdb, int length);
PDB *FindHetatmResidueSpec(PDB *pdb, char *resspec);
PDB *FindResidueSpec(PDB *pdb, char *resspec);
PDB *FindNextResidue(PDB *pdb);
PDB *DupePDB(PDB *in);
BOOL CopyPDBCoords(PDB *out, PDB *in);
void CalcCellTrans(VEC3F UnitCell, VEC3F CellAngles, 
                     VEC3F *xtrans, VEC3F *ytrans, VEC3F *ztrans);
int GetCrystPDB(FILE *fp, VEC3F *UnitCell, VEC3F *CellAngles,
                  char *spacegroup,
                  REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
void WriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                     char *spacegroup,
                     REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4]);
PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                      char *chain2, int resnum2, char *insert2);

PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert);
PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert);
PDB *FindAtomInRes(PDB *pdb, char *atnam);
BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1,
                 int resnum2, char insert2);
BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2);
BOOL AtomNameMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL AtomNameRawMatch(char *atnam, char *spec, BOOL *ErrorWarn);
BOOL LegalAtomSpec(char *spec);
BOOL RepOneSChain(PDB *pdb, char *ResSpec, char aa, char *ChiTable,
                    char *RefCoords);
void EndRepSChain(void);
char **ReadSeqresPDB(FILE *fp, int *nchains);
PDB *SelectCaPDB(PDB *pdb);
char *FixAtomName(char *name, REAL occup);

void FreeWholePDB(WHOLEPDB *wpdb);
WHOLEPDB *ReadWholePDB(FILE *fpin);
WHOLEPDB *ReadWholePDBAtoms(FILE *fpin);
BOOL AddCBtoGly(PDB *pdb);
BOOL AddCBtoAllGly(PDB *pdb);
PDB *StripGlyCB(PDB *pdb);
PDB *RemoveAlternates(PDB *pdb);
PDB *BuildAtomNeighbourPDBList(PDB *pdb, PDB *pRes, REAL NeighbDist);
PDB *FindAtomWildcardInRes(PDB *pdb, char *pattern);
PDB *DupeResiduePDB(PDB *in);
PDB *StripWatersPDB(PDB *pdbin, int *natom);
PDBSTRUCT *AllocPDBStructure(PDB *pdb);
PDB *FindNextChain(PDB *pdb);
void FreePDBStructure(PDBSTRUCT *pdbstruct);
void SetElementSymbolFromAtomName(char *element, char * atom_name);

#endif

/************************************************************************/
/* Deprecated functions: BuffInp.h                                      */

#ifdef _BUFFINPUT_H_DEPRECATED
#undef _BUFFINPUT_H_DEPRECATED

INBUFFER *OpenBufferedFile(char *filename, int maxstr);
BOOL ReadBufferedFile(INBUFFER *bfp, char *string, int length);
BOOL ProbeBufferedFile(INBUFFER *bfp, char *string, int length);

#endif

/************************************************************************/
/* Deprecated functions: ErrStack.h                                     */
#ifdef _ERRSTACK_H_DEPRECATED
#undef _ERRSTACK_H_DEPRECATED


void StoreError(char *routine, char *error);
void ShowErrors(void *PrintRoutine(char *), BOOL Trace);

#endif

/************************************************************************/
/* Deprecated functions: MathUtil.h                                     */

#ifdef _MATHUTIL_H_DEPRECATED
#undef _MATHUTIL_H_DEPRECATED


void CalcSD(REAL val, int action, REAL *mean, REAL *SD);
void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
               int *NValues, REAL *mean, REAL *SD);
REAL pearson(REAL *x, REAL *y, int NItem);
REAL pearson1(REAL *x, REAL *y, int NItem);

void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2);
void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2);
void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2);
REAL VecLen3(VEC3F Vec);
REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2);
REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac);
ULONG factorial(int n);
ULONG Factdiv(int n1, int n2);
ULONG NPerm(int n, int r);
ULONG NComb(int n, int r);

#endif

/************************************************************************/
/* Deprecated functions: WindIO.h                                       */

#ifdef _WINDIO_H_DEPRECATED
#undef _WINDIO_H_DEPRECATED


void screen(char *string);
void prompt(char *string);
void RePrompt(void);
void GetKybdString(char *string, int maxlen);
void PagingOn(void);
void PagingOff(void);
void WindowMode(BOOL mode);
void WindowInteractive(BOOL mode);
int YorN(char deflt);

#endif

/************************************************************************/
/* Deprecated functions: aalist.h                                       */

#ifdef _AALIST_H_DEPRECATED
#undef _AALIST_H_DEPRECATED


AA *InsertNextResiduesInAAList(AA *a, char res, int nres);
AA *InsertNextResidueInAAList(AA *a, char res);
char *BuildSeqFromAAList(AA *aa);
AA *InsertResidueInAAListAt(AA *aa, char res, int pos);
AA *InsertResiduesInAAListAt(AA *aa, char res, int nres, int pos);
AA *BuildAAList(char *seq);
int FindAAListOffsetByResnum(AA *aa, int resnum);
AA *FindAAListItemByResnum(AA *aa, int resnum);
void SetAAListFlagByResnum(AA *aa, int resnum);
char *BuildFlagSeqFromAAList(AA *aa, char ch);
int GetAAListLen(AA *aa);

#endif

/************************************************************************/
/* Deprecated functions: angle.h                                        */

#ifdef _ANGLE_H_DEPRECATED
#undef _ANGLE_H_DEPRECATED


REAL angle(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
           REAL xk, REAL yk, REAL zk);
REAL phi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj,
         REAL xk, REAL yk, REAL zk, REAL xl, REAL yl, REAL zl);
REAL simpleangle(REAL ang);
REAL TrueAngle(REAL opp, REAL adj);
BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, 
                 REAL bond, REAL theta, REAL torsion,
                 VEC3F *coords);
                 
#endif

/************************************************************************/
/* Deprecated functions: array.h                                        */

#ifdef _ARRAY_H_DEPRECATED
#undef _ARRAY_H_DEPRECATED


char **Array2D(int size, int dim1, int dim2);
void FreeArray2D(char **array, int dim1, int dim2);
char ***Array3D(int size, int dim1, int dim2, int dim3);
void FreeArray3D(char ***array, int dim1, int dim2, int dim3);

#endif

/************************************************************************/
/* Deprecated functions: cssr.h                                         */


#ifdef _CSSR_H_DEPRECATED
#undef _CSSR_H_DEPRECATED


CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title);
PDB *ReadCSSRasPDB(FILE *fp, int *natom);
void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta,
                   REAL gamma);
void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta,
                  REAL gamma);
void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
           REAL amatrx[3][3], int isw, int ncode);
void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title);

#endif

/************************************************************************/
/* Deprecated functions: fit.h                                          */

#ifdef _FIT_H_DEPRECATED
#undef _FIT_H_DEPRECATED


BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
              BOOL column);

#endif

/************************************************************************/
/* Deprecated functions: hbond.h                                        */

#ifdef _hbond_h_deprecated
#undef _hbond_h_deprecated


int  IsHBonded(PDB *res1, PDB *res2, int type);
BOOL ValidHBond(PDB *AtomH, PDB *AtomD, PDB *AtomA, PDB *AtomP);
int IsMCDonorHBonded(PDB *res1, PDB *res2, int type);
int IsMCAcceptorHBonded(PDB *res1, PDB *res2, int type);

#endif

/************************************************************************/
/* Deprecated functions: help.h                                         */

#ifdef _HELP_H_DEPRECATED
#undef _HELP_H_DEPRECATED


void Help(char *string, char *HelpFile);
void DoHelp(char *string, char *HelpFile);

#endif

/************************************************************************/
/* Deprecated functions: hpgl.h                                         */

#ifdef _HPGL_H_DEPRECATED
#undef _HPGL_H_DEPRECATED


BOOL HPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin);
void HPGLPen(int num);
void HPGLMove(REAL x, REAL y);
void HPGLDraw(REAL x, REAL y);
void HPGLSetDash(int style);
void HPGLFont(int font, REAL size);
void HPGLLText(REAL x, REAL y, char *string);
void HPGLCBText(REAL x, REAL y, REAL offset, char *text);
void HPGLROffText(REAL x, REAL y, REAL offset, char *text);
void HPGLLCText(REAL x, REAL y, char *text);
void HPGLCTText(REAL x, REAL y, REAL offset, char *text);
void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, 
               REAL TitleSize, char *label, int LabelFont, REAL LabelSize);
void HPGLEnd(void);
void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase);

#endif

/************************************************************************/
/* Deprecated functions: matrix.h                                       */

#ifdef _MATRIX_H_DEPRECATED
#undef _MATRIX_H_DEPRECATED


void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout);
void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3]);
void invert33(REAL s[3][3], REAL ss[3][3]);
void CreateRotMat(char direction, REAL angle, REAL matrix[3][3]);
REAL VecDist(REAL *a, REAL *b, int len);

#endif

/************************************************************************/
/* Deprecated functions: parse.h                                        */

#ifdef _PARSE_H_DEPRECATED
#undef _PARSE_H_DEPRECATED


int parse(char *comline, int nkeys, KeyWd *keywords, REAL *REALparam,
          char **strparam);
int mparse(char *comline, int nkeys, MKeyWd *keywords, REAL *REALparam,
           char **strparam, int *nparams);
int match(char *comstring, char *string2, int *nletters);
int GetString(char *command, char *strparam);
int GetParam(char  *command, REAL *value, int *nletters);

#endif

/************************************************************************/
/* Deprecated functions: plotting.h                                     */

#ifdef _PLOTTING_H_DEPRECATED
#undef _PLOTTING_H_DEPRECATED


BOOL AMInitPlot(char *filename, char *title, int dest, REAL OutXSize, 
                REAL OutYSize, REAL OutXOff, REAL OutYOff,
                char *AltFont, REAL xmargin, REAL ymargin,
                REAL DataXMin, REAL DataYMin, REAL DataXMax,
                REAL DataYMax);
void AMSetPen(int dest, int pen);
void AMMove(int dest, REAL x, REAL y);
void AMDraw(int dest, REAL x, REAL y);
void AMSetLineStyle(int dest, int style);
void AMEndLine(int dest);
void AMSetFont(int dest, char *PSFontName, REAL FontSize);
void AMText(int dest, REAL x, REAL y, char *text);
void AMCBText(int dest, REAL x, REAL y, char *text);
void AMRText(int dest, REAL x, REAL y, REAL offset, char *text);
void AMLCText(int dest, REAL x, REAL y, char *text);
void AMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text);
void AMEndPlot(int dest);
int  PS2HPGLFont(char *font);
char *SimplifyText(char *string);

#endif

/************************************************************************/
/* Deprecated functions: ps.h                                           */

#ifdef _PS_H_DEPRECATED
#undef _PS_H_DEPRECATED


BOOL PSInit(char *FName, char *creator, char *AltFont);
void PSThick(REAL thickness);
void PSMove(REAL X, REAL Y);
void PSDraw(REAL X, REAL Y);
void PSSetDash(char *linepatt);
void PSClearDash(void);
void PSStroke(void);
void PSFont(char *fontname, REAL size);
void PSLText(REAL X, REAL Y, char *label);
void PSCBText(REAL X, REAL Y, REAL Offset, char *label);
void PSROffText(REAL X, REAL Y, REAL offset, char *label);
void PSLCText(REAL X, REAL Y, char *label);
void PSCTText(REAL X, REAL Y, REAL Offset, char *label);
void PSVText(REAL x, REAL y, REAL xoff, char *text, char *font, REAL size,
             char *label, char *lfont, REAL lsize);
void PSShowText(char *text);
void PSEnd(void);
char *PSCorrectCase(char *font);

#endif

/************************************************************************/
/* Deprecated functions: safemem.h                                      */

#ifdef _SAFEMEM_H_DEPRECATED
#undef _SAFEMEM_H_DEPRECATED


void *safemalloc(int nbytes);
BOOL safefree(void *ptr);
void safeleaks(void);

#endif

/************************************************************************/
/* Deprecated macros: seq.h                                             */

#ifdef _SEQ_H_DEPRECATED
#undef _SEQ_H_DEPRECATED

/*    ***CRAIG*** These won't give deprecation messages! ***CRAIG***    */
#define PDB2Seq(x)          blDoPDB2Seq((x), FALSE, FALSE, FALSE)
#define PDB2SeqX(x)         blDoPDB2Seq((x), TRUE,  FALSE, FALSE)
#define PDB2SeqNoX(x)       blDoPDB2Seq((x), FALSE, FALSE, TRUE)
#define PDB2SeqXNoX(x)      blDoPDB2Seq((x), TRUE,  FALSE, TRUE)

#define PDBProt2Seq(x)      blDoPDB2Seq((x), FALSE, TRUE, FALSE)
#define PDBProt2SeqX(x)     blDoPDB2Seq((x), TRUE,  TRUE, FALSE)
#define PDBProt2SeqNoX(x)   blDoPDB2Seq((x), FALSE, TRUE, TRUE)
#define PDBProt2SeqXNoX(x)  blDoPDB2Seq((x), TRUE,  TRUE, TRUE)

/************************************************************************/
/* Deprecated functions: seq.h                                          */

char throne(char *three);
char thronex(char *three);
char *onethr(char one);
char *DoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX);
int SplitSeq(char *LinearSeq, char **seqs);
int ReadSimplePIR(FILE *fp, int  maxres, char **seqs);
int ReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
            SEQINFO *seqinfo, BOOL *punct, BOOL *error);
int ReadRawPIR(FILE *fp, char **seqs, int maxchain, BOOL upcase,
               SEQINFO *seqinfo, BOOL *error);
int align(char *seq1, int  length1, char *seq2, int  length2, 
          BOOL verbose, BOOL identity, int  penalty, 
          char *align1, char *align2, int  *align_len);
int affinealign(char *seq1, int  length1, char *seq2, int  length2, 
                BOOL verbose, BOOL identity, int  penalty, int penext,
                char *align1, char *align2, int  *align_len);
int CalcMDMScore(char resa, char resb);
int affinealignuc(char *seq1, int  length1, char *seq2, int  length2, 
                  BOOL verbose, BOOL identity, int  penalty, int penext,
                  char *align1, char *align2, int  *align_len);
int CalcMDMScoreUC(char resa, char resb);
BOOL ReadMDM(char *mdmfile);
int ZeroMDM(void);
char DNAtoAA(char *dna);
int TrueSeqLen(char *sequence);
int KnownSeqLen(char *sequence);
BOOL NumericReadMDM(char *mdmfile);
int NumericCalcMDMScore(int resa, int resb);
int NumericAffineAlign(int *seq1, int length1, int *seq2, int length2, 
                       BOOL verbose, BOOL identity, int penalty,
                       int penext, int *align1, int *align2, 
                       int *align_len);
void WeightMDMScore(char resa, char resb, REAL weight);

#endif

/************************************************************************/
/** \endcond                                                            */
/************************************************************************/

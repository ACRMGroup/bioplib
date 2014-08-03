/************************************************************************/
/**

   \file       seq.h
   
   \version    V2.12
   \date       31.07.14
   \brief      Header file for sequence handling
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1991-2014
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
-  V2.0  11.03.94 Original V2 release
-  V2.1  11.05.94 Added DNAtoAA() & TrueSeqLen() prototypes
-  V2.2  13.05.93 Added KnownSeqLen() prototype
-  V2.3  28.02.95 Added ReadRawPIR()
-  V2.4  25.07.95 Added the gBioplibSeqNucleicAcid external for throne()
-  V2.5  11.07.96 Added CalcMDMScore()
-  V2.6  17.09.95 Added ZeroMDM()
-  V2.7  26.08.97 Added macro interfaces to new DoPDB2Seq()
-  V2.8  08.03.00 Added Numeric***() alignment routines
-  V2.9  02.10.00 Modified DoPDB2Seq()
-  V2.10 27.02.07 Added CalcMDMScoreUC() and affinealignuc()
-  V2.11 07.07.14 Use bl prefix for functions. Renamed PDB2Seq macros to
                  blPDB2Seq and added PDB2Seq defines to deprecated.h 
                  By: CTP
-  V2.12 31.07.14 Updated deprecation: Removed deprecated.h, added 
                  prototypes for renamed functions and defines for PDB2Seq
                  macros. By: CTP

*************************************************************************/
#ifndef _SEQ_H
#define _SEQ_H

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
 */
#define ALLOCSIZE 80  /* ReadPIR() uses this as a chunk size for 
                         allocating memory
                      */

typedef struct
{
   BOOL fragment,
        paren,
        DotInParen,
        NonExpJoin,
        UnknownPos,
        Incomplete,
        Truncated,
        Juxtapose;
   char code[16],
        name[160],
        source[160];
}  SEQINFO;

extern BOOL gBioplibSeqNucleicAcid;

#define blPDB2Seq(x)         blDoPDB2Seq((x), FALSE, FALSE, FALSE)
#define blPDB2SeqX(x)        blDoPDB2Seq((x), TRUE,  FALSE, FALSE)
#define blPDB2SeqNoX(x)      blDoPDB2Seq((x), FALSE, FALSE, TRUE)
#define blPDB2SeqXNoX(x)     blDoPDB2Seq((x), TRUE,  FALSE, TRUE)

#define blPDBProt2Seq(x)     blDoPDB2Seq((x), FALSE, TRUE, FALSE)
#define blPDBProt2SeqX(x)    blDoPDB2Seq((x), TRUE,  TRUE, FALSE)
#define blPDBProt2SeqNoX(x)  blDoPDB2Seq((x), FALSE, TRUE, TRUE)
#define blPDBProt2SeqXNoX(x) blDoPDB2Seq((x), TRUE,  TRUE, TRUE)

char blThrone(char *three);
char blThronex(char *three);
char *blOnethr(char one);
char *blDoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX);
int blSplitSeq(char *LinearSeq, char **seqs);
int blReadSimplePIR(FILE *fp, int  maxres, char **seqs);
int blReadPIR(FILE *fp, BOOL DoInsert, char **seqs, int maxchain, 
              SEQINFO *seqinfo, BOOL *punct, BOOL *error);
int blReadRawPIR(FILE *fp, char **seqs, int maxchain, BOOL upcase,
                 SEQINFO *seqinfo, BOOL *error);
int blAlign(char *seq1, int  length1, char *seq2, int  length2, 
            BOOL verbose, BOOL identity, int  penalty, 
            char *align1, char *align2, int  *align_len);
int blAffinealign(char *seq1, int  length1, char *seq2, int  length2, 
                  BOOL verbose, BOOL identity, int  penalty, int penext,
                  char *align1, char *align2, int  *align_len);
int blCalcMDMScore(char resa, char resb);
int blAffinealignuc(char *seq1, int  length1, char *seq2, int  length2, 
                    BOOL verbose, BOOL identity, int  penalty, int penext,
                    char *align1, char *align2, int  *align_len);
int blCalcMDMScoreUC(char resa, char resb);
BOOL blReadMDM(char *mdmfile);
int blZeroMDM(void);
char blDNAtoAA(char *dna);
int blTrueSeqLen(char *sequence);
int blKnownSeqLen(char *sequence);
BOOL blNumericReadMDM(char *mdmfile);
int blNumericCalcMDMScore(int resa, int resb);
int blNumericAffineAlign(int *seq1, int length1, int *seq2, int length2, 
                         BOOL verbose, BOOL identity, int penalty,
                         int penext, int *align1, int *align2, 
                         int *align_len);



/************************************************************************/
/* Deprecated macros: seq.h                                             */

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
/************************************************************************/

#endif

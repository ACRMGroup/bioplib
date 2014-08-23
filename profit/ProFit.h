/*************************************************************************

   Program:    ProFit
   File:       ProFit.h
   
   Version:    V3.1
   Date:       31.03.09
   Function:   Protein Fitting program. Includes and defines.
   
   Copyright:  SciTech Software / UCL 1992-2009
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain.

   It may not be copied or made available to third parties, but may be
   freely used by non-profit-making organisations who have obtained it
   directly from the author or by FTP.

   You are requested to send EMail to the author to say that you are 
   using this code so that you may be informed of future updates.

   The code may not be made available on other FTP sites without express
   permission from the author.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If
   someone else breaks this code, the author doesn't want to be blamed
   for code that does not work! You may not distribute any
   modifications, but are encouraged to send them to the author so
   that they may be incorporated into future versions of the code.

   Such modifications become the property of Dr. Andrew C.R. Martin and
   SciTech Software though their origin will be acknowledged.

   The code may not be sold commercially or used for commercial purposes
   without prior permission from the author.
   
**************************************************************************

   Description:
   ============
   ProFit is a least squares fitting program for proteins.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  25.09.92 Original
   V0.5  08.10.93 Various tidying for Unix & book routines
   V0.6  05.01.94 Modified help and data defines for unix getenv()
   V0.7  24.11.94 Skipped
   V0.8  17.07.95 Removed all windowing stuff
   V1.0  18.07.95 First official release (at last!).
   V1.1  20.07.95 Added WEIGHT command support and translation vector
                  output from MATRIX command
   V1.2  22.07.95 Added GAPPEN command
   V1.3  31.07.95 Skipped
   V1.4  14.08.95 Skipped
   V1.5  21.08.95 Skipped
   V1.6  20.11.95 Added READALIGNMENT command
   V1.6e 31.05.96 Added BVALUE command
   V1.6f 13.06.96 Added BWEIGHT command; changed gDoWeights initialisation
   V1.6g 18.06.96 Removed MODE_SEQUENTIAL and MODE_RESNUM. Replaced
                  by ZONE_* versions from pdb.h
   V1.7  23.07.96 Added MAXATSPEC
   V1.7b 11.11.96 gUseBVal is now handled as an int
   V1.7c 18.11.96 Added IGNOREMISSING option
   V1.7d 20.12.96 Added gNFittedCoor and NFITTED command
   V1.8  07.05.98 Skipped for release
   V2.0  01.03.01 More commands; various things made arrays for multiple
                  structure fitting
   V2.1  28.03.01 Parameter for ITERATE and added CENTRE command
   V2.2  20.12.01 Skipped for release
   V2.3  01.12.04 Skipped for release
   V2.4  03.06.05 Skipped for release
   V2.5  07.06.05 Skipped for release
   V-.-  28.03.08 Added gCZoneList[]. Increased NCOMM for new commands. 
   V-.-  02.04.08 Added globals for handling PDB headers and footers.
   V-.-  02.04.08 Increased NCOMM for new command. 
   V-.-  07.04.08 Added distance cutoff for including atom pairs in RMSd.
   V-.-  15.04.08 Added WHOLEPDB linked lists.
   V-.-  01.05.08 Added gOccRank to control reading of low-occupancy atoms.
   V-.-  02.05.08 Headers and footers now handled by WHOLEPDB.
   V2.6  28.05.08 Removed unused globals and defines for headers and 
                  footers.
   V-.-  04.06.08 Added gMatchSymAtoms and gSymType[] for matching symmetrical
                  atom pairs(eg CD1 - CD2 and CE1 - CE2 in Tyr)
                  Increased NCOMM for new command to set symmetric matching.
   V-.-  16.07.08 Added gap extension parameter, gGapPenExt, used in alignment
                  functions. Set default parameters for alignment to:
                  gGapPen = 10 and gGapPenExt to 2.
   V-.-  18.07.08 Increased NCOMM for new command.
   V-.-  30.07.08 Included bioplib/matrix.h.
   V2.6  07.08.08 Added gMultiVsRef.
   V2.6  20.08.08 Added gTwistAngle and gTwistMatrix used for rotating and 
                  refitting a mobile structure.
   V2.6  21.08.08 Added gRotateRefit flag to switch on rotate and refit.
   V2.6  23.10.08 Added gWtAverage flag to allow old use of old weighting 
                  system.
   V3.0  06.11.08 Release Version.
   V3.0  07.11.08 Added reference number,gMultiRef, for multistructure.
   V3.0  14.11.08 Added GNU Readline Library support.
   V3.0  04.02.09 Replaced gRotateRefit with ROTATE_REFIT #define
   V3.0  18.02.09 Added #include for bioplib/aalist.h
   V3.0  04.02.09 Removed ROTATE_REFIT #define - Now compile option.
   V3.1  31.03.09 Skipped for release

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "bioplib/MathType.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/parse.h"
#include "bioplib/fit.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"
#include "bioplib/help.h"
#include "bioplib/array.h"
#include "bioplib/matrix.h"
#include "bioplib/aalist.h"

#ifdef READLINE_SUPPORT
#undef NEWLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

/************************************************************************/
/* Defines
*/
#define NUMTYPES       50     /* Number of atom types                   */
#define MAXATSPEC       8     /* Max length of a single atom spec       */
#define NCOMM          53     /* Number of keywords                     */
#define MAXNUMPARAM     2     /* Max number of numeric parameters       */
#define MAXSTRPARAM     2     /* Max number of string parameters        */
#define MAXSTRLEN     160     /* Max length of returned string          */
#define MAXCHAIN      160     /* Max number of chains in an align file  */
#define MAXBUFF       256     /* Buffer for filenames etc               */
#define MAXSTRUC     1000     /* Max number of structures to fit        */
#define MAXITER      1000     /* Max allowed number of iterations for
                                 updating zones                         */
#define MAXMULTIITER  100     /* Max allowed number of iterations for
                                 multiple structures                    */
#define DEF_MAXEQUIVDISTSQ 9.0/* Max distance before a pair is added to
                                 equivalence list in iterative mode     */
#define ITER_STOP    0.01     /* Change in RMSD for zone convergence    */
#define MULTI_ITER_STOP 0.001 /* Change in RMSD for multi struc 
                                 convergence                            */
#define STRUC_REFERENCE 0     /* Flags for which structure to be loaded */
#define STRUC_MOBILE    1

#define ATOM_FITTING    0     /* Flags for ValidAtom()                  */
#define ATOM_RMS        1

#define WEIGHT_NONE     0     /* B-value weighting schemes              */
#define WEIGHT_BVAL     1
#define WEIGHT_INVBVAL  2

#define SYMM_ATM_PAIRS 11     /* Number of symmetric atom pairs         */

typedef FILE file;

#define PDBDISTSQ(a, b) ((((a)->x - (b)->x) * ((a)->x - (b)->x)) + \
                         (((a)->y - (b)->y) * ((a)->y - (b)->y)) + \
                         (((a)->z - (b)->z) * ((a)->z - (b)->z)))

#define HELPFILE "ProFit.help"   /* Help file                           */
#define MDMFILE  "mdm78.mat"     /* Mutation data matrix                */

/************************************************************************/
/* Structure definitions
*/
typedef struct zonestruct
{
   struct zonestruct *next;
   int               start1,
                     stop1,
                     start2,
                     stop2,
                     mode;
   char              chain1,
                     chain2,
                     startinsert1,
                     startinsert2,
                     stopinsert1,
                     stopinsert2;
}  ZONE;

/* Type definition to store a X,Y coordinate pair in the matrix         */
typedef struct
{
   
   int x, y;
   
}
XY;

/************************************************************************/
/* Prototype definitions
*/
#include "protos.h"

/************************************************************************/
/* Globals
*/
#ifdef MAIN /*----------------------------------------------------------*/
char   gFitAtoms[NUMTYPES][MAXATSPEC], /* Atom types to be fitted       */
       gRMSAtoms[NUMTYPES][MAXATSPEC], /* Atom types for RMS calculation*/
       gRefFilename[MAXBUFF],          /* Reference filename            */
       gMobFilename[MAXSTRUC][MAXBUFF],/* Mobile filename               */
       *gRefSeq = NULL,                /* Sequences                     */
       *gMobSeq[MAXSTRUC],
       gSymType[SYMM_ATM_PAIRS][4][5]; /* Symmetric Atom Types          */
WHOLEPDB  *gRefWPDB = NULL,            /* WHOLEPDB linked lists         */
          *gMobWPDB[MAXSTRUC],
          *gFitWPDB[MAXSTRUC];
PDB       *gRefPDB = NULL,             /* PDB linked lists              */
          *gMobPDB[MAXSTRUC],
          *gFitPDB[MAXSTRUC];
MKeyWd gKeyWords[NCOMM];            /* Array to store keywords          */
char   *gStrParam[MAXSTRPARAM];     /* Array for returned strings       */
REAL   gNumParam[MAXNUMPARAM],      /* Array for returned numbers       */
       gRotMat[MAXSTRUC][3][3],     /* Rotation matrix                  */
       *gWeights      = NULL,       /* Weights array                    */
       gBValue        = 10000.0,    /* Max BVal to consider in fitting  */
       gMaxEquivDistSq= DEF_MAXEQUIVDISTSQ, /* Max distance before a pair
                                       is added to equivalence list in 
                                       iterative mode                   */
       gDistCutoff    = 0.0,        /* Distance cutoff for including atom
                                       pairs when calculating RMSd */
       gTwistAngle    = 42.0,       /* Rotation angle for refit routine */
       gRotMatTwist[3][3];          /* Rotation matrix for refit routine*/

ZONE   *gZoneList[MAXSTRUC],        /* List of zones                    */
       *gRZoneList[MAXSTRUC],
       *gCZoneList[MAXSTRUC];
int    gCurrentMode   = ZONE_MODE_RESNUM, /* Numbering mode             */
       gUserRMSZone   = FALSE,      /* User has specified things for RMS*/
       gUserRMSAtoms  = FALSE,
       gUserFitZone   = FALSE,      /* User has specified fit zone      */
       gFitted        = FALSE,      /* Structures fitted                */
       gNOTFitAtoms   = FALSE,      /* NOT atom selections              */
       gNOTRMSAtoms   = FALSE,
       gHetAtoms      = FALSE,      /* Include het atoms?               */
       gIterate       = FALSE,      /* Iterative fitting?               */
       gDoWeights     = WEIGHT_NONE,/* Weight by BVal column?           */
       gGapPen        = 10,         /* Align Gap penalty                */
       gGapPenExt     = 2,          /* Align Gap Extension penalty      */
       gUseBVal       = 0,          /* Use BValue cutoff on atom sel    */
       gIgnoreMissing = 0,          /* Ignore missing atoms             */
       gNFittedCoor   = 0,          /* Number of coordinates fitted     */
       gMultiCount    = 0,          /* Number of strucs in multi mode   */
       gQuiet         = FALSE,      /* Stop warning messages            */
       gCentre        = FALSE,      /* Leave structure centred at origin*/
       gLimit[2],                   /* Limit range from alignment       */
       gReadHeader    = FALSE,      /* Read/Write PDB Headers           */
       gUseDistCutoff = FALSE,      /* Use distance cutoff for including
                                       atom pairs when calculating RMSd */
       gOccRank       = 1,          /* Occupancy ranking (>= 1)         */
       gMatchSymAtoms = FALSE,      /* Match Symmetrical Atoms          */
       gMultiVsRef    = FALSE,      /* Set multi RMSD calculations      */
       gWtAverage     = TRUE,       /* Weighted averaging in multi mode */
       gMultiRef      = 0;          /* Mobile used as multi reference   */
COOR   *gRefCoor      = NULL,       /* Coordinate arrays                */
       *gMobCoor[MAXSTRUC];
VEC3F  gRefCofG,                    /* CofG of fitted region            */
       gMobCofG[MAXSTRUC];
#else       /*----------------------------------------------------------*/
extern char    gFitAtoms[NUMTYPES][MAXATSPEC],
               gRMSAtoms[NUMTYPES][MAXATSPEC],
               gRefFilename[MAXBUFF],
               gMobFilename[MAXSTRUC][MAXBUFF],
               *gRefSeq,
               *gMobSeq[MAXSTRUC],
               gSymType[SYMM_ATM_PAIRS][4][5];
extern WHOLEPDB  *gRefWPDB,
                 *gMobWPDB[MAXSTRUC],
                 *gFitWPDB[MAXSTRUC];
extern PDB     *gRefPDB,
               *gMobPDB[MAXSTRUC],
               *gFitPDB[MAXSTRUC];
extern MKeyWd  gKeyWords[NCOMM];
extern char    *gStrParam[MAXSTRPARAM];
extern REAL    gNumParam[MAXNUMPARAM],
               gRotMat[MAXSTRUC][3][3],
               *gWeights,
               gBValue,
               gMaxEquivDistSq,
               gDistCutoff,
               gTwistAngle,
               gRotMatTwist[3][3];
extern ZONE    *gZoneList[MAXSTRUC],
               *gRZoneList[MAXSTRUC],
               *gCZoneList[MAXSTRUC];
extern int     gCurrentMode,
               gUserRMSZone,
               gUserRMSAtoms,
               gUserFitZone,
               gFitted,
               gNOTFitAtoms,
               gNOTRMSAtoms,
               gHetAtoms,
               gIterate,
               gDoWeights,
               gGapPen,
               gGapPenExt,
               gUseBVal,
               gIgnoreMissing,
               gNFittedCoor,
               gMultiCount,
               gQuiet,
               gCentre,
               gLimit[2],
               gReadHeader,
               gUseDistCutoff,
               gOccRank,
               gMatchSymAtoms,
               gMultiVsRef,
               gWtAverage,
               gMultiRef;

extern COOR    *gRefCoor,
               *gMobCoor[MAXSTRUC];
extern VEC3F   gRefCofG,
               gMobCofG[MAXSTRUC];
#endif      /*----------------------------------------------------------*/



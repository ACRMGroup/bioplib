/************************************************************************/
/**

   \file       general.h
   
   \version    V1.22
   \date       10.11.17
   \brief      Header file for general purpose routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2017
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
-  V1.0  11.05.94 Original    By: ACRM
-  V1.1  24.08.94 Added OpenStdFiles()
-  V1.2  22.09.94 Added OpenFile()
-  V1.3  17.07.95 Added countchar()
-  V1.4  11.09.95 Added fgetsany()
-  V1.5  18.10.95 Moved YorN() to WindIO.h
-  V1.6  06.11.95 Added StoreString(), InStringList() and FreeStringList()
-  V1.7  15.12.95 Added QueryStrStr()
-  V1.8  09.07.96 Added IndxReal()
-  V1.9  18.09.96 Added padchar()
-  V1.11 13.06.00 Added strcatalloc()
-  V1.12 30.05.02 Added WrapString(), WrapPrint(), RightJustify(), 
                  GetWordNC() and getfield()
-  V1.13 07.07.14 Use bl prefix for functions By: CTP
-  V1.14 31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.15 14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP
-  V1.16 10.03.15 Added blSplitStringOnCommas()  By: ACRM
-  V1.17 12.03.15 Added blSplitStringOnChars(), blCheckProgName()
-  V1.18 26.03.15 Added blStrncat()
-  V1.19 28.04.15 Added blCollapseSpaces()
-  V1.20 14.05.15 Added blStrdup()
-  V1.21 26.06.15 Added FREESTRINGLIST() macro
-  V1.22 10.11.17 Added blRemoveSpaces()

*************************************************************************/
#ifndef _GENERAL_H
#define _GENERAL_H

#include <stdio.h>
#include "SysDefs.h"
#include "MathType.h"

typedef struct _stringlist
{
   struct _stringlist *next;
   char               *string;
}  STRINGLIST;

#define FREESTRINGLIST(l) do {                          \
      STRINGLIST *_s = NULL;                            \
      for(_s = (l); _s != NULL; NEXT(_s)) {             \
         if(_s->string != NULL) free(_s->string);       \
      }                                                 \
      FREELIST((l), STRINGLIST);                        \
   }                                                    \
   while(0)


void blStringToLower(char *string1, char *string2);
void blStringToUpper(char *string1, char *string2);
char *blKillLeadSpaces(char *string);
void blKillLine(FILE *fp);
void blSetExtn(char *File, char *Ext);
int blChindex(char *string, char ch);
void blWord(char *string1, char *string2);
void blWordN(char *string1, char *string2, int  MaxChar);
void blPadterm(char *string, int length); /* defined in cssr.h */
void blPadchar(char *string, int length, char ch);
BOOL blCheckExtn(char *string, char *ext);
char *blFtostr(char *str, int maxlen, REAL x, int precision);

void blGetFilestem(char *filename, char *stem);
int blUpstrcmp(char *word1, char *word2);
int blUpstrncmp(char *word1, char *word2, int ncomp);
char *blGetWord(char *buffer, char *word, int maxsize);
char **blSplitStringOnCommas(char *string, int minItemLen);
char **blSplitStringOnChars(char *string);
BOOL blOpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out);
FILE *blOpenFile(char *filename, char *envvar, char *mode, BOOL *noenv);
int blCountchar(char *string, char ch);
char *blFgetsany(FILE *fp);
char *blStrcatalloc(char *instr, char *catstr);
char *blStrncat(char *out, const char *in, size_t len);

STRINGLIST *blStoreString(STRINGLIST *StringList, char *string);
BOOL blInStringList(STRINGLIST *StringList, char *string);
void blFreeStringList(STRINGLIST *StringList);

char *blQueryStrStr(char *string, char *substring);

void blIndexReal(REAL *arrin, int *indx, int n);

FILE *blOpenOrPipe(char *filename);
int blCloseOrPipe(FILE *fp);

BOOL blWrapString(char *in, char *out, int maxlen);
BOOL blWrapPrint(FILE *out, char *string);
void blRightJustify(char *string);
char *blGetWordNC(char *buffer, char *word, int maxlen);
void blGetfield(char *buffer, int start, int width, char *str);
BOOL blCheckProgName(char *name, char *expected);
char *blCollapseSpaces(char *inText);
char *blRemoveSpaces(char *inText);

/************************************************************************/
/* Include deprecated functions                                         */
#define _GENERAL_H_DEPRECATED
#include "deprecated.h"
/************************************************************************/

#endif

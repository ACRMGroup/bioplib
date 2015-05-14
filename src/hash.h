/************************************************************************/
/**

   \file       hash.h
   
   \version    V1.0
   \date       14.05.15
   \brief      Defines for using hash functions
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2015
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
-  V1.0  14.05.15 Original    By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdarg.h>
#include "SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define DEF_HASHSIZE     1009    /* Must be a prime number              */
#define HASHTYPE_STRING  0
#define HASHTYPE_INT     1
#define HASHTYPE_CHAR    2
#define HASHTYPE_DOUBLE  3
#define HASHTYPE_POINTER 4

/* Internal hash structure for actual entries                           */
typedef struct _hashTab
{
   struct _hashTab *next;
   char            *key;
   char            *string;
   BPTR            ptrValue;
   double          doubleValue;
   int             intValue;
   int             type;
   char            charValue;
}  _HASHENTRY;

/* End user hash structure                                              */
typedef struct
{
   _HASHENTRY    **table;
   ULONG         size;
}  HASHTABLE;


/************************************************************************/
/* Prototypes
*/
HASHTABLE *blInitializeHash(ULONG hashsize);
void      blFreeHash(HASHTABLE *hashtable);
char      *blStrdup(char *instr);
void      blFreeHashKeyList(char **keylist);
char      **blGetHashKeyList(HASHTABLE *hashtable);
BOOL      blSetHashValue(HASHTABLE *hashtable, char *key, int type, ...);
BPTR      blGetHashValue(HASHTABLE *hashtable, char *key, int *type);
BOOL      blDumpHash(FILE *out, HASHTABLE *hashtable);
int       blGetHashValueInt(HASHTABLE *hashtable, char *key);
double    blGetHashValueDouble(HASHTABLE *hashtable, char *key);
char      blGetHashValueChar(HASHTABLE *hashtable, char *key);
char      *blGetHashValueString(HASHTABLE *hashtable, char *key);
BPTR      blGetHashValuePointer(HASHTABLE *hashtable, char *key);
BOOL      blHashKeyDefined(HASHTABLE *hashtable, char *key);
void      blDeleteHashKey(HASHTABLE *hashtable, char *key);

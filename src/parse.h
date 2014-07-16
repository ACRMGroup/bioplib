/************************************************************************/
/**

   \file       Parse.h
   
   \version    V1.9
   \date       07.07.14
   \brief      Include file for the command parser
   
   \copyright  SciTech Software 1991-2014
   \author     Andrew C. R. Martin
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

   Here are defined the MAKEKEY macro, STRING and NUMBER defines, the
   KeyWd structure and return values for the parser.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  11.07.90 Original
-  V1.1  08.12.92 Defines prototypes
-  V1.2  16.06.93 Added memory check to MAKEKEY
-  V1.3-1.6       Skipped
-  V1.7  01.03.94 Added mparse()
-  V1.8  11.03.94 Skipped
-  V1.9  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
#ifndef _PARSE_H
#define _PARSE_H

/************************************************************************/
/* Includes
*/
#include "MathType.h"
#include "deprecated.h"

/************************************************************************/
/* Defines
*/

#define STRING         1    /* Defines used for the MAKEKEY macro       */
#define NUMBER         0

#define PARSE_ERRC    -1    /* Error return values from parse()         */
#define PARSE_ERRP    -2
#define PARSE_COMMENT -3

/************************************************************************/
/* Type definitions
*/
typedef struct              /* Used to store keywords for parse()       */
{
   char  *name;
   int   string, nparam;
}  KeyWd;

typedef struct              /* Used to store keywords for mparse()      */
{
   char  *name;
   int   string, minparam, maxparam;
}  MKeyWd;

/************************************************************************/
/* Macros
*/
/* Create a keyword for parse()                                         */
#define MAKEKEY(x,w,v,z) \
        (x).name = (char *)malloc((strlen(w)+2) * sizeof(char)); \
        if((x).name != NULL) strcpy((x).name,w); \
                              (x).string = v; \
                              (x).nparam = z

/* Create a keyword for mparse()                                       */
#define MAKEMKEY(x,w,v,mn,mx) \
        (x).name = (char *)malloc((strlen(w)+2) * sizeof(char)); \
        if((x).name != NULL) strcpy((x).name,w); \
                              (x).string   = v; \
                              (x).minparam = mn; \
                              (x).maxparam = mx

/************************************************************************/
/* Prototypes
*/
int blParse(char *comline, int nkeys, KeyWd *keywords, REAL *REALparam,
          char **strparam);
int blMparse(char *comline, int nkeys, MKeyWd *keywords, REAL *REALparam,
          char **strparam, int *nparams);
int blMatch(char *comstring, char *string2, int *nletters);
int blGetString(char *command, char *strparam);
int blGetParam(char  *command, REAL *value, int *nletters);

#endif

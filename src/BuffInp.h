/************************************************************************/
/**

   \file       BuffInp.h
   
   \version    V1.3
   \date       14.08.14
   \brief      Header file for BuffInp.c
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2014
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

-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.3  14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP

*************************************************************************/
#ifndef _BUFFINPUT_H
#define _BUFFINPUT_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "SysDefs.h"

/************************************************************************/
/* Defines
*/
typedef struct
{
   FILE *fp;
   char *buffer;
   int  nlines,
        maxstr;
}  INBUFFER;

/************************************************************************/
/* Prototypes
*/
INBUFFER *blOpenBufferedFile(char *filename, int maxstr);
BOOL blReadBufferedFile(INBUFFER *bfp, char *string, int length);
BOOL blProbeBufferedFile(INBUFFER *bfp, char *string, int length);

/************************************************************************/
/* Include deprecated functions                                         */
#define _BUFFINPUT_H_DEPRECATED
#include "deprecated.h"
/************************************************************************/

#endif

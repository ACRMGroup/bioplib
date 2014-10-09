/************************************************************************/
/**

   \file       WindIO.h
   
   \version    V1.6
   \date       14.08.14
   \brief      Header for window/normal interface routines

   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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
-  V1.4  07.07.14 Use bl prefix for functions By: CTP
-  V1.5  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.6  14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP

*************************************************************************/
#ifndef _WINDIO_H
#define _WINDIO_H

#include "SysDefs.h"

/* Prototypes */
void blScreen(char *string);
void blPrompt(char *string);
void blRePrompt(void);
void blGetKybdString(char *string, int maxlen);
void blPagingOn(void);
void blPagingOff(void);
void blWindowMode(BOOL mode);
void blWindowInteractive(BOOL mode);
int blYorN(char deflt);

/************************************************************************/
/* Include deprecated functions                                         */
#define _WINDIO_H_DEPRECATED
#include "deprecated.h" 
/************************************************************************/


#endif

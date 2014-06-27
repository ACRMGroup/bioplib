/************************************************************************/
/**

   \file       ps.h
   
   \version    V1.11R
   \date       23.06.94
   \brief      Include file for PostScript routine
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-4
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

   Include file for using with PSRoutines.c
   Variables are defined only from the main program. Otherwise they
   are referenced as external.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  06.02.91 Original
-  V1.7  25.02.91 Fixed prototypes and definition of PSFile.
-  V1.10 07.05.92 Changed all prototypes to doubles
-  V1.11 23.06.94 Made gPSFile a global

*************************************************************************/
#ifndef _PS_H
#define _PS_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "SysDefs.h"
#include "MathType.h"

/************************************************************************/
/* Globals
*/
#ifdef _PS_MAIN   /* ------------------ _PS_MAIN ---------------------- */
   REAL   PSxpicsize    = 5.0,
          PSypicsize    = 5.0,
          PSxoffset     = 1.0,
          PSyoffset     = 2.0;
   FILE   *gPSFile      = NULL;
#else             /* ----------------- Not _PS_MAIN ------------------- */
   extern REAL    PSxpicsize,
                  PSypicsize,
                  PSxoffset,
                  PSyoffset;
   extern FILE    *gPSFile;
#endif            /* -------------------------------------------------- */

/************************************************************************/
/* Prototypes
*/
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

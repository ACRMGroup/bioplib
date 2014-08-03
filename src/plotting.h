/************************************************************************/
/**

   \file       plotting.h
   
   \version    V1.2
   \date       31.07.14
   \brief      Include file for using plotting routines
   
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
   
-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP

*************************************************************************/
#ifndef _PLOTTING_H
#define _PLOTTING_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "SysDefs.h"
#include "MathType.h"

#include "hpgl.h"
#include "ps.h"

/************************************************************************/
/* Defines
*/
#define DEST_SCREEN  0
#define DEST_PS      1
#define DEST_HPGL    2

/************************************************************************/
/* Prototypes
*/
BOOL blAMInitPlot(char *filename, char *title, int dest, REAL OutXSize, 
                  REAL OutYSize, REAL OutXOff, REAL OutYOff,
                  char *AltFont, REAL xmargin, REAL ymargin,
                  REAL DataXMin, REAL DataYMin, REAL DataXMax,
                  REAL DataYMax);
void blAMSetPen(int dest, int pen);
void blAMMove(int dest, REAL x, REAL y);
void blAMDraw(int dest, REAL x, REAL y);
void blAMSetLineStyle(int dest, int style);
void blAMEndLine(int dest);
void blAMSetFont(int dest, char *PSFontName, REAL FontSize);
void blAMText(int dest, REAL x, REAL y, char *text);
void blAMCBText(int dest, REAL x, REAL y, char *text);
void blAMRText(int dest, REAL x, REAL y, REAL offset, char *text);
void blAMLCText(int dest, REAL x, REAL y, char *text);
void blAMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text);
void blAMEndPlot(int dest);
int  blPS2HPGLFont(char *font);
char *blSimplifyText(char *string);

/************************************************************************/
/* Deprecated functions: plotting.h                                     */
/** \cond deprecated                                                    */

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

/* \endcond                                                             */
/************************************************************************/

#endif

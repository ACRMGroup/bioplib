/************************************************************************/
/**

   \file       plotting.h
   
   \version    V1.0R
   \date       01.03.94
   \brief      Include file for using plotting routines
   
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

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

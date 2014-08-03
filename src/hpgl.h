/************************************************************************/
/**

   \file       hpgl.h
   
   \version    V1.2
   \date       31.07.14
   \brief      Include file for hpgl
   
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
-  V1.0  25.03.91 Original
-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP

*************************************************************************/
#ifndef _HPGL_H
#define _HPGL_H

/************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "SysDefs.h"
#include "MathType.h"

/************************************************************************/
/* Prototypes
*/
BOOL blHPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin);
void blHPGLPen(int num);
void blHPGLMove(REAL x, REAL y);
void blHPGLDraw(REAL x, REAL y);
void blHPGLSetDash(int style);
void blHPGLFont(int font, REAL size);
void blHPGLLText(REAL x, REAL y, char *string);
void blHPGLCBText(REAL x, REAL y, REAL offset, char *text);
void blHPGLROffText(REAL x, REAL y, REAL offset, char *text);
void blHPGLLCText(REAL x, REAL y, char *text);
void blHPGLCTText(REAL x, REAL y, REAL offset, char *text);
void blHPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, 
                 REAL TitleSize, char *label, int LabelFont, REAL LabelSize);
void blHPGLEnd(void);
void blHPGLShowText(char *text, BOOL orientation, int XBase, int YBase);

/************************************************************************/
/* Deprecated functions: hpgl.h                                         */
/** \cond deprecated                                                    */

BOOL HPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin);
void HPGLPen(int num);
void HPGLMove(REAL x, REAL y);
void HPGLDraw(REAL x, REAL y);
void HPGLSetDash(int style);
void HPGLFont(int font, REAL size);
void HPGLLText(REAL x, REAL y, char *string);
void HPGLCBText(REAL x, REAL y, REAL offset, char *text);
void HPGLROffText(REAL x, REAL y, REAL offset, char *text);
void HPGLLCText(REAL x, REAL y, char *text);
void HPGLCTText(REAL x, REAL y, REAL offset, char *text);
void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, 
               REAL TitleSize, char *label, int LabelFont, REAL LabelSize);
void HPGLEnd(void);
void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase);

/* \endcond                                                             */
/************************************************************************/

#endif

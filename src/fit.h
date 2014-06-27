/************************************************************************/
/**

   \file       fit.h
   
   \version    V1.1R
   \date       01.03.94
   \brief      Include file for least squares fitting
   
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
-  V1.0  04.02.91 Original
-  V1.1  08.12.92 Removed qikfit() prototype as is static

*************************************************************************/
#ifndef _FIT_H
#define _FIT_H

#include "MathType.h"
#include "SysDefs.h"

/* Prototypes for functions defined in fit.c                            */
BOOL matfit(COOR *x1, COOR *x2, REAL rm[3][3], int n, REAL *wt1, 
            BOOL column);

#endif


/************************************************************************/
/**

   \file       MathUtil.h
   
   \version    V1.7
   \date       15.05.15
   \brief      Prototypes, etc. for maths utility routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2015
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
-  V1.0  01.03.94 Original
-  V1.1  18.06.96 Added vector routines
-  V1.2  10.09.96 Added combperm.c routines
-  V1.3  06.10.98 Added VecAdd3()
-  V1.4  07.07.14 Use bl prefix for functions By: CTP
-  V1.5  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.6  14.08.14 Moved deprecated function prototypes to deprecated.h 
                  By: CTP
-  V1.7  15.05.15 Added blIsPrime() and blFindNextPrime()  By: ACRM

*************************************************************************/
#ifndef _MATHUTIL_H
#define _MATHUTIL_H

#include <math.h>
#include "MathType.h"
#include "SysDefs.h"

/* Prototypes */
void blCalcSD(REAL val, int action, REAL *mean, REAL *SD);
void blCalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                 int *NValues, REAL *mean, REAL *SD);
REAL blPearson(REAL *x, REAL *y, int NItem);
REAL blPearson1(REAL *x, REAL *y, int NItem);

void blCrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2);
void blVecSub3(VEC3F *Out, VEC3F In1, VEC3F In2);
void blVecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2);
REAL blVecLen3(VEC3F Vec);
REAL blDistPtVect(VEC3F Point, VEC3F End1, VEC3F End2);
REAL blPointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac);
ULONG blFactorial(int n);
ULONG blFactdiv(int n1, int n2);
ULONG blNPerm(int n, int r);
ULONG blNComb(int n, int r);
BOOL  blIsPrime(ULONG input);
ULONG blFindNextPrime(ULONG input, BOOL above);


/************************************************************************/
/* Include deprecated functions                                         */
#define _MATHUTIL_H_DEPRECATED
# include "deprecated.h" 
/************************************************************************/

#endif

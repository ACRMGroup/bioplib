/************************************************************************/
/**

   \file       MathUtil.h
   
   \version    V1.3R
   \date       06.10.98
   \brief      Prototypes, etc. for maths utility routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-8
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

*************************************************************************/
#ifndef _MATHUTIL_H
#define _MATHUTIL_H

#include <math.h>
#include "MathType.h"
#include "SysDefs.h"

void CalcSD(REAL val, int action, REAL *mean, REAL *SD);
void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
               int *NValues, REAL *mean, REAL *SD);
REAL pearson(REAL *x, REAL *y, int NItem);
REAL pearson1(REAL *x, REAL *y, int NItem);

void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2);
void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2);
void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2);
REAL VecLen3(VEC3F Vec);
REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2);
REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac);
ULONG factorial(int n);
ULONG factdiv(int n1, int n2);
ULONG NPerm(int n, int r);
ULONG NComb(int n, int r);

#endif


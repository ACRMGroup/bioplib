/************************************************************************/
/**

   \file       matrix.h
   
   \version    V1.6R
   \date       27.09.95
   \brief      Include file for matrix operations
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1995
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
/* Includes
*/

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#ifndef _MATRIX_H
#define _MATRIX_H
#include "MathType.h"

void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout);
void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3]);
void invert33(REAL s[3][3], REAL ss[3][3]);
void CreateRotMat(char direction, REAL angle, REAL matrix[3][3]);
REAL VecDist(REAL *a, REAL *b, int len);

#endif

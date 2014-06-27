/************************************************************************/
/**

   \file       MathType.h
   
   \version    V1.0R
   \date       30.08.94
   \brief      Type definitions for maths
   
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
#ifndef _MATHTYPE_H
#define _MATHTYPE_H

/* This is for compilers running on machines such as Amigas, Macs and
   older Sun workstations using 680X0 series processors with maths
   coprocessors. This assumes that the symbol _M68881 is defined when
   the compiler is run to use the maths coprocessor and that a file
   called m68881.h is to be included to make full use of the coprocessor
*/
#ifdef _M68881
#include <m68881.h>
#endif

/* Note, that if this is changed to float, all I/O routines using type
   REAL will need %lf's changing to %f's
*/
typedef double REAL;

typedef struct
{  REAL x, y, z;
}  VEC3F;

typedef VEC3F COOR;

/* Define PI if not done                                                */
#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

#endif

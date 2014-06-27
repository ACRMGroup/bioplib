/************************************************************************/
/**

   \file       array.h
   
   \version    V1.5R
   \date       30.05.02
   \brief      Include file for 2D/3D array functions
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2002
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
-  V1.4  18.03.94
-  V1.5  30.05.02 Added 3D functions

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

#ifndef _ARRAY_H
#define _ARRAY_H

char **Array2D(int size, int dim1, int dim2);
void FreeArray2D(char **array, int dim1, int dim2);

char ***Array3D(int size, int dim1, int dim2, int dim3);
void FreeArray3D(char ***array, int dim1, int dim2, int dim3);

#endif

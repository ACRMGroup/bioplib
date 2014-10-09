/************************************************************************/
/**

   \file       WriteCrystPDB.c
   
   \version    V1.1
   \date       07.07.14
   \brief      
   
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
-  V1.0R 12.10.05 Original
-  V1.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include "MathType.h"

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
/*>void blWriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                      char *spacegroup,
                      REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
   -------------------------------------------------------------------
*//**

   \param[in]     *fp            Output file pointet
   \param[in]     UnitCell       The unit cell dimensions
   \param[in]     CellAngles     The unit cell angles
   \param[in]     *spacegroup    The crystal's space group
   \param[in]     OrigMatrix     The origin matrix
   \param[in]     ScaleMatrix    The scale matrix

   Write crystal parameters (unit cell, space group, origin and scale
   matrices) to a PDB file.

-  12.10.95 Original    By: ACRM
-  17.10.95 Corrected %lf to %f in fprintf()s
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blWriteCrystPDB(FILE *fp, VEC3F UnitCell, VEC3F CellAngles,
                   char *spacegroup,
                   REAL OrigMatrix[3][4], REAL ScaleMatrix[3][4])
{
   int i;
   
   CellAngles.x *= 180.0/PI;
   CellAngles.y *= 180.0/PI;
   CellAngles.z *= 180.0/PI;

   fprintf(fp,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %14s\n",
           UnitCell.x, UnitCell.y, UnitCell.z,
           CellAngles.x, CellAngles.y, CellAngles.z, spacegroup);
   for(i=0; i<3; i++)
   {
      fprintf(fp,"ORIGX%1d    %10.6f%10.6f%10.6f%15.5f\n", i+1,
              OrigMatrix[i][0],
              OrigMatrix[i][1],
              OrigMatrix[i][2],
              OrigMatrix[i][3]);
   }
   for(i=0; i<3; i++)
   {
      fprintf(fp,"SCALE%1d    %10.6f%10.6f%10.6f%15.5f\n", i+1,
              ScaleMatrix[i][0],
              ScaleMatrix[i][1],
              ScaleMatrix[i][2],
              ScaleMatrix[i][3]);
   }
}


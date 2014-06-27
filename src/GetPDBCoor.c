/************************************************************************/
/**

   \file       GetPDBCoor.c
   
   \version    V1.3
   \date       14.03.96
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-6
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
-  V1.0  01.03.94 Original release
-  V1.1  11.03.94 Fixed bug in calls to matfit(). Had not been changed 
                  to reflect modification in MatMult3_33().
-  V1.2  14.03.94 Fixed FitPDB(); wasn't filling in the output matrix
-  V1.3  14.03.96 Added FitCaPDB()
                  Changed FitPDB() and FitCaCbPDB() to use 
                  ApplyMatrixPDB() rather than RotatePDB() since the PDB
                  linked lists are already at the origin

*************************************************************************/
/* Includes
*/
#include <stdlib.h>

#include "pdb.h"
#include "macros.h"

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
/*>int GetPDBCoor(PDB *pdb, COOR **coor)
   -------------------------------------
   Input:   PDB  *pdb    PDB linked list
   Output:  COOR **coor  Array of coordinate structure pointers
   Returns: int          Number of coordinates copied

   Get the coordinates out of a PDB linked list into an array of type COOR
   The COOR array is allocated for you

   17.06.93 Original    By: ACRM
   22.06.93 Corrected return value
*/
int GetPDBCoor(PDB *pdb, COOR **coor)
{
   PDB   *p;
   int   NAtoms;
   
   /* Check valid pdb linked list                                       */
   if(pdb == NULL) return(0);
   
   /* Count the atoms in the PDB linked list                            */
   for(p=pdb,NAtoms=0; p!=NULL; NEXT(p))  NAtoms++;
   
   /* Allocate memory for coordinate array                              */
   *coor = (COOR *)malloc(NAtoms * sizeof(COOR));
   if(*coor == NULL) return(0);
   
   /* Copy the coordinates into the array                               */
   for(p=pdb,NAtoms=0; p!=NULL; NEXT(p),NAtoms++)
   {
      (*coor)[NAtoms].x = p->x;
      (*coor)[NAtoms].y = p->y;
      (*coor)[NAtoms].z = p->z;
   }

   return(NAtoms);
}



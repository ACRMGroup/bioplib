/************************************************************************/
/**

   \file       FitPDB.c
   
   \version    V1.4
   \date       07.07.14
   \brief      Fit two PDB linked lists. Also a weighted fit and support
               routines
   
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
-  V1.4  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "MathType.h"
#include "SysDefs.h"
#include "macros.h"
#include "fit.h"
#include "pdb.h"

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
/*>BOOL blFitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
   --------------------------------------------------------
*//**

   \param[in]     *ref_pdb     Reference PDB linked list
   \param[in,out] *fit_pdb     Mobile PDB linked list
   \param[out]    rm           Rotation matrix (May be input as NULL).
   \return                     Success

   Fits two PDB linked lists. Actually fits fit_pdb onto ref_pdb and also
   returns the rotation matrix. This may be NULL if these data are not
   required.

-  17.06.93 Original based on code from ProFit  By: ACRM
-  11.03.94 Changed call to matfit(). Corrected to normal matrix.
-  14.03.94 Actually fills in the rotation matrix (!). Restores original
            data if fitting failed.
-  14.03.96 Changed to use ApplyMatrixPDB() rather than RotatePDB() since
            we are already at the origin
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blFitPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3];
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_CofG,
         fit_CofG;
   int   NCoor       = 0,
         i, j;
   BOOL  RetVal;

   /* Get the CofG of the reference structure                           */
   blGetCofGPDB(ref_pdb, &ref_CofG);
   blGetCofGPDB(fit_pdb, &fit_CofG);

   /* Move them both to the origin                                      */
   blOriginPDB(ref_pdb);
   blOriginPDB(fit_pdb);
   
   /* Create coordinate arrays                                          */
   NCoor = blGetPDBCoor(ref_pdb, &ref_coor);
   if(blGetPDBCoor(fit_pdb, &fit_coor) != NCoor) 
   {
      /* Free and return if arrays don't match                          */
      if(ref_coor) free(ref_coor);
      if(fit_coor) free(fit_coor);
      return(FALSE);
   }
   
   /* Can't fit with fewer than 3 coordinates                           */
   if(NCoor < 3)
   {
      if(ref_coor) free(ref_coor);
      if(fit_coor) free(fit_coor);
      return(FALSE);
   }
   
   /* Everything OK, go ahead with the fitting                          */
   RetVal = blmatfit(ref_coor,fit_coor,RotMat,NCoor,NULL,FALSE);
   
   /* Now we can rotate the rotation list                               */
   if(RetVal)
   {
      blApplyMatrixPDB(fit_pdb, RotMat);
      blTranslatePDB(fit_pdb, ref_CofG);
      blTranslatePDB(ref_pdb, ref_CofG);
   }
   else
   {
      blTranslatePDB(fit_pdb, fit_CofG);
      blTranslatePDB(ref_pdb, ref_CofG);
   }

   /* Free the coordinate arrays                                        */
   free(ref_coor);
   free(fit_coor);

   /* Fill in the rotation matrix for output, if required               */
   if(rm!=NULL)
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }

   return(RetVal);
}


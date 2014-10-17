/************************************************************************/
/**

   \file       FitCaCbPDB.c
   
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
-  V1.4  07.07.14 Use bl prefix for functions. By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Fitting
   #FUNCTION  blFitCaCbPDB()
   Does a weighted fitting of 2 PDB linked lists. The CA and CB are given
   a weight of 1.0 while the other atoms are given a weight of 1.0/natom
   in the residue.
   Thus for N,CA,C,CB backbone only, this will be N and C with weights 
   of 0.25
*/
/************************************************************************/
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
/*>BOOL blFitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
   ------------------------------------------------------------
*//**

   \param[in]     *ref_pdb   Reference PDB linked list
   \param[in,out] *fit_pdb   PDB linked list to be fitted
   \param[out]    rm         Rotation matrix
   \return                   Success

   Does a weigthed fitting of 2 PDB linked lists. The CA and CB are given
   a weight of 1.0 while the other atoms are given a weight of 1.0/natom
   in the residue.
   Thus for N,CA,C,CB backbone only, this will be N and C with weights 
   of 0.25

-  21.06.93 Original    By: ACRM
-  22.06.93 Added i increment (!). Corrected return of rm
-  11.03.94 Changed call to matfit(). Corrected to normal matrix.
-  15.03.94 Defines the rotation matrix locally and only copies for
            output if rm is not NULL
-  14.03.96 Changed to use ApplyMatrixPDB() rather than RotatePDB() since
            we are already at the origin
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blFitCaCbPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3],
         *weights    = NULL;
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_CofG;
   int   NCoor       = 0,
         natom       = 0,
         resnum      = 0,
         i           = 0,
         j           = 0;
   char  insert      = ' ';
   BOOL  RetVal      = TRUE;
   PDB   *p          = NULL,
         *q          = NULL,
         *start      = NULL;

   /* Get the CofG of the reference structure                           */
   blGetCofGPDB(ref_pdb, &ref_CofG);

   /* Move them both to the origin                                      */
   blOriginPDB(ref_pdb);
   blOriginPDB(fit_pdb);
   
   /* Create coordinate arrays                                          */
   NCoor = blGetPDBCoor(ref_pdb, &ref_coor);
   if((NCoor == 0) || (blGetPDBCoor(fit_pdb, &fit_coor) != NCoor))
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

   /* Allocate memory for weigths array                                 */
   if((weights = (REAL *)malloc(NCoor * sizeof(REAL))) == NULL)
   {
      if(ref_coor) free(ref_coor);
      if(fit_coor) free(fit_coor);
      return(FALSE);
   }

   /* Create the weight array                                           */
   resnum = ref_pdb->resnum;
   insert = ref_pdb->insert[0];
   start  = ref_pdb;
   natom  = 0;
   i      = 0;

   for(p=ref_pdb; p!=NULL; NEXT(p))
   {
      if(p->resnum != resnum || p->insert[0] != insert)
      {
         /* Start of next residue                                       */
         for(q=start; q!=p; NEXT(q))
         {
            weights[i] = (REAL)(1.0/(REAL)natom);
            if(!strncmp(q->atnam,"CA  ",4)) weights[i] = (REAL)1.0;
            if(!strncmp(q->atnam,"CB  ",4)) weights[i] = (REAL)1.0;
            i++;
         }
         natom = 0;
         start = p;
      }
      natom++;
   }

   /* End of structure                                                  */
   for(q=start; q!=NULL; NEXT(q))
   {
      weights[i] = (REAL)(1.0/(REAL)natom);
      if(!strncmp(q->atnam,"CA  ",4)) weights[i] = (REAL)1.0;
      if(!strncmp(q->atnam,"CB  ",4)) weights[i] = (REAL)1.0;
      i++;
   }
   
   /* Everything OK, go ahead with the fitting                          */
   RetVal = blMatfit(ref_coor,fit_coor,RotMat,NCoor,weights,FALSE);
   
   /* Now we can rotate the rotation list                               */
   if(RetVal)
      blApplyMatrixPDB(fit_pdb, RotMat);

   /* Translate both structures back to reference CofG                  */
   blTranslatePDB(fit_pdb, ref_CofG);
   blTranslatePDB(ref_pdb, ref_CofG);

   /* Free the coordinate and weights arrays                            */
   free(ref_coor);
   free(fit_coor);
   free(weights);

   /* Fill in the rotation matrix for output, if required               */
   if(rm!=NULL)
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }

   return(RetVal);
}


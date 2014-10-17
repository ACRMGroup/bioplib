/************************************************************************/
/**

   \file       FitCaPDB.c
   
   \version    V1.6
   \date       19.08.14
   \brief      Fit two PDB linked lists. Also a weighted fit and support
               routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2009
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
   - V1.0 01.03.94 Original release
   - V1.1 11.03.94 Fixed bug in calls to matfit(). Had not been
                   changed to reflect modification in MatMult3_33().
   - V1.2 14.03.94 Fixed FitPDB(); wasn't filling in the output matrix
   - V1.3 14.03.96 Added FitCaPDB() Changed FitPDB() and FitCaCbPDB()
                   to use ApplyMatrixPDB() rather than RotatePDB()
                   since the PDB linked lists are already at the
                   origin
   - V1.4 28.01.09 Initialize RetVal - this randomly worked in 32bit
                   but broke in 64bit
   - V1.5 07.07.14 Use bl prefix for functions By: CTP
   - V1.6 19.08.14 Fixed calls to renamed function:
                   blSelectAtomsPDBAsCopy() By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Fitting
   #FUNCTION  blFitCaPDB()
   Fits two PDB linked lists using only the CA atoms. 
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
/*>BOOL blFitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
   --------------------------------------------------------
*//**

   \param[in]     *ref_pdb     Reference PDB linked list
   \param[in,out] *fit_pdb     Mobile PDB linked list
   \param[out]    rm           Rotation matrix (May be input as NULL).
   \return                      Success

   Fits two PDB linked lists using only the CA atoms. 

   Actually fits fit_pdb onto ref_pdb and also returns the rotation 
   matrix. This may be NULL if these data are not required.

-  14.03.96 Original based on FitPDB()   By: ACRM
-  28.01.09 Initialize RetVal to TRUE!
-  07.07.14 Use bl prefix for functions By: CTP
-  19.08.14 Added AsCopy suffix to calls to blSelectAtomsPDB() By: CTP
*/
BOOL blFitCaPDB(PDB *ref_pdb, PDB *fit_pdb, REAL rm[3][3])
{
   REAL  RotMat[3][3];
   COOR  *ref_coor   = NULL,
         *fit_coor   = NULL;
   VEC3F ref_ca_CofG,
         fit_ca_CofG,
         tvect;
   int   NCoor       = 0,
         i, j,
         natoms;
   BOOL  RetVal = TRUE;
   PDB   *ref_ca_pdb = NULL,
         *fit_ca_pdb = NULL;
   char  *sel[2];

   /* First extract only the CA atoms                                   */
   SELECT(sel[0], "CA  ");
   if(sel[0]==NULL)
      return(FALSE);
   if( (ref_ca_pdb = blSelectAtomsPDBAsCopy(ref_pdb, 1, sel, &natoms))
       == NULL )
      RetVal = FALSE;
   if( (fit_ca_pdb = blSelectAtomsPDBAsCopy(fit_pdb, 1, sel, &natoms))
       == NULL )
      RetVal = FALSE;
   free(sel[0]);
   
   /* If we succeeded in building our CA PDB linked lists...            */
   if(RetVal)
   {
      /* Get the CofG of the CA structures and the original mobile      */
      blGetCofGPDB(ref_ca_pdb, &ref_ca_CofG);
      blGetCofGPDB(fit_ca_pdb, &fit_ca_CofG);
      
      /* Move them both to the origin                                   */
      blOriginPDB(ref_ca_pdb);
      blOriginPDB(fit_ca_pdb);
      
      /* Create coordinate arrays, checking numbers match               */
      NCoor = blGetPDBCoor(ref_ca_pdb, &ref_coor);
      if(blGetPDBCoor(fit_ca_pdb, &fit_coor) != NCoor)
      {
         RetVal = FALSE;
      }
      else
      {
         /* Can't fit with fewer than 3 coordinates                     */
         if(NCoor < 3)
         {
            RetVal = FALSE;
         }
         else
         {
            /* Everything OK, go ahead with the fitting                 */
            if(!blMatfit(ref_coor,fit_coor,RotMat,NCoor,NULL,FALSE))
            {
               RetVal = FALSE;
            }
            else
            {
               /* Apply the operations to the true coordinates          */
               tvect.x = (-fit_ca_CofG.x);
               tvect.y = (-fit_ca_CofG.y);
               tvect.z = (-fit_ca_CofG.z);
               blTranslatePDB(fit_pdb, tvect);
               blApplyMatrixPDB(fit_pdb, RotMat);
               blTranslatePDB(fit_pdb, ref_ca_CofG);
            }
         }
      }
   }
   
   /* Free the coordinate arrays and CA PDB linked lists                */
   if(ref_coor)   free(ref_coor);
   if(fit_coor)   free(fit_coor);
   if(ref_ca_pdb) FREELIST(ref_ca_pdb, PDB);
   if(fit_ca_pdb) FREELIST(fit_ca_pdb, PDB);
         
   /* Fill in the rotation matrix for output, if required               */
   if(RetVal && (rm!=NULL))
   {
      for(i=0; i<3; i++)
         for(j=0; j<3; j++)
            rm[i][j] = RotMat[i][j];
   }

   return(RetVal);
}


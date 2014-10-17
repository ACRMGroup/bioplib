/************************************************************************/
/**

   \file       CalcRMSPDB.c
   
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
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Calculations
   #FUNCTION  blCalcRMSPDB()
   Calculate the RMS deviation between two pre-fitted PDB linked lists.
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
/*>REAL blCalcRMSPDB(PDB *pdb1, PDB *pdb2)
   ---------------------------------------
*//**

   \param[in]     *pdb1   First PDB linked list
   \param[in]     *pdb2   Second PDB linked list
   \return                RMS deviation

   Calculate the RMS deviation between two fitted PDB linked lists. 
   No fitting is done. The two lists must contain equivalent 
   structures (same atom types in same order). No checks are made 
   on this.

-  11.03.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
REAL blCalcRMSPDB(PDB *pdb1, PDB *pdb2)
{
   PDB *p, *q;
   int count = 0;
   REAL dist = (REAL)0.0,
        rms;
   
   for(p=pdb1, q=pdb2; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
   {
      dist += DISTSQ(p,q);
      count++;
   }

   rms = (REAL)((count)?sqrt((double)(dist/(REAL)count)):0.0);
   return(rms);
}



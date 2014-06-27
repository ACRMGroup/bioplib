/************************************************************************/
/**

   \file       NPerm.c
   
   \version    V1.0
   \date       10.09.96
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996
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
-  V1.0  10.09.96 Original

*************************************************************************/
/* Includes
*/
#include "SysDefs.h"
#include "MathUtil.h"

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
/*>ULONG NPerm(int n, int r)
   -------------------------
   Calculates number of permutations of n items in r groups
   Returns 0 if a numeric overflow occurs.

   09.09.96 Original   By: ACRM
*/
ULONG NPerm(int n, int r)
{
   return(factdiv(n,(n-r)));
}



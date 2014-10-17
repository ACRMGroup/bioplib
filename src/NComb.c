/************************************************************************/
/**

   \file       NComb.c
   
   \version    V1.1
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996-2014
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
-  V1.1  07.07.14 Use bl prefix for functions By: CTP


*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Miscellaneous
   #ROUTINE  blNComb()
   Calculates number of combinations of n items in r groups
*/
/************************************************************************/
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
/*>ULONG blNComb(int n, int r)
   ---------------------------
*//**

   Calculates number of combinations of n items in r groups
   Returns 0 if a numeric overflow occurs.
   
-  09.09.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
ULONG blNComb(int n, int r)
{
   ULONG f;
   f = blFactorial(r);
   
   return((f>0L)?(blNPerm(n,r)/blFactorial(r)):0L);
}

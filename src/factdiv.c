/************************************************************************/
/**

   \file       factdiv.c
   
   \version    V1.1
   \date       07.07.14
   \brief      
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1996-2014
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
/* Includes
*/
#include "SysDefs.h"

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
/*>ULONG blFactdiv(int n1, int n2)
   -------------------------------
*//**

   Calculates the factorial of one number divided by the factorial of
   another (smaller) number.
   Returns 0 on numeric overflow or if n2 > n1

-  09.09.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
ULONG blFactdiv(int n1, int n2)
{
   int   i;
   ULONG ret  = 1L,
         prev = 0L;

   if(n2 > n1)
      return(0);
   
   for(i=n2+1; i<=n1; i++)
   {
      ret *= (ULONG)i;
      if(ret < prev)
         return(0);
      prev = ret;
   }
   
   return(ret);
}



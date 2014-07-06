/************************************************************************/
/**

   \file       CalcSD.c
   
   \version    V1.3R
   \date       22.06.94
   \brief      Calculate mean and standard deviation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1990-3
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

   This routine calculates the mean and standard deviation from a set
   of numbers. The routine is called with each value to be sampled
   and the action required is specified.

**************************************************************************

   Usage:
   ======


      #include "MathUtil.h" before using these routines

      CalcSD(val,action,mean,SD)
      --------------------------
      Input:   val     int       The value to be sampled
      Input:   action  short     0: Sample the value
                                 1: Calculate & return mean and SD
                                 2: Clear the sample lists
      Output:  mean    *REAL     The returned mean
      Output:  SD      *REAL     The returned standard deviation

      The output values are only set when action=1


**************************************************************************

   Revision History:
   =================
-  V1.0  30.03.90 Original    By: ACRM
-  V1.1  17.06.93 Modified for book
-  V1.2  01.03.94 Added CalcExtSD()
-  V1.3  22.06.94 Fixed for just one value

*************************************************************************/
/* Includes
*/
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
/*>void CalcSD(REAL val, int action, REAL *mean, REAL *SD)
   -------------------------------------------------------
*//**

   Calculate the mean and standard deviation from a set of numbers. 
   The routine is called with each value to be sampled and the action 
   required is specified:

   \param[in]     val          The value to be sampled
   \param[in]     action       0: Sample the value
                               1: Calculate & return mean and SD
                               2: Clear the sample lists
   \param[out]    mean         The returned mean
   \param[out]    SD           The returned standard deviation

   The output values are only set when action=1

-  30.03.90 Original    By: ACRM
-  17.06.93 Modified for book
-  22.06.94 Fixed for 0 values supplied
*/
void CalcSD(REAL val, int action, REAL *mean, REAL *SD)
{
   static REAL  SxSq = 0.0,
                Sx   = 0.0;
   static int   NSDValues = 0;
    
   switch(action)
   {
   case 0:
       NSDValues++;
       SxSq += (val * val);
       Sx   += val;
       break;
        
   case 1:
       *mean = *SD = (REAL)0.0;
       if(NSDValues > 0)
          *mean = Sx / NSDValues;
       if(NSDValues > 1)
          *SD   = sqrt((SxSq - (Sx * Sx)/NSDValues)/(NSDValues - 1));
       break;
        
   case 2:
       SxSq = 0.0;
       Sx   = 0.0;
       NSDValues = 0;
   }
}


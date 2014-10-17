/************************************************************************/
/**

   \file       CalcExtSD.c
   
   \version    V1.4
   \date       07.07.14
   \brief      Calculate mean and standard deviation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1990-2014
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

      void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                     int *NValues, REAL *mean, REAL *SD)
      ----------------------------------------------------------
      Input:   val     int       The value to be sampled
      Input:   action  short     0: Sample the value
                                 1: Calculate & return mean and SD
                                 2: Clear the sample lists
      Output:  mean    *REAL     The returned mean
      Output:  SD      *REAL     The returned standard deviation
      I/O:     Sx      *REAL     Sum of values
      I/O:     SxSq    *REAL     Sum of values squared
      I/O:     NValues *int      Number of values

      The output values are only set when action==1

**************************************************************************

   Revision History:
   =================
-  V1.0  30.03.90 Original    By: ACRM
-  V1.1  17.06.93 Modified for book
-  V1.2  01.03.94 Added CalcExtSD()
-  V1.3  22.06.94 Fixed for just one value
-  V1.4  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Statistics
   #FUNCTION  blCalcExtSD()
   Calculate the mean and standard deviation from a set of numbers. 
   The routine is called with each value to be sampled and then again
   to obtain the results.
   This is the same as blCalcSD() except that the variables used to
   accumulate totals are kept outside the function instead of being 
   static within the function
*/
/************************************************************************/
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
/*>void blCalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                  int *NValues, REAL *mean, REAL *SD)
   ------------------------------------------------------------
*//**

   \param[in]     val          The value to be sampled
   \param[in]     action       0: Sample the value
                               1: Calculate & return mean and SD
                               2: Clear the sample lists
   \param[out]    mean         The returned mean
   \param[out]    SD           The returned standard deviation
   \param[in,out] Sx           Sum of values
   \param[in,out] SxSq         Sum of values squared
   \param[in,out] NValues      Number of values

   Calculate the mean and standard deviation from a set of numbers. 
   The routine is called with each value to be sampled and the action 
   required is specified:

   The output values are only set when action==1

   This is the same as blCalcSD() except that the Sx, SxSq and NValues
   variables are kept outside the function instead of being static
   within the function

-  13.10.93 Original based on CalcSD   By: ACRM
-  22.06.94 Fixed for only one value supplied
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blCalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, 
                 int *NValues, REAL *mean, REAL *SD)
{
   switch(action)
   {
   case 0:
       (*NValues)++;
       *SxSq += (val * val);
       *Sx   += val;
       break;
        
   case 1:
       *mean = *SD = (REAL)0.0;
       if(*NValues > 0)
          *mean = (*Sx) / (*NValues);
       if(*NValues > 1)
          *SD   = sqrt((*SxSq - ((*Sx) * (*Sx)) / (*NValues)) /
                       (*NValues - 1));
       break;
        
   case 2:
       *SxSq = 0.0;
       *Sx   = 0.0;
       *NValues = 0;
   }
}


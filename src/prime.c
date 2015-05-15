/************************************************************************/
/**

   \file       prime.c
   
   \version    V1.0
   \date       15.05.15
   \brief      Routines to identify prime numbers
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2015
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
-  V1.0   15.05.15  Original   By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP Maths

   #FUNCTION blIsPrime()
   Tests whether a number is prime

   #FUNCTION blFindNextPrime()
   Finds the next prime

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
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
/*>ULONG blFindNextPrime(ULONG primeNum, BOOL above)
----------------------------------------------------
*//**
   \param[in]    primeNum   Input number (prime or not)
   \param[in]    above      Always return a number above the input
   \return                  The next prime

   Returns a prime above (or equal to if 'above' is FALSE) the input
   number

-  15.05.15  Original   By: ACRM
*/
ULONG blFindNextPrime(ULONG primeNum, BOOL above)
{
   /* Check we have some chance of finding a prime
      Largest 64-bit ULONG prime is 
         18446744073709551557 = ULONG_MAX - 58
      Largest 32-bit ULONG prime is
         4294967291 = ULONG - 4
    */
   if(ULONG_MAX == 4294967295UL)   /* 32-bit                            */
   {
      if(primeNum > (ULONG_MAX-4))
      return(0);
   }
   else                            /* 64-bit                            */
   {
      if(primeNum > (ULONG_MAX-59))
      return(0);
   }

   if(above)
      primeNum++;
 
   /* If the number is even and not 2, make it odd                      */
   if(((primeNum%2) == 0) && (primeNum != 2))
      primeNum++;
   
   /* Look at each odd number until we find a prime                     */
   while(!blIsPrime(primeNum))
      primeNum+=2;

   return(primeNum);
}
 

/************************************************************************/
/*>BOOL blIsPrime(ULONG input)
------------------------------
*//**
   \param[in]   input    Number to test
   \return               Is this number a prime?

   Tests whether the input number is a prime

-  15.05.15  Original   By: ACRM
*/
BOOL blIsPrime(ULONG input)
{
   BOOL prime = TRUE;

   /* numbers between 
         4,294,967,293 and 4,294,967,295 (32bit) 
      and
         18,446,744,073,709,551,613 and 18,446,744,073,709,551,615 (64bit)
      are not prime
   */

   if(input > (ULONG_MAX-2))
      return(FALSE);
   
   if(input == 2)
      return(TRUE);
   
   /* Return if it's 0, 1 or even                                       */
   if(input%2 == 0 || input <= 1)
   {
      prime = FALSE;
   } 
   else 
   {
      ULONG i;

      /* Test by dividing by each odd number from 3 up to the square
         root of the input number
      */
      for(i=3; i<=(ULONG)sqrt((double)input); i+=2)
      {
         if(input%i == 0)
         {
            prime = FALSE;
            break;
         }
      }
   }
   return prime;
}

#ifdef TEST
/************************************************************************/
int main(int argc, char **argv)
{
    ULONG input =  1;
    ULONG nextPrime;
 
    if(argc!=2)
       return(1);

    if(sscanf(argv[1], "%lu", &input))
    {       
        nextPrime = blFindNextPrime(input, FALSE);
 
        if(input == nextPrime)
        {
           printf("The input number was prime %lu\n", 
                  input);
        }
        else if(nextPrime == 0)
        {
           printf("The input number was out of range %lu\n", 
                  input);
        }
        else
        {
           printf("The prime number following %lu is %lu\n", 
                  input, nextPrime);
        }
    }
    return(0);
}
#endif

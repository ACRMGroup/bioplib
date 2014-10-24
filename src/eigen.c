/************************************************************************/
/**

   \file       eigen.c
   
   \version    V1.0
   \date       03.10.14
   \brief      Calculates Eigen values and Eigen vectors for a 
               symmetric matrix
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2014
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
-  V1.0   03.10.14   Original

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Matrices
   #FUNCTION blEigen()
   Calculates the eigenvalues and eigenvectors of a REAL symmetric matrix
   Note that this routine destroys the values above the diagonal of the
   matrix.
*/

/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <math.h>
#include "MathType.h"
#include "array.h"
#include "eigen.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXITERATION 50  /* Maximum number of iterations                */

/* Test whether a number x is smaller than double precision can cope with
   when added to a number y
*/
#define TESTSMALL(x, y) ((fabs(y)+(x)) == fabs(y))

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void PerformJacobiRotation(int ip, int iq, REAL g, int n, 
                                  REAL **matrix, REAL **eigenVectors, 
                                  REAL *eigenValues, REAL *ta_pq);


/************************************************************************/
/*>int blEigen(REAL **matrix, REAL **eigenVectors, REAL *eigenValues, 
               int matrixSize)
   ------------------------------------------------------------------
*//**
   \param[in]  **matrix         Symmetric matrix
   \param[in]  matrixSize       Dimension of the matrix
   \param[out] **eigenVectors   The eigen vectors
   \param[out] *eigenValues     The eigen values
   \return                      If >0, the number of Jacobi rotations
                                performed
                                If <0, error:
                                EIGEN_NOMEMORY - memory allocation
                                EIGEN_NOCONVERGE - didn't converge
                                   in 50 iterations

   Calculates the eigenvalues and eigenvectors of a REAL symmetric matrix
   Note that this routine destroys the values above the diagonal of the
   matrix.

-  02.10.14  Original   By: ACRM
*/
int blEigen(REAL **matrix, REAL **eigenVectors, REAL *eigenValues, 
            int matrixSize)

{
   int  column         = 0, 
        row            = 0, 
        iteration      = 0, 
        nRotations     = 0;
   REAL threshold      = 0.0, 
        offDiagonalSum = 0.0, 
        *diagonal      = NULL, 
        *ta_pq         = NULL;

   diagonal = (REAL *)malloc(matrixSize*sizeof(REAL));
   ta_pq    = (REAL *)malloc(matrixSize*sizeof(REAL));

   if((diagonal == NULL) || (ta_pq == NULL))
   {
      if(diagonal != NULL) free(diagonal);
      if(ta_pq    != NULL) free(ta_pq);
      return(EIGEN_NOMEMORY);
   }

   /* Initialize the eigenVectors matrix to the identity matrix.        */
   for(row=0; row<matrixSize; row++) 
   {
      for(column=0; column<matrixSize; column++)
         eigenVectors[row][column]=0.0;
      eigenVectors[row][row]=1.0;
   }

   /* Initialize diagonal and eigenValues to the diagonal of matrix and
      initialize ta_{pq} to accumulate terms of the form in NumRec
      Equation 11.1.14  
   */
   for(row=0; row<matrixSize; row++) 
   {
      diagonal[row]    = matrix[row][row];
      eigenValues[row] = matrix[row][row];
      ta_pq[row]       = 0.0;
   }

   /* Iterate until the off-diagonal sum is reduced to zero             */
   for(iteration=0; iteration<MAXITERATION; iteration++) 
   {
      /* Sum the off-diagonal elements as in NumRec Equation 11.1.19    */
      offDiagonalSum = 0.0;
      for(row=0; row<matrixSize-1; row++) 
      {
         for(column=row+1; column<matrixSize; column++)
         {
            offDiagonalSum += fabs(matrix[row][column]);
         }
      }

      /* If we have converged such that the off-diagonal sum is zero, then
         we have our solution, so we return
      */
      if(offDiagonalSum == 0.0) 
      {
         free(ta_pq);
         free(diagonal);
         return(nRotations);
      }

      /* Use a larger threshold for the first 4 iterations. This is 
         NumRec Equation 11.1.25
      */
      if(iteration < 4)
      {
         threshold = 0.2 * offDiagonalSum / (matrixSize*matrixSize);
      }
      else
      {
         threshold = 0.0;
      }
      
      for(row=0; row<matrixSize-1; row++) 
      {
         for(column=row+1; column<matrixSize; column++) 
         {
            REAL offDiagonal = 100.0 * fabs(matrix[row][column]);

            /* After the first four iterations, we skip the rotation if 
               the off-diagonal element is small
            */
            if(iteration > 4                            && 
               TESTSMALL(offDiagonal, eigenValues[row]) &&
               TESTSMALL(offDiagonal, eigenValues[column]))
            {
               matrix[row][column]=0.0;
            }
            else if(fabs(matrix[row][column]) > threshold) 
            {
               PerformJacobiRotation(row, column, offDiagonal, matrixSize,
                                     matrix, eigenVectors, eigenValues,
                                     ta_pq);

               /* Increment the number of rotations                     */
               nRotations++;
            }
         }
      }

      /* Update eigenValues with the value of ta_{pq}, and 
         reinitialize ta_{pq}
      */
      for(row=0; row<matrixSize; row++) 
      {
         diagonal[row]    += ta_pq[row];
         eigenValues[row]  = diagonal[row];
         ta_pq[row]        = 0.0;
      }
   }

   free(ta_pq);
   free(diagonal);
   return(EIGEN_NOCONVERGE);
}


/************************************************************************/
/*>static void PerformJacobiRotation(int row, int column, REAL g, 
                                     int matrixSize, 
                                     REAL **matrix, REAL **eigenVectors, 
                                     REAL *eigenValues, REAL *ta_pq)
   ---------------------------------------------------------------------
*//**
   \param[in]     row             Index of matrix row
   \param[in]     column          Index of matrix column
   \param[in]     offDiagonal     Off-diagonal value
   \param[in]     matrixSize      Size of the matrix     
   \param[in,out] matrix          The REAL matrix
   \param[in,out] eigenVectors    The eigen vectors
   \param[in,out] eigenValues     The eigen values
   \param[in,out] ta_pq           ta_{pq}

   Does the actual work of performing a Jacobi rotation

-  03.10.14  Original   By: ACRM
*/
static void PerformJacobiRotation(int row, int column, REAL offDiagonal, 
                                  int matrixSize, 
                                  REAL **matrix, REAL **eigenVectors, 
                                  REAL *eigenValues, REAL *ta_pq)
{
   REAL deltaEigenValue = eigenValues[column] - eigenValues[row];
   REAL tTimesElement,
        t,
        theta, 
        tau, 
        tOverSqrtTSq, 
        oneOverSqrtTSq;
   int  j;

   if(TESTSMALL(offDiagonal, deltaEigenValue))
   {
      t = (matrix[row][column]) / deltaEigenValue;
   }
   else 
   {
      /* Calculate theta as defined in NumRec Equation 11.1.8           */
      theta = 0.5 * deltaEigenValue / matrix[row][column];
      /* Calculate t as defined in NumRec Equation 11.1.10              */
      t     = 1.0 / (fabs(theta) + sqrt(1.0+theta*theta));
      if(theta < 0.0) t = -t;
   }
   
   tTimesElement        = t * matrix[row][column];

   ta_pq[row]          -= tTimesElement;
   ta_pq[column]       += tTimesElement;
   eigenValues[row]    -= tTimesElement;
   eigenValues[column] += tTimesElement;
   matrix[row][column]  = 0.0;

   oneOverSqrtTSq       = 1.0 / sqrt(1 + (t*t));
   tOverSqrtTSq         = t * oneOverSqrtTSq;
   tau                  = tOverSqrtTSq / (1.0 + oneOverSqrtTSq);
   
   /* Perform rotations for 0<=j<row                                    */
   for(j=0; j<row; j++) 
   { 
      REAL temp1        = matrix[j][row];
      REAL temp2        = matrix[j][column];
      matrix[j][row]    = temp1 - tOverSqrtTSq*(temp2 + temp1*tau);
      matrix[j][column] = temp2 + tOverSqrtTSq*(temp1 - temp2*tau);
   }
   
   /* Perform rotations for row<j<column                                */
   for(j=row+1; j<column; j++) 
   { 
      REAL temp1        = matrix[row][j];
      REAL temp2        = matrix[j][column];
      matrix[row][j]    = temp1 - tOverSqrtTSq*(temp2 + temp1*tau);
      matrix[j][column] = temp2 + tOverSqrtTSq*(temp1 - temp2*tau);
   }
   
   /* Perform rotations for column<j<matrixSize                         */
   for(j=column+1; j<matrixSize; j++) 
   { 
      REAL temp1        = matrix[row][j];
      REAL temp2        = matrix[column][j];
      matrix[row][j]    = temp1 - tOverSqrtTSq*(temp2 + temp1*tau);
      matrix[column][j] = temp2 + tOverSqrtTSq*(temp1 - temp2*tau);
   }
   
   /* And rotate the eigen vectors                                      */
   for(j=0; j<matrixSize; j++) 
   {
      REAL temp1              = eigenVectors[j][row];
      REAL temp2              = eigenVectors[j][column];
      eigenVectors[j][row]    = temp1 - tOverSqrtTSq*(temp2 + temp1*tau);
      eigenVectors[j][column] = temp2 + tOverSqrtTSq*(temp1 - temp2*tau);
   }
}

/************************************************************************/
/**

   \file       regression.c
   
   \version    V1.0
   \date       08.10.14
   \brief      Routines for finding the best fit line through a set of 
               points
   
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

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
-  V1.0   08.10.14  Original   By: ACRM  Based on code by Abhi Raghavan

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Matrices
   #FUNCTION blCalculateCovarianceMatrix()
   Find the covariance matrix for a given matrix
 
   #GROUP    Maths
   #SUBGROUP Geometry (3D)
   #FUNCTION blCalculateBestFitLine()
   Calculates a best fit line through a set of coordinates in 3D. 
   Results are returned as a vector and a point through which the vector
   must pass.

   #FUNCTION blFindCentroid()
   Calculates the centroid of an array of points. (Like blGetCofGPDB() 
   but for a coordinate array instead of a PDB linked list)
*/

/************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>
#include "pdb.h"
#include "array.h"
#include "macros.h"
#include "MathType.h"
#include "SysDefs.h"
#include "regression.h"
#include "eigen.h"

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
/*>BOOL blCalculateBestFitLine(REAL **coordinates, int numberOfPoints,
                               int numberOfDimensions, REAL *centroid,
                               REAL *eigenVector)
   -----------------------------------------------------------------
*//**
   \param[in]   coordinates        2D Array containing the coordinates
   \param[in]   numberOfPoints     The number of points in the array
   \param[in]   numberOfDimensions The number of dimensions (usually 3)
   \param[out]  centroid           The centroid of the points
   \param[out]  eigenVector        The best Eigen vector 
   \return                         Success?

   Calculates a best fit line through a set of coordinates in 3D. 
   The results are returned in 'eigenVector'

   The code calculates the centroid and ouputs this. The covariance 
   matrix is then calculated and the Eigen vectors and values of the 
   covariance matrix are then calculated. The Eigenvalue with the largest
   Eigenvector represents the best fit line passing through the centroid.

-  07.10.14 Original   By: ACRM  Based on code by Abhi Raghavan
*/
BOOL blCalculateBestFitLine(REAL **coordinates, int numberOfPoints,
                            int numberOfDimensions, REAL *centroid,
                            REAL *eigenVector)
{
   REAL *eigenValues           = NULL,
        **covarianceMatrix     = NULL,
        **eigenVectorMatrix    = NULL;
   BOOL retValue               = TRUE;
   int  largestEigenValueIndex = 0,
        i                      = 0;

   /* Allocate memory                                                   */
   eigenValues = (REAL *)malloc(numberOfDimensions * sizeof(REAL));
   covarianceMatrix = (REAL **)Array2D(sizeof(REAL), 
                                       numberOfDimensions,
                                       numberOfDimensions);
   eigenVectorMatrix = (REAL **)Array2D(sizeof(REAL), 
                                        numberOfDimensions,
                                        numberOfDimensions);
   
   if((eigenValues == NULL) || (eigenVectorMatrix == NULL) ||
      (covarianceMatrix == NULL))
   {
      retValue = FALSE;
   }
   else
   {
      /* Find the centroid of the coordinates                           */
      blFindCentroid(coordinates, numberOfPoints, 3, centroid);
      
      /* Calculate the covariance matrix                                */
      if(!blCalculateCovarianceMatrix(coordinates, numberOfPoints, 
                                      numberOfDimensions, 
                                      covarianceMatrix))
      {
         retValue = FALSE;
      }
      else
      {
         REAL largestEigenValue;

         /* Find the eigen values and vectors.                          */
         if(blEigen(covarianceMatrix, eigenVectorMatrix, eigenValues,
                    numberOfDimensions) < 0)
         {
            retValue = FALSE;
         }
         else
         {
            /* Find the eigen vector corresponding to the largest Eigen 
               value.
            */
            largestEigenValue      = eigenValues[0];
            largestEigenValueIndex = 0;
            for(i=1; i<numberOfDimensions; i++)
            {
               if(eigenValues[i] > largestEigenValue)
               {
                  largestEigenValue      = eigenValues[i];
                  largestEigenValueIndex = i;
               }
            }
            
            for(i=0;i<numberOfDimensions;i++)
            {
               eigenVector[i] =
                  eigenVectorMatrix[i][largestEigenValueIndex];
            }
         }
      }
   }
   
   /* Free allocated memory                                             */
   if(eigenValues)       free(eigenValues);
   if(eigenVectorMatrix) FreeArray2D((char **)eigenVectorMatrix, 
                                     numberOfDimensions, 
                                     numberOfDimensions);
   if(covarianceMatrix)  FreeArray2D((char **)covarianceMatrix, 
                                     numberOfDimensions, 
                                     numberOfDimensions);

   return(retValue);
}


/************************************************************************/
/*>void blFindCentroid(REAL **matrix, int numX, int numY,
                       REAL *centroid)
   ------------------------------------------------------
*//**
   \param[in]  matrix      Coordinate array
   \param[in]  numX        The x-dimension of the array
   \param[in]  numY        The y-dimension of the array
   \param[out] centroid    The centroid of the points

   Calculates the centroid of an array of points

-  07.10.14  Original   By: ACRM
*/
void blFindCentroid(REAL **matrix, int numX, int numY, 
                    REAL *centroid)
{
   int point,
       dim;

   for(dim=0; dim<numY; dim++)
      centroid[dim] = (REAL)0.0;
   
   for(point=0; point<numX; point++)
   {
      for(dim=0; dim<3; dim++)
      {
         centroid[dim] += matrix[point][dim];
      }
   }

   for(dim=0; dim<3; dim++)
   {
      centroid[dim] /= numX;
   }
}


/************************************************************************/
/*>BOOL blCalculateCovarianceMatrix(REAL **matrix, int numX, int numY,
                                    REAL **cov)
   -----------------------------------------------------------------
*//**
   \param[in]   matrix           The matrix
   \param[in]   numX             x-dimension of the matrix
   \param[in]   numY             y-dimension of the matrix
   \param[out]  cov              The covariance matrix
   \return                       Success

   Find the covariance matrix for a given matrix

-  07.10.14 Original   By: ACRM  (based on code by Abhi Raghavan)
*/
BOOL blCalculateCovarianceMatrix(REAL **matrix, int numX, int numY, 
                                 REAL **cov)
{
   int  i     = 0,
        j     = 0,
        k     = 0;
   REAL total = 0.0,
        *centroid;

   if((centroid=(REAL *)malloc(numY * sizeof(REAL)))==NULL)
      return(FALSE);

   /* Calculate the centroid of the matrix                              */
   blFindCentroid(matrix, numX, numY, centroid);

   /* Calculate the covariance matrix                                   */
   for(i=0; i<numY; i++)
   {
      for(j=i; j<numY; j++)
      {
         total=0.0;

         for(k=0; k<numX; k++)
         {
            total += ( (matrix[k][i] - centroid[i]) * 
                       (matrix[k][j] - centroid[j]) );
         }

         /* cov[i][j] = total/(REAL)numX;                               */
         /* cov[j][i] = total/(REAL)numX;                               */

         cov[i][j] = total;
         cov[j][i] = total;
      }
   }

   /* Free memory and return                                            */
   free(centroid);

   return(TRUE);
}



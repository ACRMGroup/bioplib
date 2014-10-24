#ifndef _REGRESSION_H
#define _REGRESSION_H

BOOL blCalculateBestFitLine(double **coordinates,
                            int numberOfPoints,
                            int numberOfDimensions,
                            double *centroid,
                            double *eigenVector);
void blFindCentroid(REAL **coordinates, int numberOfPoints, 
                    int numberOfDimensions, REAL *centroid);
BOOL blCalculateCovarianceMatrix(REAL **x, int numX,
                                 int numY, REAL **cov);

/************************************************************************/
/* Include deprecated functions                                         */
#define _REGRESSION_H_DEPRECATED
#include "deprecated.h"
/************************************************************************/

#endif

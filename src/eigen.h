#ifndef _EIGEN_H
#define _EIGEN_H 1
int blEigen(REAL **M, REAL **Vectors, REAL *lambda, int n);
#define EIGEN_NOMEMORY   (-1)
#define EIGEN_NOCONVERGE (-2)

/************************************************************************/
/* Include deprecated functions                                         */
#define _EIGEN_H_DEPRECATED
#include "deprecated.h"
/************************************************************************/


#endif

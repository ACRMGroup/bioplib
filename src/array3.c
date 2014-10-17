/************************************************************************/
/**

   \file       array3.c
   
   \version    V1.1
   \date       07.07.14
   \brief      Allocate and free 3D arrays
   
   \copyright  (c) Dr. Andrew C. R. Martin, University of Reading, 2002-14
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

   Creates a 3D array where the first dimension is a set of pointers. This
   is better for passing into subroutines than the conventional C method
   of simply declaring:

      TYPE  matrix[10][10][10];

   which, when passed to a fuunction, loses the concept of dimensions
   unless the matrix is explicitly defined with these dimension in the
   function.
   
   This routine creates an array of pointers to 1-D arrays and can thus be
   passed to functions successfully.

**************************************************************************

   Usage:
   ======
   matrix = (TYPE ***)Array3D(sizeof(TYPE), nrows, ncolumns, nplanes);
   
\code
   matrix = (float **)Array3D(sizeof(float), 10, 10, 10);
\endcode

   Returns NULL (having freed any allocated memory) if there is a problem.

**************************************************************************

   Revision History:
   =================
-  V1.0  30.05.02 Original
-  V1.1  07.07.14 Include array.h. Remove FreeArray3D() prototype. 
         Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP Array handling

   #ROUTINE blArray3D()
   Create a 3D array of elements of size `size' with dimensions `dim1' 
   rows by `dim2' columns by `dim3' planes

   #ROUTINE blFreeArray3D()
   Frees a 3D array with dimensions `dim1' rows by `dim2' columns by
   `dim3' planes.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "array.h"

/************************************************************************/
/* Defines and macros
*/
#ifndef NULL
#define NULL ((void *)0)
#endif


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>char ***blArray3D(int size, int dim1, int dim2, int dim3)
   ---------------------------------------------------------
*//**

   \param[in]     size    Size of an array element
   \param[in]     dim1    First dimension (number of rows)
   \param[in]     dim2    Second dimension (number of columns)
   \param[in]     dim3    Third dimension (number of planes)
   \return                Array of pointers. Must be cast to required 
                          type

   Create a 3D array of elements of size `size' with dimensions `dim1' 
   rows by `dim2' columns by `dim3' planes

-  30.05.02 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
char ***blArray3D(int size, int dim1, int dim2, int dim3)
{
   char  ***array  = NULL;
   int   i, j;
   
   /* Allocate memory for the outer dimension array                     */
   if((array = (char ***)malloc(dim1 * sizeof(char **))) == NULL)
      return(NULL);
      
   /* Set all positions to NULL                                         */
   for(i=0; i<dim1; i++)   array[i] = NULL;

   /* Allocate memory for each array in the second dimension            */
   for(i=0; i<dim1; i++)
   {
      /* If allocation fails, jump to badexit                           */
      if((array[i] = (char **)malloc(dim2 * sizeof(char *))) == NULL)
         goto badexit;

      /* Set all positions to NULL                                      */
      for(j=0; j<dim2; j++)   array[i][j] = NULL;

      /* Allocate memory for each array in the third dimension          */
      for(j=0; j<dim2; j++)
      {
         /* If allocation fails, jump to badexit                        */
         if((array[i][j] = (char *)malloc(dim3 * size)) == NULL)
            goto badexit;
      }
   }
   
   return(array);
   
badexit:
   blFreeArray3D(array, dim1, dim2, dim3);

   return(NULL);
}

/************************************************************************/
/*>void blFreeArray3D(char ***array, int dim1, int dim2, int dim3)
   ---------------------------------------------------------------
*//**

   \param[in]     array Array of pointers to be freed
   \param[in]     dim1  First dimension (number of rows)
   \param[in]     dim2  Second dimension (number of columns)
   \param[in]     dim3  Third dimension (number of planes)

   Frees a 3D array with dimensions `dim1' rows by `dim2' columns by
   `dim3' planes.

-  30.05.02 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blFreeArray3D(char ***array, int dim1, int dim2, int dim3)
{
   int   i, j;
   
   if(array!=NULL)
   {
      for(i=0; i<dim1; i++)
      {
         if(array[i] != NULL)
         {
            for(j=0; j<dim2; j++)
            {
               if(array[i][j] != NULL)
                  free(array[i][j]);
            }
            free(array[i]);
         }
      }
      
      free(array);
   }
}

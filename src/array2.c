/************************************************************************/
/**

   \file       array2.c
   
   \version    V1.5
   \date       07/07.14
   \brief      Allocate and free 2D arrays
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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

   Creates a 2D array where the first dimension is a set of pointers. This
   is better for passing into subroutines than the conventional C method
   of simply declaring:

      TYPE  matrix[10][10];

   which, when passed to a function, loses the concept of dimensions
   unless the matrix is explicitly defined with these dimension in the
   function.
   
   This routine creates an array of pointers to 1-D arrays and can thus be
   passed to functions successfully.

**************************************************************************

   Usage:
   ======
   matrix = (TYPE **)Array2D(sizeof(TYPE), nrows, ncolumns);
   
\code
   matrix = (float **)Array2D(sizeof(float), 10, 10);
\endcode

   Returns NULL (having freed any allocated memory) if there is a problem.

**************************************************************************

   Revision History:
   =================
-  V1.0  07.10.92 Original
-  V1.1  29.01.93 Added includes of sysdefs.h & malloc.h for MS-DOS
-  V1.2  16.06.93 Includes stdlib.h rather than malloc.h
-  V1.3  01.03.94 Corrected other include file usage
-  V1.4  18.03.94 Added NULL definition for systems which don't define
                  it in stdlib.h
-  V1.5  07.07.14 Include array.h Use bl prefix for functions By: CTP


*************************************************************************/
/* Includes
*/
#include <stdlib.h>

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
/*>char **blArray2D(int size, int dim1, int dim2)
   --------------------------------------------
*//**

   \param[in]     size    Size of an array element
   \param[in]     dim1    First dimension (number of rows)
   \param[in]     dim2    Second dimension (number of columns)
   \return                Array of pointers. Must be cast to required 
                          type

   Create a 2D array of elements of size `size' with dimensions `dim1' 
   rows by `dim2' columns.

-  07.10.92 Original
-  12.07.93 Tidied and commented
-  07.07.14 Use bl prefix for functions By: CTP

*/
char **blArray2D(int size, 
                 int dim1, 
                 int dim2)
{
   char  **array  = NULL;
   int   i;
   
   /* Allocate memory for the outer dimension array                     */
   if((array = (char **)malloc(dim1 * sizeof(char *))) == NULL)
      return(NULL);
      
   /* Set all positions to NULL                                         */
   for(i=0; i<dim1; i++)   array[i] = NULL;

   /* Allocate memory for each array in the second dimension            */
   for(i=0; i<dim1; i++)
   {
      /* If allocation fails, jump to badexit                           */
      if((array[i] = (char *)malloc(dim2 * size)) == NULL)
         goto badexit;
   }
   
   return(array);
   
badexit:
   for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
   free(array);
   return(NULL);
}

/************************************************************************/
/*>void blFreeArray2D(char **array, int dim1, int dim2)
   ----------------------------------------------------
*//**

   \param[in]     array Array of pointers to be freed
   \param[in]     dim1  First dimension (number of rows)
   \param[in]     dim2  Second dimension (number of columns)

   Frees a 2D array with dimensions `dim1' rows by `dim2' columns.

-  07.10.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blFreeArray2D(char   **array,
                   int    dim1, 
                   int    dim2)
{
   int   i;
   
   if(array)
   {
      for(i=0; i<dim1; i++)   if(array[i]) free(array[i]);
      free(array);
   }
}

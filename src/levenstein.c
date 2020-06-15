/************************************************************************/
/**

   Program:    
   \file       levenstein.c
   
   \version    V1.0
   \date       15.06.20
   \brief      Calculate a Levenshtein distance
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2020
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
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

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Fast and low-memory method to calculate the Levenshtein distance
   between two strings. Maintains only two rows of the matrix rather
   than the whole thing.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

/************************************************************************/
/* Defines and macros
*/
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MIN3
#define MIN3(c,d,e) (MIN(MIN((c),(d)),(e)))
#endif

#define SWAPINTPTR(a, b)  {                     \
   int *t = a;                                  \
   a = b;                                       \
   b = t;                                       \
   }
   

/************************************************************************/
/*>int blLevenshteinDistance(char *columnString, char *rowString)
   --------------------------------------------------------------
*//**
   \input   columnString   A NULL-terminated string pointer   
   \input   rowString      A NULL-terminated string pointer
   \return                 Levenshtein distance (<0 for error)

   Calculates the Levenshtein edit distance using the iterative
   method with two matrix rows. See 
   https://en.m.wikipedia.org/wiki/Levenshtein_distance

   The method is memory efficient as it uses only two rows from the
   complete matrix - to minimize memory usage, ensure that rowString
   is the shorter string.

   However, it is suboptimal. The amount of memory required may be
   reduced to one row and one (index) word of overhead, for better
   cache locality. Thus speed comes from using very little memory,
   often keeping the buffer entirely in cache, and reducing the amount
   of work by skipping any prefix and postfix that won't add to the
   cost. You only really need to know the three values when you want
   to update a cell in the matrix and you can keep two of them in a
   buffer, while keeping the third value in a fixed location.

   See
   https://web.archive.org/web/20180612143641/https://bitbucket.org/clearer/iosifovich/

   and the C++ implementation of that version at
   https://github.com/rljacobson/Levenshtein/blob/master/damlev.cpp

-  15.06.20  Original   By: ACRM
*/
int blLevenshteinDistance(char *columnString, char *rowString)
{
   int columnSize   = strlen(columnString);
   int rowSize      = strlen(rowString);
   int *rowPrevious = NULL;
   int *rowCurrent  = NULL;
   int i, j, result;
   
   /* Allocate memory for the two rows                                  */
   if((rowPrevious = (int *)malloc((rowSize+1) * sizeof(int)))==NULL)
      return(-1);
   
   if((rowCurrent  = (int *)malloc((rowSize+1) * sizeof(int)))==NULL)
   {
      free(rowPrevious);
      return(-1);
   }
   
   /* Initialize rowPrevious (the previous row of distances).  This
      row is A[0][i] (where A[][] would be the full matrix). The edit
      distance for an empty columnString is just the number of
      characters to delete from rowString
   */
   for(i=0; i<=rowSize; i++)
   {
      rowPrevious[i] = i;
   }
   
   for(i=0; i<columnSize; i++)
   {
      /* Calculate rowCurrent (current row distances) from the
         previous row (rowPrevious). The first element of rowCurrent
         is A[i+1][0], the edit distance comes from deleting (i+1)
         chars from columnString to match an empty rowString
      */
      rowCurrent[0] = i + 1;

      /* use the formula to fill in the rest of the row                 */
      for(j=0; j<rowSize; j++)
      {
         /* Calculate costs for A[i+1][j+1]                             */
         int deletionCost     = rowPrevious[j+1] + 1;
         int insertionCost    = rowCurrent[j] + 1;
         int substitutionCost = 0;
         
         if(columnString[i] == rowString[j])
         {
            substitutionCost = rowPrevious[j];
         }
         else
         {
            substitutionCost = rowPrevious[j] + 1;
         }
         
         rowCurrent[j+1] = MIN3(deletionCost,
                                insertionCost,
                                substitutionCost);
      }

      /* Copy rowCurrent (current row) to rowPrevious (previous row)
         for the next iteration. Since the data in rowCurrent are
         always invalidated, for efficiency, the swap is done just by
         swapping the pointers.
      */
      SWAPINTPTR(rowPrevious, rowCurrent);
   }

   /* After the final swap, the results of rowCurrent are now in 
      rowPrevious 
   */ 
   result = rowPrevious[rowSize];
   free(rowPrevious);
   free(rowCurrent);
   
   return(result);
}


/************************************************************************/
#ifdef TEST
#include <stdio.h>
#include "stredit.h"
int main(int argc, char **argv)
{
   static char *one = "kitten";
   static char *two = "sitting";

   int dist = blLevenshteinDistance(one, two);
   printf("Distance: %d\n", dist);

   return(0);
}
#endif



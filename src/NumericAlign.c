/************************************************************************/
/**

   \file       NumericAlign.c
   
   \version    V1.3
   \date       07.07.14
   \brief      Perform Needleman & Wunsch sequence alignment on two
               sequences encoded as numeric symbols.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin / University of Reading 1993-2014
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

   Note, the code herein is very heavily based on code written by Dr. 
   Andrew C.R. Martin while self-employed. Some modifications were made to
   that original code while employed at University College London. This
   version which handles sequences encoded as arrays of numbers rather
   than as character arrays was modified from the original version(s)
   while employed at Reading University.

   Description:
   ============

   A simple Needleman & Wunsch Dynamic Programming alignment of 2 
   sequences encoded as numeric symbols.  
   A window is not used so the routine may be a bit slow on long 
   sequences.

**************************************************************************

   Usage:
   ======

   First call NumericReadMDM() to read the mutation data matrix, then call
   NumericAffineAlign() to align the sequences.

**************************************************************************

   Revision History:
   =================
-  V1.0  08.03.00 A modified version of align.c V2.12 written by ACRM
                  from 19.06.90 to 06.03.00
-  V1.1  28.09.00 Fixed bug at end of alignment if one sequence finishes
                  first
-  V1.2  06.02.03 Fixed for new version of GetWord()
-  V1.3  07.07.14 Use bl prefix for functions By: CTP


*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP Alignment
   #ROUTINE  blNumericReadMDM()
   Read mutation data matrix into static global arrays for number-encoded
   sequences

   #ROUTINE  blNumericCalcMDMScore()
   Calculate score from static globally stored mutation data matrix for
   number-encoded sequences

   #ROUTINE  blNumericAffineAlign()
   Perform simple N&W alignment using sequences encodede as arrays of
   numeric tokens
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "SysDefs.h"
#include "macros.h"
#include "array.h"
#include "general.h"
#include "seq.h"

/************************************************************************/
/* Defines and macros
*/
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define DATAENV "DATADIR"   /* Environment variable or assign           */

#define MAXBUFF 2048

/* Type definition to store a X,Y coordinate pair in the matrix         */
typedef struct
{
   int x, y;
}  XY;


/************************************************************************/
/* Globals
*/
static int  **sMDMScore;
static int  sMDMSize = 0;

/************************************************************************/
/*
BOOL NumericReadMDM(char *mdmfile);
int NumericCalcMDMScore(int resa, int resb);
int NumericAffineAlign(int  *seq1, 
                       int  length1, 
                       int  *seq2, 
                       int  length2, 
                       BOOL verbose, 
                       BOOL identity, 
                       int  penalty, 
                       int  penext,
                       int *align1, 
                       int *align2,
                       int  *align_len);
*/

/************************************************************************/
/* Prototypes
*/
static int  NumericSearchForBest(int **matrix, int length1, int length2, 
                          int *BestI, int *BestJ, int *seq1, int *seq2, 
                          int *align1, int *align2);
static int  NumericTraceBack(int **matrix, XY **dirn, int length1, 
                             int length2, int *seq1, int *seq2, 
                             int *align1, int *align2, 
                             int *align_len);



/************************************************************************/
/*>static int NumericSearchForBest(int **matrix, int length1, 
                                   int length2, int *BestI, int *BestJ, 
                                   int *seq1, int *seq2, int *align1, 
                                   int *align2)
   ----------------------------------------------------------------
*//**
   \param[in]     **matrix   N&W matrix
   \param[in]     length1    Length of first sequence
   \param[in]     length2    Length of second sequence
   \param[in]     *BestI     x position of highest score
   \param[in]     *BestJ     y position of highest score
   \param[in]     *seq1      First sequence
   \param[in]     *seq2      Second sequence
   \param[out]    *align1    First sequence with end aligned correctly
   \param[out]    *align2    Second sequence with end aligned correctly
   \return                     Alignment length thus far

   Searches the outside of the matrix for the best score and starts the
   alignment by putting in any starting 0s as insertion symbols.

   Identical to align.c/SearchForBest(), but uses int arrays rather than
   characters

-  08.03.00 Original based on align.c/SearchForBest() 08.10.92 By: ACRM
*/
static int NumericSearchForBest(int  **matrix, 
                                int  length1, 
                                int  length2, 
                                int  *BestI, 
                                int  *BestJ,
                                int  *seq1, 
                                int  *seq2, 
                                int  *align1, 
                                int  *align2)
{
   int   ai, 
         besti,   bestj, 
         i,       j;
   
   /* Now search the outside of the matrix for the highest scoring cell */
   ai    = 0;
   besti = 0;
   for(i = 1; i < length1; i++) 
   {
      if(matrix[i][0] > matrix[besti][0]) besti = i;
   }
   bestj = 0;
   for(j = 1; j < length2; j++)
   {
      if(matrix[0][j] > matrix[0][bestj]) bestj = j;
   }
   if(matrix[besti][0] > matrix[0][bestj])
   {
      *BestI = besti;
      *BestJ = 0;
      for(i=0; i<*BestI; i++)
      {
         align1[ai] = seq1[i];
         align2[ai++] = 0;
      }
   }
   else
   {
      *BestI = 0;
      *BestJ = bestj;
      for(j=0; j<*BestJ; j++)
      {
         align1[ai] = 0;
         align2[ai++] = seq2[j];
      }
   }
   return(ai);
}




/************************************************************************/
/*>BOOL blNumericReadMDM(char *mdmfile)
   ------------------------------------
*//**

   \param[in]     *mdmfile    Mutation data matrix filename
   \return                      Success?
   
   Read mutation data matrix into static global arrays. The matrix may
   have comments at the start introduced with a ! in the first column.
   The matrix must be complete (i.e. a triangular matrix will not
   work). A line describing the residue types must appear, and may
   be placed before or after the matrix itself

   Identical to align.c/ReadMDM() but reads into a different static
   2D array and doesn't read a symbol identifier line from the file
   as the symbols are numeric and always start from 1 (0 is used as
   the insert character)

-  08.03.00 Original based on align.c/ReadMDM() 26.07.95 By: ACRM
-  06.02.03 Fixed for new version of GetWord()
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blNumericReadMDM(char *mdmfile)
{
   FILE *mdm = NULL;
   int  i, j;
   char buffer[MAXBUFF],
        word[16],
        *p;
   BOOL noenv;

   if((mdm=blOpenFile(mdmfile, DATAENV, "r", &noenv))==NULL)
   {
      return(FALSE);
   }

   /* Read over any comment lines                                       */
   while(fgets(buffer,MAXBUFF,mdm))
   {
      TERMINATE(buffer);
      KILLLEADSPACES(p,buffer);
      if(strlen(p) && p[0] != '!')
         break;
   }

   /* See how many fields there are in the buffer                       */
   for(p = buffer, sMDMSize = 0; p!=NULL; sMDMSize++)
      p = blGetWord(p, word, 16);


   /* Allocate memory for the MDM and the AA List                       */
   if((sMDMScore = (int **)blArray2D(sizeof(int),sMDMSize,sMDMSize))==NULL)
      return(FALSE);

   /* Fill the matrix with zeros                                        */
   for(i=0; i<sMDMSize; i++)
   {
      for(j=0; j<sMDMSize; j++)
      {
         sMDMScore[i][j] = 0;
      }
   }
   
   i=0;
   do
   {
      TERMINATE(buffer);
      KILLLEADSPACES(p, buffer);
      if(strlen(p))
      {
         blGetWord(buffer, word, 16);
         if(sscanf(word,"%d",&j))    /* A row of numbers                */
         {
            for(p = buffer, j = 0; p!=NULL && j<sMDMSize; j++)
            {
               p = blGetWord(p, word, 16);
               sscanf(word,"%d",&(sMDMScore[i][j]));
            }
            i++;
         }
      }
   } while(fgets(buffer,MAXBUFF,mdm));
   
   fclose(mdm);
   
   return(TRUE);
}

/************************************************************************/
/*>int blNumericCalcMDMScore(int resa, int resb)
   ---------------------------------------------
*//**

   \param[in]     resa      First token  
   \param[in]     resb      Second token  
   \return                  score

   Calculate score from static globally stored mutation data matrix

   Identical to align.c/CalcMDMScore(), but uses a different static score
   array and takes integer parameters. These are used as direct lookups
   into the score array rather than being searched.

-  08.03.00 Original based on align.c/CalcMDMScore() 11.07.96 By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blNumericCalcMDMScore(int resa, int resb)
{
   int        i,j;
   static int NWarn = 0;
   BOOL       Warned = FALSE;

   i = resa-1;
   j = resb-1;
   
   if(i>=sMDMSize) 
   {
      if(NWarn < 10)
         printf("Token %d not found in matrix\n",resa);
      else if(NWarn == 10)
         printf("More token not found in matrix...\n");
      Warned = TRUE;
   }
   if(j>=sMDMSize) 
   {
      if(NWarn < 10)
         printf("Token %d not found in matrix\n",resb);
      else if(NWarn == 10)
         printf("More tokens not found in matrix...\n");
      Warned = TRUE;
   }
   
   if(Warned)
   { 
      NWarn++;
      return(0);
   }

   return(sMDMScore[i][j]);
}                               

/************************************************************************/
/*>int blNumericAffineAlign(int *seq1, int length1, int *seq2, int length2, 
                            BOOL verbose, BOOL identity, int penalty, 
                            int penext, int *align1, int *align2, 
                            int *align_len)
   ---------------------------------------------------------------------
*//**
   \param[in]     *seq1         First sequence of tokens
   \param[in]     length1       First sequence length
   \param[in]     *seq2         Second sequence of tokens
   \param[in]     length2       Second sequence length
   \param[in]     verbose       Display N&W matrix
   \param[in]     identity      Use identity matrix
   \param[in]     penalty       Gap insertion penalty value
   \param[in]     penext        Extension penalty
   \param[out]    *align1       Sequence 1 aligned
   \param[out]    *align2       Sequence 2 aligned
   \param[out]    *align_len    Alignment length
   \return                         Alignment score (0 on error)
            
   Perform simple N&W alignment of seq1 and seq2. No window is used, so
   will be slow for long sequences.

   The sequences come as integer arrays containing numeric tokens

   Note that you must allocate sufficient memory for the aligned 
   sequences.
   The easy way to do this is to ensure that align1 and align2 are
   of length (length1+length2).

   Identical to align.c/affinealign(), but uses integer arrays

-  08.03.00 Original based on align.c/affinealign() 06.03.00 By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blNumericAffineAlign(int  *seq1, 
                         int  length1, 
                         int  *seq2, 
                         int  length2, 
                         BOOL verbose, 
                         BOOL identity, 
                         int  penalty, 
                         int  penext,
                         int  *align1, 
                         int  *align2,
                         int  *align_len)
{
   XY    **dirn   = NULL;
   int   **matrix = NULL,
         maxdim,
         i,    j,    k,    l,
         i1,   j1,
         dia,  right, down,
         rcell, dcell, maxoff,
         match = 1,
         thisscore,
         gapext,
         score;
   
   maxdim = MAX(length1, length2);
   
   /* Initialise the score matrix                                       */
   if((matrix = (int **)blArray2D(sizeof(int), maxdim, maxdim))==NULL)
      return(0);
   if((dirn   = (XY **)blArray2D(sizeof(XY), maxdim, maxdim))==NULL)
      return(0);
      
   for(i=0;i<maxdim;i++)
   {
      for(j=0;j<maxdim;j++)
      {
         matrix[i][j] = 0;
         dirn[i][j].x = -1;
         dirn[i][j].y = -1;
      }
   }
    
   /* Fill in scores up the right hand side of the matrix               */
   for(j=0; j<length2; j++)
   {
      if(identity)
      {
         if(seq1[length1-1] == seq2[j]) matrix[length1-1][j] = match;
      }
      else
      {
         matrix[length1-1][j] = blNumericCalcMDMScore(seq1[length1-1], 
                                                      seq2[j]);
      }
   }

   /* Fill in scores along the bottom row of the matrix                 */
   for(i=0; i<length1; i++)
   {
      if(identity)
      {
         if(seq1[i] == seq2[length2-1]) matrix[i][length2-1] = match;
      }
      else
      {
         matrix[i][length2-1] = blNumericCalcMDMScore(seq1[i], 
                                                      seq2[length2-1]);
      }
   }

   i = length1 - 1;
   j = length2 - 1;
   
   /* Move back along the diagonal                                      */
   while(i > 0 && j > 0)
   {
      i--;
      j--;

      /* Fill in the scores along this row                              */
      for(i1 = i; i1 > -1; i1--)
      {
         dia   = matrix[i1+1][j+1];

         /* Find highest score to right of diagonal                     */
         rcell = i1+2;
         if(i1+2 >= length1)  right = 0;
         else                 right = matrix[i1+2][j+1] - penalty;
         
         gapext = 1;
         for(k = i1+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j+2;
         if(j+2 >= length2)  down = 0;
         else                down   = matrix[i1+1][j+2] - penalty;
         
         gapext = 1;
         for(l = j+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i1+1][l] - (penalty + gapext*penext);

            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }
         
         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i1][j] = dia;
            dirn[i1][j].x = i1+1;
            dirn[i1][j].y = j+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i1][j] = right;
               dirn[i1][j].x = rcell;
               dirn[i1][j].y = j+1;
            }
            else
            {
               matrix[i1][j] = down;
               dirn[i1][j].x = i1+1;
               dirn[i1][j].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i1] == seq2[j]) matrix[i1][j] += match;
         }
         else
         {
            matrix[i1][j] += blNumericCalcMDMScore(seq1[i1],seq2[j]);
         }
      }

      /* Fill in the scores in this column                              */
      for(j1 = j; j1 > -1; j1--)
      {
         dia   = matrix[i+1][j1+1];
         
         /* Find highest score to right of diagonal                     */
         rcell = i+2;
         if(i+2 >= length1)   right = 0;
         else                 right = matrix[i+2][j1+1] - penalty;

         gapext = 1;
         for(k = i+3; k<length1; k++, gapext++)
         {
            thisscore = matrix[k][j1+1] - (penalty + gapext*penext);
            
            if(thisscore > right) 
            {
               right = thisscore;
               rcell = k;
            }
         }

         /* Find highest score below diagonal                           */
         dcell = j1+2;
         if(j1+2 >= length2)  down = 0;
         else                 down = matrix[i+1][j1+2] - penalty;

         gapext = 1;
         for(l = j1+3; l<length2; l++, gapext++)
         {
            thisscore = matrix[i+1][l] - (penalty + gapext*penext);
            
            if(thisscore > down) 
            {
               down = thisscore;
               dcell = l;
            }
         }

         /* Set score to best of these                                  */
         maxoff = MAX(right, down);
         if(dia >= maxoff)
         {
            matrix[i][j1] = dia;
            dirn[i][j1].x = i+1;
            dirn[i][j1].y = j1+1;
         }
         else
         {
            if(right > down)
            {
               matrix[i][j1] = right;
               dirn[i][j1].x = rcell;
               dirn[i][j1].y = j1+1;
            }
            else
            {
               matrix[i][j1] = down;
               dirn[i][j1].x = i+1;
               dirn[i][j1].y = dcell;
            }
         }
       
         /* Add the score for a match                                   */
         if(identity)
         {
            if(seq1[i] == seq2[j1]) matrix[i][j1] += match;
         }
         else
         {
            matrix[i][j1] += blNumericCalcMDMScore(seq1[i],seq2[j1]);
         }
      }
   } 
   
   score = NumericTraceBack(matrix, dirn, length1, length2,
                            seq1, seq2, align1, align2, align_len);

   if(verbose)
   {
      printf("Matrix:\n-------\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("%3d ",matrix[i][j]);
         }
         printf("\n");
      }

      printf("Path:\n-----\n");
      for(j=0; j<length2;j++)
      {
         for(i=0; i<length1; i++)
         {
            printf("(%3d,%3d) ",dirn[i][j].x,dirn[i][j].y);
         }
         printf("\n");
      }
   }
    
   blFreeArray2D((char **)matrix, maxdim, maxdim);
   blFreeArray2D((char **)dirn,   maxdim, maxdim);
    
   return(score);
}

/************************************************************************/
/*>static int NumericTraceBack(int **matrix, XY **dirn, 
                               int length1, int length2, 
                               int *seq1, int *seq2, int *align1, 
                               int *align2, int *align_len)
   ----------------------------------------------------------------
*//**
   \param[in]     **matrix   N&W matrix
   \param[in]     **dirn     Direction Matrix
   \param[in]     length1    Length of first sequence
   \param[in]     length2    Length of second sequence
   \param[in]     *seq1      First sequence
   \param[in]     *seq2      Second sequence
   \param[out]    *align1    First sequence aligned
   \param[out]    *align2    Second sequence aligned
   \param[out]    *align_len Aligned sequence length
   \return                     Alignment score

   Does the traceback to find the aligment.

   Identical to align.c/TraceBack(), but uses integer arrays.

-  08.03.00 Original based on align.c/TraceBack() 06.03.00 By: ACRM
-  28.09.00 Fixed bug at end of alignment if one sequence finishes
            first
*/
static int NumericTraceBack(int  **matrix, 
                            XY   **dirn,
                            int  length1, 
                            int  length2, 
                            int *seq1, 
                            int *seq2, 
                            int *align1, 
                            int *align2, 
                            int  *align_len)
{
   int   i,    j, 
         ai, 
         BestI,BestJ;
   XY    nextCell;

   ai = NumericSearchForBest(matrix, length1, length2, &BestI, &BestJ, 
                             seq1, seq2, align1, align2);

   /* Now trace back to find the alignment                              */
   i            = BestI;
   j            = BestJ;
   align1[ai]   = seq1[i];
   align2[ai++] = seq2[j];

   while(i < length1-1 && j < length2-1)
   {
      nextCell.x = dirn[i][j].x;
      nextCell.y = dirn[i][j].y;
      if((nextCell.x == i+1) && (nextCell.y == j+1))
      {
         /* We are inheriting from the diagonal                         */
         i++;
         j++;
      }
      else if(nextCell.y == j+1)
      {
         /* We are inheriting from the off-diagonal inserting a gap in
            the y-sequence (seq2)
         */
         i++;
         j++;
         while((i < nextCell.x) && (i < length1-1))
         {
            align1[ai] = seq1[i++];
            align2[ai++] = 0;
         }
      }
      else if(nextCell.x == i+1)
      {
         /* We are inheriting from the off-diagonal inserting a gap in
            the x-sequence (seq1)
         */
         i++;
         j++;
         while((j < nextCell.y) && (j < length2-1))
         {
            align1[ai] = 0;
            align2[ai++] = seq2[j++];
         }
      }
      else
      {
         /* Cockup!                                                     */
         fprintf(stderr,"align.c/TraceBack() internal error\n");
      }
      
      align1[ai]   = seq1[i];
      align2[ai++] = seq2[j];
   }

   /* If one sequence finished first, fill in the end with insertions   */
   if(i < length1-1)
   {
      for(j=i+1; j<length1; j++)
      {
         align1[ai]   = seq1[j];
         align2[ai++] = 0;
      }
   }
   else if(j < length2-1)
   {
      for(i=j+1; i<length2; i++)
      {
         align1[ai]   = 0;
         align2[ai++] = seq2[i];
      }
   }
   
   *align_len = ai;
   
   return(matrix[BestI][BestJ]);
}




            
      
#ifdef DEMO   
int main(int argc, char **argv)
{
   int  seq1[] = {1, 3, 1, 3, 7, 9, 5, 6},
        seq2[] = {1, 3, 1, 3, 5, 6},
        align1[100],
        align2[100];
   int  score, al_len, i;
   

   NumericReadMDM("numtopmat.mat"); 
   
   score = NumericAffineAlign(seq1, 8, seq2, 6, 
                              TRUE, FALSE,
                              5, 0, align1, align2, &al_len);

   align1[al_len] = '\0';
   align2[al_len] = '\0';

   for(i=0;i<al_len;i++)
      printf("%2d=", align1[i]);
   printf("\n");
   
   for(i=0;i<al_len;i++)
      printf("%2d=", align2[i]);
   printf("\n");
   
   return(0);
}
#endif


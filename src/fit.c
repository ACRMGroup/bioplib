/************************************************************************/
/**

   \file       fit.c
   
   \version    V1.7
   \date       17.07.14
   \brief      Perform least squares fitting of coordinate sets
   
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

   This code performs least squares fitting of two coordinate set using
   the method of McLachlan as modified by Sutcliffe.

**************************************************************************

   Usage:
   ======
   Passed two coordinate arrays both centred around the origin and,
   optionally, an array of weights, returns a rotation matrix.

**************************************************************************

   Revision History:
   =================
-  V1.0  04.02.91 Original
-  V1.1  01.06.92 ANSIed and static'd
-  V1.2  08.12.92 Changed abs() to ABS() using macros.h. Includes stdio.h
-  V1.3  11.02.94 Changed column flag to BOOL
-  V1.4  03.06.97 Corrected documentation
-  V1.5  03.04.09 Initialize clep in qikfit() By: CTP
-  V1.6  07.07.14 Use bl prefix for functions By: CTP
-  V1.7  17.07.14 Removed unused varables  By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Coordinate Fitting
   #SUBGROUP Fitting based on coordinate arrays
   #ROUTINE  blMatfit()
   Fit coordinate array x2 to x1 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if wt1 is not NULL.
   If column is set the matrix will be returned column-wise rather 
   than row-wise.
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdio.h>

#include "MathType.h"
#include "fit.h"
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define SMALL  1.0e-20     /* Convergence cutoffs                       */
#define SMALSN 1.0e-10

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column);

/************************************************************************/
/*>BOOL blMatfit(COOR *x1, COOR *x2, REAL rm[3][3], int n,
                 REAL *wt1, BOOL column)
   -------------------------------------------------------
*//**

   \param[in]     *x1         First (fixed) array of coordinates
   \param[in]     *x2         Second (mobile) array of coordinates
   \param[in]     n           Number of coordinates
   \param[in]     *wt1        Weight array or NULL
   \param[in]     column      TRUE: Output a column-wise matrix (as used
                                 by FRODO)
                              FALSE: Output a standard row-wise matrix.
   \param[out]    rm          Returned rotation matrix
   \return                    TRUE:  success
                              FALSE: error

   Fit coordinate array x2 to x1 both centred around the origin and of 
   length n. Optionally weighted with the wt1 array if wt1 is not NULL.
   If column is set the matrix will be returned column-wise rather 
   than row-wise.

-  04.02.91 Original based on code by Mike Sutcliffe
-  01.06.92 ANSIed & doc'd
-  17.06.93 various changes for release (including parameters)
-  11.03.94 column changed to BOOL
-  25.11.02 Corrected header!
-  07.07.14 Use bl prefix for functions By: CTP

*/
BOOL blMatfit(COOR    *x1,        /* First coord array    */
              COOR    *x2,        /* Second coord array   */
              REAL    rm[3][3],   /* Rotation matrix      */
              int     n,          /* Number of points     */
              REAL    *wt1,       /* Weight array         */
              BOOL    column)     /* Column-wise output   */
{
   int  i,j;
   REAL umat[3][3];

   
   if(n<2)
   {
      return(FALSE);
   }

   if(wt1)
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += wt1[j] * x1[j].x * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].x * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += wt1[j] * x1[j].y * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].y * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += wt1[j] * x1[j].z * x2[j].x;
                  umat[i][1] += wt1[j] * x1[j].z * x2[j].y;
                  umat[i][2] += wt1[j] * x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   } 
   else
   {
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++) umat[i][j] = 0.0;

         for(j=0;j<n;j++)
         {
            switch(i)
            {
               case 0:
                  umat[i][0] += x1[j].x * x2[j].x;
                  umat[i][1] += x1[j].x * x2[j].y;
                  umat[i][2] += x1[j].x * x2[j].z;
                  break;
               case 1:
                  umat[i][0] += x1[j].y * x2[j].x;
                  umat[i][1] += x1[j].y * x2[j].y;
                  umat[i][2] += x1[j].y * x2[j].z;
                  break;
               case 2:
                  umat[i][0] += x1[j].z * x2[j].x;
                  umat[i][1] += x1[j].z * x2[j].y;
                  umat[i][2] += x1[j].z * x2[j].z;
                  break;
            }
         }
      }
   }
   qikfit(umat,rm,column);

   return(TRUE);
}
   
/************************************************************************/
/*>static void qikfit(REAL umat[3][3], REAL rm[3][3], BOOL column)
   ---------------------------------------------------------------
*//**

   \param[in]     umat           The U matrix
   \param[in]     column         TRUE: Create a column-wise matrix
                                 (other way round from normal).
   \param[out]    rm             The output rotation matrix
  
   Does the actual fitting for matfit().
-  04.02.91 Original based on code by Mike Sutcliffe
-  01.06.92 ANSIed & doc'd
-  11.03.94 column changed to BOOL
-  03.04.09 Initialize clep for fussy compliers. By: CTP
-  17.07.14 Removed unused variables  By: ACRM
*/
static void qikfit(REAL  umat[3][3],
                   REAL  rm[3][3],
                   BOOL  column)
{
   
   REAL  rot[3][3],
         turmat[3][3],
         c[3][3],
         coup[3],
         dir[3],
         step[3],
         v[3],
         rtsum,rtsump,
         stp,stcoup,
         ud,tr,ta,cs,sn,ac,
         delta,
         gfac,
         cle,
         clep = 0.0;
   int   i,j,k,l,m,
         jmax,
         ncyc,
         nsteep,
         nrem;

   /* Rotate repeatedly to reduce couple about initial direction to zero.
      Clear the rotation matrix
   */
   for(l=0;l<3;l++)
   {
      for(m=0;m<3;m++)
         rot[l][m] = 0.0;
      rot[l][l] = 1.0;
   }

   /* Copy vmat[][] (sp) into umat[][] (dp)                             */
   jmax = 30;
   rtsum = umat[0][0] + umat[1][1] + umat[2][2];
   delta = 0.0;

   for(ncyc=0;ncyc<jmax;ncyc++)
   {
      /* Modified CG. For first and every NSTEEP cycles, set previous
         step as zero and do an SD step
      */
      nsteep = 3;
      nrem = ncyc-nsteep*(int)(ncyc/nsteep);

      if(!nrem)
      {
         for(i=0;i<3;i++) step[i]=0.0;
         clep = 1.0;
      }
      
      /* Couple                                                         */
      coup[0] = umat[1][2]-umat[2][1];
      coup[1] = umat[2][0]-umat[0][2];
      coup[2] = umat[0][1]-umat[1][0];
      cle     = sqrt(coup[0]*coup[0] + coup[1]*coup[1] + coup[2]*coup[2]);

      /* Gradient vector is now -coup                                   */
      gfac = (cle/clep)*(cle/clep);

      /* Value of rtsum from previous step                              */
      rtsump = rtsum;
      clep   = cle;
      if(cle < SMALL) break;

      /* Step vector conjugate to  previous                             */
      stp = 0.0;
      for(i=0;i<3;i++)
      {
         step[i]=coup[i]+step[i]*gfac;
         stp   += (step[i] * step[i]);
      }
      stp = 1.0/sqrt(stp);
         
      /* Normalised step                                                */
      for(i=0;i<3;i++) dir[i] = stp*step[i];

      /* Couple resolved along step direction                           */
      stcoup = coup[0]*dir[0] + coup[1]*dir[1] + coup[2]*dir[2];

      /* Component of UMAT along direction                              */
      ud = 0.0;
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            ud += umat[l][m]*dir[l]*dir[m];


      tr = umat[0][0]+umat[1][1]+umat[2][2]-ud;
      ta = sqrt(tr*tr + stcoup*stcoup);
      cs=tr/ta;
      sn=stcoup/ta;
         
      /* If cs<0 then posiiton is unstable, so don't stop               */
      if((cs>0.0) && (ABS(sn)<SMALSN)) break;
            
      /* Turn matrix for correcting rotation:

         Symmetric part
      */
      ac = 1.0-cs;
      for(l=0;l<3;l++)
      {
         v[l] = ac*dir[l];
         for(m=0;m<3;m++)
            turmat[l][m] = v[l]*dir[m];
         turmat[l][l] += cs;
         v[l]=dir[l]*sn;
      }

      /* Asymmetric part                                                */
      turmat[0][1] -= v[2];
      turmat[1][2] -= v[0];
      turmat[2][0] -= v[1];
      turmat[1][0] += v[2];
      turmat[2][1] += v[0];
      turmat[0][2] += v[1];

      /* Update total rotation matrix                                   */
      for(l=0;l<3;l++)
      {
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*rot[k][m];
         }
      }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            rot[l][m] = c[l][m];

      /* Update umat tensor                                             */
      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
         {
            c[l][m] = 0.0;
            for(k=0;k<3;k++)
               c[l][m] += turmat[l][k]*umat[k][m];
         }

      for(l=0;l<3;l++)
         for(m=0;m<3;m++)
            umat[l][m] = c[l][m];

      rtsum = umat[0][0] + umat[1][1] + umat[2][2];
      delta = rtsum - rtsump;

      /* If no improvement in this cycle then stop                      */
      if(ABS(delta)<SMALL) break;

      /* Next cycle                                                     */
   }

   /* Copy rotation matrix for output                                   */
   if(column)
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[j][i] = rot[i][j];
   }
   else
   {
      for(i=0;i<3;i++)
         for(j=0;j<3;j++)
            rm[i][j] = rot[i][j];
   }
}

/*************************************************************************

   Program:    ProFit
   File:       fitting.c
   
   Version:    V3.1
   Date:       31.03.09
   Function:   Protein Fitting program. 
   
   Copyright:  SciTech Software / UCL 1992-2009
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain.

   It may not be copied or made available to third parties, but may be
   freely used by non-profit-making organisations who have obtained it
   directly from the author or by FTP.

   You are requested to send EMail to the author to say that you are 
   using this code so that you may be informed of future updates.

   The code may not be made available on other FTP sites without express
   permission from the author.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If
   someone else breaks this code, the author doesn't want to be blamed
   for code that does not work! You may not distribute any
   modifications, but are encouraged to send them to the author so
   that they may be incorporated into future versions of the code.

   Such modifications become the property of Dr. Andrew C.R. Martin and
   SciTech Software though their origin will be acknowledged.

   The code may not be sold commercially or used for commercial purposes
   without prior permission from the author.
   
**************************************************************************

   Description:
   ============
   These routines perform the actual fitting and RMS calculation.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  25.09.92 Original
   V0.5  08.10.93 Various tidying for Unix & chaned for booklib 
   V0.6  05.01.94 More tidying
   V0.7  24.11.94 Skipped
   V0.8  17.07.95 Changed screen() calls to printf()
                  Multiple chains now work correctly.
   V1.0  18.07.95 Insert codes now work.
                  First official release (at last!).
   V1.1  20.07.95 Added WEIGHT command support and translation vector
                  output from MATRIX command
   V1.2  22.07.95 Skipped
   V1.3  31.07.95 Skipped
   V1.4  14.08.95 Fixed bug in CalcRMS() which was skipping the last 
                  residue when printing RMS by res.
   V1.5  21.08.95 Skipped
   V1.6  20.11.95 Skipped
   V1.6c 13.12.95 Extra info printed when zones mismatch
   V1.6e 31.05.96 Added test on B-values
   V1.6f 13.06.96 Added weight inverting for BWEIGHT command
   V1.6g 18.06.96 Moved FindZone() to FindZonePDB() in Bioplib
   V1.7  23.07.96 Supports atom wildcards. Some comment tidying.
   V1.7a 07.11.96 Added -ve Bvalues to mean ignore < cutoff
   V1.7b 11.11.96 Checks actual value of gUseBVal
   V1.7c 18.11.96 Added option to ignore missing atoms
   V1.7d 20.12.96 Added setting of gNFittedCoor
   V1.7e 27.06.97 Allows WRITE and RESIDUE to output to a pipe
   V1.7f 03.07.97 Added break into CreateFitArrays() to fix core dump
                  on bad multiple-occupancy PDB files
   V1.8  07.05.98 Skipped for release
   V2.0  01.03.01 Now supports multiple structure fitting and iterative
                  zone updating
   V2.1  28.03.01 Parameter for ITERATE and added CENTRE command
   V2.2  20.12.01 Skipped for release
   V2.3  01.12.04 Skipped for release
   V2.4  03.06.05 Skipped for release
   V2.5  07.06.05 Skipped for release
   V2.6  18.03.08 Added CentreOnZone() By: CTP
   V2.6  27.10.08 Added NoFitStructures() and DoNoFit().
   V3.0  06.11.08 Release Version
   V3.0  16.02.09 Rewrote CalculateRotationMatrix().
   V3.1  31.03.09 Skipped for release

*************************************************************************/
/* Includes
*/
#include "ProFit.h"




/************************************************************************/
/*>void FitStructures(void)
   ------------------------
   Fits the 2 structures using the currently defined ranges and displays
   the RMSd.

   28.09.92 Framework
   29.09.92 Various subsidiary bits added, still doesn't actually do
            fitting
   30.09.92 Added NOT option and calls to do actual fitting.
   17.07.95 Changed screen() to printf()
   18.07.95 Added initialisation of inserts in zones
   12.01.01 Moved ShowRMS() out of DoFitting() into here
   15.01.01 Added iteration of fitting zones
   01.02.01 Added multi-structure fitting and iteration of multiple
            structures
   15.02.01 Was iterating over multiple structures even when there
            were only 2.
   03.04.08 Added parameter to ShowRMS()
   10.06.08 disabled error check for multiple chains when using 
            iterative fitting. By: CTP
   29.08.08 added a final iteration for multi during which the coordinates
            for the averaged reference structure are NOT updated. 
   05.09.08 Final iteration for multiple structure fitting does not 
            update fitting zones when using iterative zone updating. 
*/
void FitStructures(void)
{
   ZONE  *z1,
         *z2;
   int   atmnum,
         NCoor,
         strucnum,
         niter;
   REAL  rmstot,
         rmsprev = (-100.0),
         deltaRMS;
   BOOL  final = FALSE;
   
   
   gFitted = FALSE;
   
   if(!gRefFilename[0])
   {
      printf("   Error==> Reference structure undefined.\n");
      return;
   }
   if(!gMobFilename[0][0])
   {
      printf("   Error==> Mobile structure undefined.\n");
      return;
   }

   /* 10.06.08 Error check for multiple chains disabled as this is now
      supported By: CTP          
   */
/***   
   if(gIterate)
   {
      // Check for numbers of chains 
      if(countchar(gRefSeq,'*') > 0)
      {
         printf("   Error==> Structures must have only one chain \
for iterative zones\n");
         return;
      }
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         if(countchar(gMobSeq[strucnum],'*') > 0)
         {
            printf("   Error==> Structures must have only one \
chain for iterative zones\n");
            return;
         }
      }
   }
***/
   
   if(!gQuiet)
   {
      printf("   Fitting structures...\n");
   }

   /* First copy the zones for display to match those for fitting       */
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      if(gRZoneList[strucnum] != NULL)
      {
         FREELIST(gRZoneList[strucnum],ZONE);
         gRZoneList[strucnum] = NULL;
      }
      for(z1=gZoneList[strucnum]; z1!=NULL; NEXT(z1))
      {
         /* Allocate an entry in RMS zone list                          */
         if(gRZoneList[strucnum])
         {
            /* Move to end of zone list                                 */
            z2=gRZoneList[strucnum];
            LAST(z2);
            ALLOCNEXT(z2,ZONE);
         }
         else
         {
            INIT(gRZoneList[strucnum],ZONE);
            z2 = gRZoneList[strucnum];
         }
         
         if(z2==NULL)
         {
            printf("   Error==> No memory for RMS zone!\n");
         }
         else
         {
            /* Add this zone to the RMS zone list                       */
            z2->chain1       = z1->chain1;
            z2->start1       = z1->start1;
            z2->startinsert1 = z1->startinsert1;
            z2->stop1        = z1->stop1;
            z2->stopinsert1  = z1->stopinsert1;
            z2->chain2       = z1->chain2;
            z2->start2       = z1->start2;
            z2->startinsert2 = z1->startinsert2;
            z2->stop2        = z1->stop2;
            z2->stopinsert2  = z1->stopinsert2;
            z2->mode         = z1->mode;
         }
      }
   }

   /* Now copy the atoms for RMS calculation                            */
   gNOTRMSAtoms = gNOTFitAtoms;
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      strcpy(gRMSAtoms[atmnum],gFitAtoms[atmnum]);

   if(gMultiCount > 1)
   {
      /* Keep looping counting the iterations                           */
      for(niter=0; ; niter++)
      {
         printf ("   Multi-structure fit iteration %d\n", niter);
         
         rmstot = (REAL)0.0;
         
         /* Loop through the structures we are fitting                  */
         for(strucnum=0; strucnum<gMultiCount; strucnum++)
         {
            /* Set up arrays for fitting                                */
            if((NCoor=CreateFitArrays(strucnum))!=0)
            {
               /* Reset the convergence criterion                       */
               CheckForConvergence(0, strucnum);
               
               /* Perform the fit                                       */
               if(DoFitting(NCoor, strucnum))
               {
                  while(gIterate && !final)
                  {
                     if((NCoor = UpdateFitArrays(strucnum))!=0)
                     {
                        if(!DoFitting(NCoor, strucnum))
                           return;
                        if(CheckForConvergence(NCoor, strucnum))
                           break;
                     }
                     else
                     {
                        break;
                     }
                  }
                  
                  /* Find RMS - Do not update during final iteration    */
                  rmstot += ShowRMS(FALSE,NULL,strucnum,!final,FALSE);
                  if(gIterate && !gQuiet)
                  {
                     printf("      (Over %d equivalenced CA-atoms)\n",
                            NCoor);
                  }
               }
            }
         }
         deltaRMS = (rmstot - rmsprev);
         rmsprev = rmstot;
         
         /* If we've converged or done too many iterations, do a final  */
         /* iteration then break out.                                   */
         if(final) break;

         if((ABS(deltaRMS) < MULTI_ITER_STOP) ||
            (niter > MAXMULTIITER))
            final = TRUE;
      }
   }
   else
   {
      /* Set up arrays for fitting                                      */
      if((NCoor=CreateFitArrays(FALSE))!=0)
      {
         /* Reset the convergence criterion                             */
         CheckForConvergence(0,0);
         
         /* Perform the fit                                             */
         DoFitting(NCoor, 0);

         while(gIterate)
         {
            printf("Iterating fit zones\n");
            if((NCoor = UpdateFitArrays(0))!=0)
            {
               DoFitting(NCoor, 0);
               if(CheckForConvergence(NCoor, 0))
                  break;
            }
            else
            {
               break;
            }
         }

         ShowRMS(FALSE,NULL,0,FALSE,FALSE);
         if(gIterate && !gQuiet)
         {
            printf("      (Over %d equivalenced CA-atoms)\n",
                   NCoor);
         }
      }
   }
   
   return;
}


/************************************************************************/
/*>void NoFitStructures(void)
   --------------------------
   Based on FitStructures(). Sets fitting zones but doesn't perform a fit 
   of the structures. Called by NOFIT. Used when wanting to perform
   RMS,MATRIX,etc... without fitting first.

   27.10.08 Original based on FitSructures() By: CTP
*/
void NoFitStructures(void)
{
   ZONE  *z1,
         *z2;
   int   atmnum,
         strucnum;

   /* Set gFitted and gNFittedCoor                                      */
   gFitted      = TRUE;
   gNFittedCoor = 0;
   
   if(!gRefFilename[0])
   {
      printf("   Error==> Reference structure undefined.\n");
      return;
   }
   if(!gMobFilename[0][0])
   {
      printf("   Error==> Mobile structure undefined.\n");
      return;
   }
   
   if(!gQuiet)
   {
      if(gMultiCount == 1)
         printf("   Mobile structure marked as fitted...\n");
      else 
         printf("   Mobile structures marked as fitted...\n");
   }

   /* First copy the zones for display to match those for fitting       */
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      if(gRZoneList[strucnum] != NULL)
      {
         FREELIST(gRZoneList[strucnum],ZONE);
         gRZoneList[strucnum] = NULL;
      }
      for(z1=gZoneList[strucnum]; z1!=NULL; NEXT(z1))
      {
         /* Allocate an entry in RMS zone list                          */
         if(gRZoneList[strucnum])
         {
            /* Move to end of zone list                                 */
            z2=gRZoneList[strucnum];
            LAST(z2);
            ALLOCNEXT(z2,ZONE);
         }
         else
         {
            INIT(gRZoneList[strucnum],ZONE);
            z2 = gRZoneList[strucnum];
         }
         
         if(z2==NULL)
         {
            printf("   Error==> No memory for RMS zone!\n");
         }
         else
         {
            /* Add this zone to the RMS zone list                       */
            z2->chain1       = z1->chain1;
            z2->start1       = z1->start1;
            z2->startinsert1 = z1->startinsert1;
            z2->stop1        = z1->stop1;
            z2->stopinsert1  = z1->stopinsert1;
            z2->chain2       = z1->chain2;
            z2->start2       = z1->start2;
            z2->startinsert2 = z1->startinsert2;
            z2->stop2        = z1->stop2;
            z2->stopinsert2  = z1->stopinsert2;
            z2->mode         = z1->mode;
         }
      }
   }

   /* Now copy the atoms for RMS calculation                            */
   gNOTRMSAtoms = gNOTFitAtoms;
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      strcpy(gRMSAtoms[atmnum],gFitAtoms[atmnum]);

   /* Copy the mobile PDB linked list to the rotation list. */
   if(gMultiCount > 1)
   {
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         DoNoFitting(strucnum);
      }
   }
   else
   {
      DoNoFitting(0);
   }
   
   return;
}


/************************************************************************/
/*>BOOL DoFitting(int NCoor, int strucnum)
   ---------------------------------------
   Does the actual fitting of the coordinate arrays.

   30.09.92 Original
   01.10.92 Corrected rotation procedure (!)
   08.10.93 Modified for new version of matfit().
            RotatePDB() -> ApplyMatrixPDB()
   17.07.95 Changed screen() to printf()
   19.07.95 Added parameter to ShowRMS()
   20.07.95 Added Weighted fitting
            Separate COOR structure for CofG so we don't corrupt the
            global version (allows printing of translation vector).
   25.07.95 Added another parameter to ShowRMS()
   13.06.96 Added B-value inverting for BWEIGHT command
   20.12.96 Added setting of gNFittedCoor
   12.01.01 gMobPDB[] and gFitPDB[] now arrays
            Moved ShowRMS() up to FitStructures()
   15.01.01 Now returns success/failure
   01.02.01 Added strucnum parameter
   20.02.01 gMobCofG now an array
            RotMat now local and a copy made into gRotMat
   13.08.08 Modified to Fit, Rotate & Refit. By: CTP
*/
BOOL DoFitting(int NCoor, int strucnum)
{
   PDB   *p,
         *q;
   VEC3F CofG;
   REAL  RotMat[3][3];
   int   i, j;

   gNFittedCoor = 0;
   
   if(NCoor < 3)
   {
      printf("   Error==> Fewer than 3 points to fit\n");
      return(FALSE);
   }
   else
   {
      /* If we are using inverse B-values, then invert the weights array*/
      if(gDoWeights==WEIGHT_INVBVAL)
      {
         int i;
         for(i=0; i<NCoor; i++)
            gWeights[i] = (REAL)1.0 / gWeights[i];
      }
      
      matfit(gRefCoor,gMobCoor[strucnum],RotMat,NCoor,
             ((gDoWeights!=WEIGHT_NONE)?gWeights:NULL),0);
     
      /* If we used inverse B-values, invert the weights back again     */
      if(gDoWeights==WEIGHT_INVBVAL)      
      {
         int i;
         for(i=0; i<NCoor; i++)
            gWeights[i] = (REAL)1.0 / gWeights[i];
      }
      
      /* Now copy the mobile PDB linked list to the rotation list       */
      if(gFitPDB[strucnum] != NULL) FREELIST(gFitPDB[strucnum], PDB);
      gFitPDB[strucnum] = NULL;
      
      for(p=gMobPDB[strucnum], q=NULL; p!=NULL; NEXT(p))
      {
         if(q==NULL)
         {
            INIT(gFitPDB[strucnum],PDB);
            q = gFitPDB[strucnum];
         }
         else
         {
            ALLOCNEXT(q,PDB);
         }

         if(q==NULL)
         {
            printf("   Error==> No memory for creating fitted \
structure.\n");
            return(FALSE);
         }
         
         /* Copy the coordinates                                        */
         CopyPDB(q, p);
      }
      
      /* Now we can rotate the rotation list                            */
      CofG.x = -1.0 * gMobCofG[strucnum].x;
      CofG.y = -1.0 * gMobCofG[strucnum].y;
      CofG.z = -1.0 * gMobCofG[strucnum].z;
      TranslatePDB(gFitPDB[strucnum], CofG);

      ApplyMatrixPDB(gFitPDB[strucnum], RotMat);

      TranslatePDB(gFitPDB[strucnum], gRefCofG);

      /* Record a copy of the rotation matrix for display with the MATRIX
         command
      */
      for(i=0; i<3; i++)
      {
         for(j=0; j<3; j++)
         {
            gRotMat[strucnum][i][j] = RotMat[i][j];
         }
      }

      gFitted = TRUE;
      gNFittedCoor = NCoor;
   }

#ifdef ROTATE_REFIT
   /* Rotate and Refit                                                  */
   /* ================                                                  */
   
   /* This was an attempt yo get around the saddle point local minimum
      problem in which fitting could convert 180degrees from the true
      minumum. This is now fixed in fit.c
   */
   if(!gIterate)
   {
      /* Set Variables                                                  */
      BOOL  multivsref = gMultiVsRef;
      REAL  rmsd_a = 0.0;
      REAL  rmsd_b = 0.0;
      
      REAL  rotmat_repos[3][3],
            rotmat_refit[3][3],
            rotmat_final[3][3];
      
      PDB   *ReFitPDB = NULL,
            *SwapPDB  = NULL;
      
      /* Set RMSD Calc vs Ref                                           */
      gMultiVsRef = TRUE;
      
      /* Calculate RMSD -A-                                             */
      rmsd_a = CalcRMS(FALSE,NULL,strucnum,FALSE,FALSE);
      
      /* Calculate Reposition Matrix                                    */
      MatMult33_33(gRotMat[strucnum],gRotMatTwist,rotmat_repos);
      
      /* Rotate Coordinates to New Position                             */
      ApplyMatrixCOOR(gMobCoor[strucnum], rotmat_repos, NCoor);
      
      /* Re-Fit                                                         */
      matfit(gRefCoor,gMobCoor[strucnum],rotmat_refit,NCoor,
             ((gDoWeights!=WEIGHT_NONE)?gWeights:NULL),0);
      
      /* Calculate Final Matrix                                         */
      MatMult33_33(rotmat_repos,rotmat_refit,rotmat_final);
      
      /* Make New Fitted PDB                                            */
      for(p=gMobPDB[strucnum], q=NULL; p!=NULL; NEXT(p))
      {
         if(q==NULL)
         {
            INIT(ReFitPDB,PDB);
            q = ReFitPDB;
         }
         else
         {
            ALLOCNEXT(q,PDB);
         }
         
         if(q==NULL)
         {
            printf("   Error==> No memory for creating fitted \
structure.\n");
            return(FALSE);
         }
         
         /* Copy the coordinates                                        */
         CopyPDB(q, p);
      }
      
      /* Now we can rotate the rotation list                            */
      TranslatePDB(ReFitPDB, CofG);
      ApplyMatrixPDB(ReFitPDB, rotmat_final);      
      TranslatePDB(ReFitPDB, gRefCofG);
      
      /* Calculate RMSD -B-                                             */
      SwapPDB = gFitPDB[strucnum];
      gFitPDB[strucnum] = ReFitPDB;
      rmsd_b = CalcRMS(FALSE,NULL,strucnum,FALSE,FALSE);
      gFitPDB[strucnum] = SwapPDB;
      
      /* Reset RMSD Calc vs Ref                                         */
      gMultiVsRef = multivsref;
      
      /* Select best result                                             */
      if(rmsd_a <= rmsd_b)
      {
         /* Free ReFitPDB                                               */
         if(ReFitPDB != NULL) FREELIST(ReFitPDB, PDB);
      }
      else
      {
         /* Free gFitPDB[strucnum] point to ReFitPDB                    */
         if(gFitPDB[strucnum] != NULL) FREELIST(gFitPDB[strucnum], PDB);
         gFitPDB[strucnum] = ReFitPDB;
         
         /* Copy New Rotation Matrix                                    */
         for(i=0; i<3; i++)
         {
            for(j=0; j<3; j++)
            {
               gRotMat[strucnum][i][j] = rotmat_final[i][j];
            }
         }
      }
   }
#endif   /* ROTATE_REFIT                                                */

   gFitted = TRUE;
   gNFittedCoor = NCoor;

   return(TRUE);
}


/************************************************************************/
/*>BOOL DoNoFitting(int strucnum)
   ------------------------------
   Sets gFitPDB if not already set. Sets rotation matrix to identity.

   27.10.08 Original based on DoFitting() By: CTP
   29.10.08 Tidied code.
   11.11.08 Simplified function - function no longer resets to starting 
            structures.
*/
BOOL DoNoFitting(int strucnum)
{
   PDB   *p,
         *q;
   int   i, j;

   gNFittedCoor = 0;

   /* Set rotation matrix to identity                                   */
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         gRotMat[strucnum][i][j] = (REAL)0.0;
      }
      gRotMat[strucnum][i][i] = (REAL)1.0;
   }
   
   /* Reset centres of geometry                                         */
   gRefCofG.x =  0.0;
   gRefCofG.y =  0.0;
   gRefCofG.z =  0.0;
   gMobCofG[strucnum].x = 0.0;
   gMobCofG[strucnum].y = 0.0;
   gMobCofG[strucnum].z = 0.0;
   
   /* Return if gFitPDB exists                                          */
   if(gFitPDB[0] != NULL)
   {
      gFitted = TRUE;
      return(TRUE);
   }
   
   /* Copy the mobile PDB linked list to the rotation list              */
   if(gFitPDB[strucnum] != NULL) 
      FREELIST(gFitPDB[strucnum], PDB);
   gFitPDB[strucnum] = NULL;
   
   for(p=gMobPDB[strucnum], q=NULL; p!=NULL; NEXT(p))
   {
      if(q==NULL)
      {
         INIT(gFitPDB[strucnum],PDB);
         q = gFitPDB[strucnum];
      }
      else
      {
         ALLOCNEXT(q,PDB);
      }
      
      if(q==NULL)
      {
         printf("   Error==> ");
         printf("No memory for creating structure.\n");
         return(FALSE);
      }
      
      /* Copy the coordinates                                           */
      CopyPDB(q, p);
   }
   
   gFitted = TRUE; 
   return(TRUE);
}


/************************************************************************/
/*>int ValidAtom(char *atnam, int mode)
   ------------------------------------
   Tests whether this atoms is in the appropriate list

   30.09.92 Original
   23.07.96 Calls AtomNameMatch() rather than strncmp; this handles
            wildcards in atom names
   15.02.01 Calls AtomNameRawMatch() instead of AtomNameMatch()
*/

int ValidAtom(char *atnam, int mode)
{
   int  DefReturn = FALSE,
        j;
   BOOL ErrorWarn;
   
   if(mode == ATOM_FITTING)
   {
      /* First check for all atoms                                      */
      if(gFitAtoms[0][0] == '*') return(TRUE);
      
      for(j=0;j<NUMTYPES;j++)
      {
         if(gFitAtoms[j][0] == '\0') break;
         
         if(gNOTFitAtoms)
         {
            DefReturn = TRUE;
            ErrorWarn = TRUE;
            if(AtomNameRawMatch(atnam,gFitAtoms[j],&ErrorWarn))
               return(FALSE);
         }
         else
         {
            ErrorWarn = TRUE;
            if(AtomNameRawMatch(atnam,gFitAtoms[j],&ErrorWarn))
               return(TRUE);
         }
      }
   }

   if(mode == ATOM_RMS)
   {
      /* First check for all atoms                                      */
      if(gRMSAtoms[0][0] == '*') return(TRUE);
      
      for(j=0;j<NUMTYPES;j++)
      {
         if(gRMSAtoms[j][0] == '\0') break;
         
         if(gNOTRMSAtoms)
         {
            DefReturn = TRUE;
            ErrorWarn = TRUE;
            if(AtomNameRawMatch(atnam,gRMSAtoms[j],&ErrorWarn))
               return(FALSE);
         }
         else
         {
            ErrorWarn = TRUE;
            if(AtomNameRawMatch(atnam,gRMSAtoms[j],&ErrorWarn))
               return(TRUE);
         }
      }
   }
   
   return(DefReturn);
}


/************************************************************************/
/*>REAL CalcRMS(BOOL ByRes, FILE *fp, int strucnum, BOOL UpdateReference,
                BOOL ByAtm)
   ----------------------------------------------------------------------
   Calculates RMS over currently defined RMS zones and atoms and prints
   it.

   30.09.92 Original
   01.10.92 Added check on NULL coordinates. Fix to finding mobile atoms.
   17.07.95 Changed screen() to printf()
   18.07.95 Removed zeroing of CoorCount which broke multi-zone fitting
            as in CreateFitArrays()
            Added initialisation of inserts in zones
            Added calls to FormatZone()
   19.07.95 Added ByRes parameter and calculation. Now prints the RMS
            and returns type void rather than the RMS
   25.07.95 Added fp parameter
   31.07.95 Added printing of number of residues if mismatch
   14.08.95 Was skipping the last residue when printing RMS by residue
   31.05.96 Added test on b-value
   18.06.96 Replaced MODE_* with ZONE_MODE_*
            Replaced FindZone() with FindZonePDB()
   06.11.96 Negative BVal cutoff interpreted as > bval
   11.11.96 Checks actual value of gUseBVal
   18.11.96 Added check on gIgnoreMissing
   11.01.01 gFitPDB now an array
   15.01.01 Now returns the RMS. Checks fp is non-null before printing
            Returns (-1.0) if error or doing by-residue RMSD.
   01.02.01 Added UpdateReference parameter
   20.02.01 -999 for start or end of structure rather than -1
   28.02.01 Fixed gRZoneList[0] to gRZoneList[strucnum]
   03.04.08 Modified to print distances between equivalenced atom pairs.
            Added parameter ByAtm to control printing. By: CTP
   07.04.08 Added include pair based on atom distance.
   16.04.08 Atom pairs outside of cutoff are marked in output.
   17.04.08 Residues partially/fully outside distance cutoff are marked.
   29.07.08 Multi-structure fitting updated. Fitting is against averaged 
            reference but simple rms/distance calculations are against 
            the first mobile structure.
   12.09.08 Added weighted average option for updating reference.
   23.10.08 Added gWtAverage flag.
   07.11.08 Simple rms/distance calculations are against the mobile 
            structure indicated by gMultiRef.
*/
REAL CalcRMS(BOOL ByRes, FILE *fp, int strucnum, BOOL UpdateReference,
             BOOL ByAtm)
{
   REAL  SumSq       = 0.0,
         rms         = 0.0;
   PDB   *refpdblist = NULL,
         *ref_start  = NULL,
         *ref_stop   = NULL,
         *fit_start  = NULL,
         *fit_stop   = NULL,
         *prevp      = NULL,
         *prevq      = NULL,
         *p,
         *q,
         *r,
         *m;
   ZONE  *z;
   char  ref_insert,
         fit_insert;
   int   ref_resnum,
         fit_resnum,
         ref_nres,
         fit_nres,
         CoorCount = 0,
         Found,
         CoorOutside = 0;


   /* If no zones have been specified, create a single all atoms zone   */
   if(gRZoneList[strucnum] == NULL)
   {
      INIT(gRZoneList[strucnum],ZONE);
      gRZoneList[strucnum]->chain1        = ' ';
      gRZoneList[strucnum]->start1        = -999;
      gRZoneList[strucnum]->startinsert1  = ' ';
      gRZoneList[strucnum]->stop1         = -999;
      gRZoneList[strucnum]->stopinsert1   = ' ';
      gRZoneList[strucnum]->chain2        = ' ';
      gRZoneList[strucnum]->start2        = -999;
      gRZoneList[strucnum]->startinsert2  = ' ';
      gRZoneList[strucnum]->stop2         = -999;
      gRZoneList[strucnum]->stopinsert2   = ' ';
      gRZoneList[strucnum]->mode          = gCurrentMode;
   }
   
   /* Step through each zone                                            */
   for(z=gRZoneList[strucnum]; z!=NULL; NEXT(z))
   {
      /* Set Reference Zone                                             */
      /*  - Multi-structure fitting will fit to an averaged structure   */
      /*  - Simple RMS or Distance calculations are against the mobile  */
      /*    structure set as the reference (gMultiRef) or the averaged  */
      /*    reference depending on gMultiVsRef.                         */
/***
      refpdblist = (gMultiCount > 1 && !UpdateReference && !gMultiVsRef) 
        ? gFitPDB[0] : gRefPDB;
***/
      refpdblist = (gMultiCount > 1 && !UpdateReference && !gMultiVsRef) 
                   ? gFitPDB[gMultiRef] : gRefPDB;
      
      /* Reference structure                                            */
      /*
        if(!FindZonePDB(gRefPDB, z->start1, z->startinsert1, 
                        z->stop1, z->stopinsert1, z->chain1, z->mode, 
                        &ref_start, &ref_stop))
      */
      if(!FindZonePDB(refpdblist, z->start1, z->startinsert1, 
                      z->stop1, z->stopinsert1, z->chain1, z->mode, 
                      &ref_start, &ref_stop))
      {
         char zone1[64],
              zone2[64];

         /* Check ranges have been found                                */
         printf("   Error==> Reference structure zone not found.\n");
         
         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         printf("      %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                              :"(Sequential numbering)"));

         return(-1.0);
      }

      /* Mobile structure                                               */
      if(!FindZonePDB(gFitPDB[strucnum], z->start2, z->startinsert2, 
                      z->stop2, z->stopinsert2, z->chain2, z->mode, 
                      &fit_start, &fit_stop))
      {
         char zone1[64],
              zone2[64];

         /* Check ranges have been found                                */
         printf("   Error==> Mobile structure zone not found.\n");

         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         printf("      %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                              :"(Sequential numbering)"));

         return(-1.0);
      }
      
      
      /* Check we have the same number of residues in each zone         */
      ref_nres     = 1;
      ref_resnum   = ref_start->resnum;
      ref_insert   = ref_start->insert[0];
      for(p=ref_start; p!=ref_stop; NEXT(p))
      {
         if(p->resnum != ref_resnum || p->insert[0] != ref_insert)
         {
            ref_nres++;
            ref_resnum = p->resnum;
            ref_insert = p->insert[0];
         }
      }

      fit_nres     = 1;
      fit_resnum   = fit_start->resnum;
      fit_insert   = fit_start->insert[0];
      for(p=fit_start; p!=fit_stop; NEXT(p))
      {
         if(p->resnum != fit_resnum || p->insert[0] != fit_insert)
         {
            fit_nres++;
            fit_resnum = p->resnum;
            fit_insert = p->insert[0];
         }
      }
      
      if(ref_nres != fit_nres)
      {
         char zone1[64],
              zone2[64];
         
         printf("   Error==> Number of residues in zone does not \
match.\n");
         /* Added 13.12.95                                              */
         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         
         printf("           %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                         :"(Sequential numbering)"));

         /* Added 31.07.95                                              */
         printf("            Reference: %d, Mobile: %d\n",
                ref_nres, fit_nres);


         return(-1.0);
      }

      /* Insert the atoms from this zone into the coordinate arrays     */
/*    Removed 17.07.95.....
//    CoorCount  = 0;   
*/

      ref_nres   = 0;
      ref_resnum = -999;
      ref_insert = ' ';

      for(p=ref_start; p!=ref_stop; NEXT(p))
      {
         /* Is this the start of a new reference residue?               */
         if(p->resnum != ref_resnum || p->insert[0] != ref_insert)
         {
            /* If we're displaying by residue and at least one residue
               has been processed, display the RMS for that (previous)
               residue. Reset other value to zero.
            */
            if(ByRes && ref_nres && CoorCount && fp!=NULL && !ByAtm)
            {
               char buff1[16],
                    buff2[16];

               SumSq /= CoorCount;
               rms    = (REAL)sqrt((double)SumSq);

               sprintf(buff1,"%c%d%c",
                       prevp->chain[0], prevp->resnum, prevp->insert[0]);
               sprintf(buff2,"%c%d%c",
                       prevq->chain[0], prevq->resnum, prevq->insert[0]);
               
               fprintf(fp,"%8s %s : %8s %s   RMS: %.3f",
                       buff1,
                       prevp->resnam,
                       buff2,
                       prevq->resnam,
                       rms);

               /* Flag Residues Outside Distance Cutoff By: CTP         */
               if(gUseDistCutoff && CoorOutside)
               {
                  if(CoorOutside < CoorCount)
                     fprintf(fp," *\n");  /* Partially within cutoff    */
                  else
                     fprintf(fp," **\n"); /* Outside cutoff             */
               }
               else 
                  fprintf(fp,"\n");        /* Inside cutoff             */
                 
               /* Zero the sum of squares and the coordinate count      */
               SumSq       = (REAL)0.0;
               CoorCount   = 0;
               CoorOutside = 0;
            }
            
            ref_nres++;
            ref_resnum = p->resnum;
            ref_insert = p->insert[0];

            fit_nres   = 0;
            fit_resnum = -999;
            fit_insert = ' ';
         
            for(q=fit_start; q!=fit_stop; NEXT(q))
            {
               if(q->resnum != fit_resnum || q->insert[0] != fit_insert)
               {
                  /* Start of a new mobile residue                      */
                  fit_nres++;
                  fit_resnum = q->resnum;
                  fit_insert = q->insert[0];
   
                  /* Consider the equivalent residues                   */
                  if(ref_nres == fit_nres)
                  {
                     /* p points to the start of a residue in reference,
                        q points to the start of a residue in mobile

                        Step through reference set
                     */
                     for(r=p;
                         r!=NULL && r->resnum==ref_resnum 
                                 && r->insert[0]==ref_insert;
                         NEXT(r))
                     {
                        /* 15.02.01 Changed from atnam to atnam_raw     */
                        if(ValidAtom(r->atnam_raw, ATOM_RMS))
                        {
                           if(r->x == 9999.0 &&
                              r->y == 9999.0 &&
                              r->z == 9999.0)
                           {
                              /* Ignore NULL atoms                      */
                              continue;
                           }

                           /* 31.05.96 Test the BVal
                              06.11.96 Modified to handle -ve 
                                       specifications
                              11.11.96 Checks actual value of gUseBVal
                           */
                           if(gUseBVal==1 || gUseBVal==2)
                           {
                              if(gBValue >= (REAL)0.0)
                              {
                                 if(r->bval > gBValue)
                                    continue;
                              }
                              else
                              {
                                 if(-(r->bval) > gBValue)
                                    continue;
                              }
                           }
                           
                           /* Find this atom in the mobile set          */
                           Found = FALSE;
                           for(m=q;
                               m!=NULL && m->resnum==fit_resnum 
                                       && m->insert[0]==fit_insert;
                               NEXT(m))
                           {
                              /* 28.02.01 Changed from atnam...         */
                              if(!strcmp(r->atnam_raw,m->atnam_raw))
                              {
                                int pair_total = 1; 
                                int i = 0;
                                 Found = TRUE;

                                 if(m->x == 9999.0 &&
                                    m->y == 9999.0 &&
                                    m->z == 9999.0)
                                 {
                                    /* Ignore NULL atoms                */
                                    continue;
                                 }


                                 /* 31.05.96 Test the BVal
                                    06.11.96 Modified to handle -ve 
                                             specifications
                                    11.11.96 Checks actual value of 
                                             gUseBVal
                                 */
                                 if(gUseBVal==1 || gUseBVal==3)
                                 {
                                    if(gBValue >= (REAL)0.0)
                                    {
                                       if(m->bval > gBValue)
                                          continue;
                                    }
                                    else
                                    {
                                       if(-(m->bval) > gBValue)
                                          continue;
                                    }
                                 }
                                 

                                 /* Auto-match Symmetrical Atoms By: CTP*/
                                 if(gMatchSymAtoms && 
                                    r->next != NULL && m->next != NULL)
                                 {
                                    int i = 0;
                                    for(i=0;i < SYMM_ATM_PAIRS ;i++)
                                    {
                                       if(gSymType[i][3][0] == TRUE &&
                                          !strcmp(gSymType[i][0],
                                                  r->resnam) && 
                                          !strcmp(gSymType[i][1],
                                                  r->atnam_raw) && 
                                          !strcmp(gSymType[i][2],
                                                  r->next->atnam_raw) &&
                                          !strcmp(gSymType[i][0],
                                                  m->resnam) && 
                                          !strcmp(gSymType[i][1],
                                                  m->atnam_raw) && 
                                          !strcmp(gSymType[i][2],
                                                  m->next->atnam_raw))
                                       {
                                          double sqdist_a = 0.0;
                                          double sqdist_b = 0.0;
                                          
                                          sqdist_a  = PDBDISTSQ(r,m);
                                          sqdist_a +=
                                             PDBDISTSQ(r->next, m->next);
                                          
                                          sqdist_b = PDBDISTSQ(r,m->next);
                                          sqdist_b +=
                                             PDBDISTSQ(r->next, m);
                                          
                                          /* Set Number Pairs to Process*/
                                          if(sqdist_a > sqdist_b)
                                             pair_total = 2;
                                          
                                          break;
                                       }
                                    }
                                 }

                                 /* Calculate RMSd  By: CTP             */
                                 for(i = 0; i < pair_total; i++)
                                 {
                                    
                                    /* Set atom pair to process         */
                                    PDB *r_atm = NULL;
                                    PDB *m_atm = NULL;
                                    
                                    if(pair_total == 1)
                                    {
                                       /* Process single pair           */
                                       r_atm = r; 
                                       m_atm = m;
                                    }
                                    else if(pair_total == 2 && i == 0)
                                    {
                                       /* Process pair one of two       */
                                       r_atm = r; 
                                       m_atm = m->next;
                                    }
                                    else
                                    {
                                       /* Process pair two of two       */
                                       r_atm = r->next; 
                                       m_atm = m;
                                    }
                                    
                                    /* Output Atom Distances            */
                                    if(ByRes && ByAtm)
                                    {
                                       double dist = 0.0;
                                       char res1[10];
                                       char res2[10];
                                       
                                       dist = 
                                          sqrt(PDBDISTSQ(r_atm,m_atm));
                                       
                                       sprintf(res1,"%c%4d%c",
                                               r_atm->chain[0],
                                               r_atm->resnum,
                                               r_atm->insert[0]);
                                       sprintf(res2,"%c%4d%c",
                                               m_atm->chain[0],
                                               m_atm->resnum,
                                               m_atm->insert[0]);
                                       
                                       fprintf(fp,"%8s %4s %4s :",res1,
                                               r_atm->resnam,
                                               r_atm->atnam_raw);
                                       fprintf(fp,"%8s %4s %4s ",res2,
                                               m_atm->resnam,
                                               m_atm->atnam_raw);
                                       fprintf(fp,"Dist: %.3f",dist);
                                       
                                       /* Flag pairs outside cutoff     */
                                       if(gUseDistCutoff &&
                                          dist > gDistCutoff )
                                          fprintf(fp," *");
                                       
                                       fprintf(fp,"\n");
                                    }

                                    /* Include pair based on atom 
                                       distance 
                                    */
                                    /* If displaying RMSd by residue then
                                       count coordinates outside cutoff.
                                    */
                                    if((gUseDistCutoff) &&
                                       (PDBDISTSQ(r_atm,m_atm) > 
                                        (gDistCutoff * gDistCutoff)))
                                    {
                                       if(ByRes)
                                       {
                                          CoorOutside++;
                                       }
                                       else 
                                       {
                                          continue;
                                       }
                                    }
                                    
                                    /* Sum Squares                      */
                                    SumSq     += PDBDISTSQ(r_atm,m_atm);
                                    CoorCount ++;
                                    
                                    
                                    /* If we are averaging the coordinates
                                       for multi-structure fitting, then
                                       update the reference atom 
                                       coordinates with the mean of that 
                                       and the mobile
                                    */
                                    
                                    /* Mean of Ref + Mobile Coordinates */
                                    if(UpdateReference && !gWtAverage)
                                    {
                                       r_atm->x = (r_atm->x+m_atm->x)/2.0;
                                       r_atm->y = (r_atm->y+m_atm->y)/2.0;
                                       r_atm->z = (r_atm->z+m_atm->z)/2.0;
                                    }
                                    
                                    /* Weighted Mean of Ref + Mobile 
                                       Coords 
                                    */
                                    if(UpdateReference && gWtAverage)
                                    {
                                       REAL MultiCount = gMultiCount;
                                       
                                       r_atm->x = ((MultiCount - 1) * 
                                                   r_atm->x 
                                                   + m_atm->x)/MultiCount;
                                       r_atm->y = ((MultiCount - 1) * 
                                                   r_atm->y 
                                                   + m_atm->y)/MultiCount;
                                       r_atm->z = ((MultiCount - 1) * 
                                                   r_atm->z 
                                                   + m_atm->z)/MultiCount;
                                    }
                                    
                                 }
                                 
                                 /* Set Pointers                        */
                                 if(pair_total == 2)
                                 {
                                    r     = r->next;
                                    m     = m->next;
                                 }
                                 prevq = m;
                                 
                              }
                           }
                           
                           if(!Found)
                           {
                              if(!gIgnoreMissing)
                              {
                                 printf("   Error==> Atoms do not match \
in residue:\n");
                                 printf("      Reference %4s %5d%c \
Mobile %4s %5d%c\n",
                                        p->resnam,p->resnum,p->insert[0],
                                        q->resnam,q->resnum,q->insert[0]);
                                 /* 28.02.01 Changed from r->atnam      */
                                 printf("      Unable to find reference \
atom %4s in mobile.\n", r->atnam_raw);
                                 return(-1.0);
                              }
                           }
                        }
                     }  /* End of loop through reference residue        */
                     break;
                  }  /* End of if() equivalent residues                 */
               }  /* End of start-of-new-mobile-residue                 */
            }  /* End of loop through mobile set                        */
         }  /* End of start-of-new-reference-residue                    */
         prevp = p;
      }  /* End if loop through reference set                           */
      
      /* For the last residue, if we're displaying by residue and we have
         some atoms, display the RMS for this residue.
      */
      if(ByRes && ref_nres && CoorCount && fp!=NULL && !ByAtm)
      {
         char buff1[16],
              buff2[16];
         
         SumSq /= CoorCount;
         rms    = (REAL)sqrt((double)SumSq);
         
         sprintf(buff1,"%c%d%c",
                 prevp->chain[0], prevp->resnum, prevp->insert[0]);
         sprintf(buff2,"%c%d%c",
                 prevq->chain[0], prevq->resnum, prevq->insert[0]);
         
         fprintf(fp,"%8s %s : %8s %s   RMS: %.3f",
                 buff1,
                 prevp->resnam,
                 buff2,
                 prevq->resnam,
                 rms);
         
         /* Flag Residues Outside Distance Cutoff By: CTP               */
         if(gUseDistCutoff && CoorOutside)
         {
            if(CoorOutside < CoorCount)
               fprintf(fp," *\n");  /* Partially within cutoff          */
            else
               fprintf(fp," **\n"); /* Outside cutoff                   */
         }
         else 
         {
            fprintf(fp,"\n");       /* Inside cutoff                    */
         }
         
         /* Zero the sum of squares and the coordinate count            */
         SumSq     = (REAL)0.0;
         CoorCount = 0;
         CoorOutside = 0;
      }
   }  /* End of loop through zones                                      */
   
   if(!ByRes && CoorCount == 0) 
      printf("   Error==> No atoms in specified zones\n");
   
   /* Calculate the RMS                                                 */
#ifdef DEBUG
   fprintf(stderr,"\nCalculating RMS on %d atoms\n",CoorCount);
   fprintf(stderr,"Sum of squares = %8.3f\n",SumSq);
#endif
   
   /* Output legend if using distance cutoff By: CTP                    */
   if(gUseDistCutoff && ByRes)
   {
      if(ByAtm)
      {
         if(fp!=NULL)
            fprintf(fp,"   %35s * Outside distance cutoff\n","");
      }
      else
      {
         if(fp!=NULL)
         {
            fprintf(fp,"   %28s *  Partially outside distance cutoff.\n",
                    "");
            fprintf(fp,"   %28s ** Fully outside distance cutoff\n",
                    "");
         }
      }
   }
   
   /* Output RMSd By: CTP                                               */
   if(!ByRes)
   {
      SumSq /= CoorCount;
      rms    = (REAL)sqrt((double)SumSq);
      
      if(fp!=NULL)
         fprintf(fp,"   RMS: %.3f\n",rms);

      return(rms);
   }

   return(-1.0);
}


/************************************************************************/
/*>void ShowNFitted(void)
   ----------------------
   Displays the number of equivalent atom pairs used in the last fitting.

   20.12.96 Original   By: ACRM
*/
void ShowNFitted(void)
{
   if(gFitted && gNFittedCoor)
   {
      printf("   Number of fitted atoms: %d\n",gNFittedCoor);
   }
   else
   {
      printf("   Warning==> Structures have not yet been fitted.\n");
   }
}


/************************************************************************/
/*>REAL ShowRMS(BOOL ByRes, char *filename, int strucnum,
                BOOL UpdateReference, BOOL ByAtm)
   ------------------------------------------------------
   Display the RMS over the currently defined zones.

   29.09.92 Framework
   30.09.92 Original
   17.07.95 Changed screen() to printf()
            Handles inserts
   19.07.95 Added ByRes parameter passed to CalcRMS()
            RMS is now printed by CalcRMS()
   25.07.95 Added filename parameter
            Opens file if specified and passes FILE pointer to CalcRMS()
   27.06.97 Changed call to fopen() to OpenOrPipe
   01.02.01 Added strucnum and UpdateReference parameters
            Now returns the RMSD
   20.02.01 -999 for start or end of structure rather than -1
   03.04.08 Added parameter to ShowRMS() The parameter, ByAtm, turns on 
            printing of atom distances by CalcRMS(). By: CTP
*/
REAL ShowRMS(BOOL ByRes, char *filename, int strucnum, 
             BOOL UpdateReference, BOOL ByAtm)
{
   ZONE   *z1, *z2;
   int    atmnum;
   FILE   *fp = stdout;
   REAL   rmsd = (REAL)(-1.0);
   

   if(gFitted)
   {
      if(filename)
      {
         if((fp=OpenOrPipe(filename))==NULL)
         {
            printf("   Warning==> unable to open file for by-residue \
RMS\n");
            fp = stdout;
         }
      }

      /* Copy zones if user hasn't specified otherwise                  */
      if(!gUserRMSZone)
      {
         /* Free the current zone list                                  */
         if(gRZoneList[strucnum] != NULL)
         {
            FREELIST(gRZoneList[strucnum],ZONE);
            gRZoneList[strucnum] = NULL;
         }

         /* Create or copy a new one                                    */
         if(gZoneList[strucnum] == NULL)
         {
            INIT(gRZoneList[strucnum],ZONE);
            gRZoneList[strucnum]->chain1       = ' ';
            gRZoneList[strucnum]->start1       = -999;
            gRZoneList[strucnum]->startinsert1 = ' ';
            gRZoneList[strucnum]->stop1        = -999;
            gRZoneList[strucnum]->stopinsert1  = ' ';
            gRZoneList[strucnum]->chain2       = ' ';
            gRZoneList[strucnum]->start2       = -999;
            gRZoneList[strucnum]->startinsert2 = ' ';
            gRZoneList[strucnum]->stop2        = -999;
            gRZoneList[strucnum]->stopinsert2  = ' ';
            gRZoneList[strucnum]->mode         = gCurrentMode;
         }
         else
         {
            for(z1=gZoneList[strucnum]; z1!=NULL; NEXT(z1))
            {
               /* Allocate an entry in RMS zone list                    */
               if(gRZoneList[strucnum])
               {
                  /* Move to end of zone list                           */
                  z2=gRZoneList[strucnum];
                  LAST(z2);
                  ALLOCNEXT(z2,ZONE);
               }
               else
               {
                  INIT(gRZoneList[strucnum],ZONE);
                  z2 = gRZoneList[strucnum];
               }
            
               if(z2==NULL)
               {
                  printf("   Error==> No memory for RMS zone!\n");
               }
               else
               {
                  /* Add this zone to the RMS zone list                 */
                  z2->chain1       = z1->chain1;
                  z2->start1       = z1->start1;
                  z2->startinsert1 = z1->startinsert1;
                  z2->stop1        = z1->stop1;
                  z2->stopinsert1  = z1->stopinsert1;
                  z2->chain2       = z1->chain2;
                  z2->start2       = z1->start2;
                  z2->startinsert2 = z1->startinsert2;
                  z2->stop2        = z1->stop2;
                  z2->stopinsert2  = z1->stopinsert2;
                  z2->mode         = z1->mode;
               }
            }
         }
      }
   
      if(!gUserRMSAtoms)
      {
         /* Copy the atoms for RMS calculation if user hasn't specified */
         for(atmnum=0; atmnum<NUMTYPES; atmnum++)
            strcpy(gRMSAtoms[atmnum],gFitAtoms[atmnum]);
      }
   
      rmsd = CalcRMS(ByRes,fp, strucnum, UpdateReference, ByAtm);

      if(fp != stdout)
         CloseOrPipe(fp);
   }
   else
   {
      printf("   Warning==> Structures have not yet been fitted.\n");
   }

   return(rmsd);
}


/************************************************************************/
/*>int CheckForConvergence(int NCoor, int strucnum)
   ------------------------------------------------
   Checks whether the RMSD has converged or we have done too many 
   iterations.

   15.01.01 Original   By: ACRM
   01.02.01 Added strucnum parameter
   03.04.08 Added parameter to CalcRMS() By: CTP
*/
int CheckForConvergence(int NCoor, int strucnum)
{
   static REAL lastRMS = (REAL)(-1.0);
   static int  niter   = 0;
   REAL   rms;
   
   if(NCoor == 0)
   {
      lastRMS = (REAL)(-1.0);
      niter   = 0;
   }
   else
   {
      rms = CalcRMS(FALSE, NULL, strucnum, FALSE, FALSE);
      if(lastRMS >= (REAL)(0.0))
      {
         REAL deltaRMS = (rms - lastRMS);
         if((ABS(deltaRMS) < ITER_STOP) ||
            (niter > MAXITER))
            return(TRUE);
      }

      lastRMS = rms;
      niter++;
   }
   
   return(FALSE);
}


/************************************************************************/
/*>int UpdateFitArrays(int strucnum)
   ---------------------------------
   Update the fit arrays for iterative fitting by updating the equivalence
   list using DP

   15.01.01 Original   By: ACRM
*/
int UpdateFitArrays(int strucnum)
{
   int      length1,
            length2,
            align_len,
            NCoor;
   char     *ref_align = NULL,
            *mob_align = NULL;
   REAL     score;
   PDB      *RefCaPDB = NULL,
            *MobCaPDB = NULL,
            **RefIndex = NULL,
            **MobIndex = NULL;
   char     *sel[2];
         
   if(!gQuiet)
   {
      printf("      Updating Fitting Zones...\n");
   }

   /* Extract the CA atoms and index them so they can be accessed by 
      offset number
   */
   SELECT(sel[0],"CA  ");
   RefCaPDB = SelectAtomsPDB(gRefPDB, 1, sel, &length1);
   MobCaPDB = SelectAtomsPDB(gFitPDB[strucnum], 1, sel, &length2);
   RefIndex = IndexPDB(RefCaPDB, &length1);
   MobIndex = IndexPDB(MobCaPDB, &length2);
   
   /* Allocate memory for alignment sequences                           */
   if((ref_align = (char *)malloc((length1+length2)*sizeof(char)))==
      NULL)
   {
      printf("   Warning==> No memory for alignment!\n");
      return(0);
   }
   if((mob_align = (char *)malloc((length1+length2)*sizeof(char)))==
      NULL)
   {
      printf("   Warning==> No memory for alignment!\n");
      free(ref_align);
      return(0);
   }
   
   /* Perform the alignment                                             */
   score = AlignOnCADistances(RefIndex, length1,
                              MobIndex, length2,
                              ref_align, mob_align, &align_len);
   if(score <= (REAL)0.0)
   {
      printf("   Error==> Unable to perform alignment!\n");
      return(0);
   }
   
   /* Clear any current fitting zones                                   */
   SetFitZone("CLEAR", strucnum);
   
   /* Now set zones based on alignment                                  */
   SetNWZones(ref_align, mob_align, align_len, RefIndex, MobIndex, 
              strucnum);

   /* Create the fitting arrays from the zones                          */
   NCoor = CreateFitArrays(strucnum);
   
   /* Free allocated memory                                             */
   free(ref_align);
   free(mob_align);
   FREELIST(RefCaPDB, PDB);
   FREELIST(MobCaPDB, PDB);
   free(RefIndex);
   free(MobIndex);
   
   return(NCoor);
}

/************************************************************************/
/*>REAL Distance(PDB *p, PDB *q)
   -----------------------------
   Calculate a distance-based score for AlignOnCADistances() used to
   update the equivalence list. Score returned is 1/distance (with
   a minimum allowed distance of 0.0001

   15.01.01 Original   By: ACRM
*/
#define TINY 0.0001
REAL Distance(PDB *p, PDB *q)
{
   REAL dist;
   
   dist = DIST(p, q);
   
   if(dist < TINY)
      dist = TINY;

   return((REAL)1.0/dist);
}

/************************************************************************/
/*>REAL AlignOnCADistances(PDB **RefIndex, int length1, 
                           PDB **MobIndex, int length2,
                           char *align1, char *align2, int *align_len)
   -------------------------------------------------------------------
   Performs a DP alignment to find the updated equivalences on the basis
   of selecting closest distances

   15.01.01 Original   By: ACRM
*/
REAL AlignOnCADistances(PDB **RefIndex, int length1, 
                        PDB **MobIndex, int length2,
                        char *align1, char *align2, int *align_len)
{
   XY    **dirn   = NULL;
   int   maxdim,
         i,    j,    k,    l,
         i1,   j1,
         rcell, dcell;
   REAL  **matrix = NULL,
         thisscore,
         gapext,
         score, maxoff,
         dia,  right, down;

   /* gap penalties are set to zero - we don't care how many gaps we
      introduce
   */
   REAL  penalty = 0.0, 
         penext  = 0.0;

   maxdim = MAX(length1, length2);
   
   /* Initialise the score matrix                                       */
   if((matrix = (REAL **)Array2D(sizeof(REAL), maxdim, maxdim))==NULL)
      return(0);
   if((dirn   = (XY **)Array2D(sizeof(XY), maxdim, maxdim))==NULL)
      return(0);
      
   for(i=0;i<maxdim;i++)
   {
      for(j=0;j<maxdim;j++)
      {
         matrix[i][j] = (REAL)0.0;
         dirn[i][j].x = -1;
         dirn[i][j].y = -1;
      }
   }
    
   /* Fill in scores up the right hand side of the matrix               */
   for(j=0; j<length2; j++)
   {
      REAL dist;
      dist = Distance(RefIndex[length1-1], MobIndex[j]);
      matrix[length1-1][j] = dist;
   }

   /* Fill in scores along the bottom row of the matrix                 */
   for(i=0; i<length1; i++)
   {
      matrix[i][length2-1] = Distance(RefIndex[i], MobIndex[length2-1]);
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
         matrix[i1][j] += Distance(RefIndex[i1],MobIndex[j]);
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
         matrix[i][j1] += Distance(RefIndex[i],MobIndex[j1]);
      }
   } 
   
   score = TraceBackDistMat(matrix, dirn, length1, length2,
                            RefIndex, MobIndex, align1, align2, 
                            align_len);

#ifdef VERBOSE
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
#endif
    
   FreeArray2D((char **)matrix, maxdim, maxdim);
   FreeArray2D((char **)dirn,   maxdim, maxdim);
    
   return(score);
}



/************************************************************************/
/*>REAL TraceBackDistMat(int **matrix, XY **dirn, 
                         int length1, int length2, 
                         PDB **RefIndex, PDB **MobIndex, 
                         char *align1, 
                         char *align2, int *align_len)
   -----------------------------------------------------
   Input:   int  **matrix   N&W matrix
            XY   **dirn     Direction Matrix
            int  length1    Length of first sequence
            int  length2    Length of second sequence
            PDB  **RefIndex      First sequence
            ODB  **MobIndex      Second sequence
   Output:  char *align1    First sequence aligned
            char *align2    Second sequence aligned
            int  *align_len Aligned sequence length
   Returns: int             Alignment score

   Does the traceback to find the aligment.

   15.01.01 Original based on TraceBack()
*/
REAL TraceBackDistMat(REAL **matrix, 
                      XY   **dirn,
                      int  length1, 
                      int  length2, 
                      PDB  **RefIndex, 
                      PDB  **MobIndex, 
                      char *align1, 
                      char *align2, 
                      int  *align_len)
{
   int   i,    j, 
         ai, 
         BestI,BestJ;
   XY    nextCell;

   ai = SearchForBestDistMat(matrix, length1, length2, &BestI, &BestJ, 
                             RefIndex, MobIndex, align1, align2);

   /* Now trace back to find the alignment                              */
   i            = BestI;
   j            = BestJ;
   align1[ai]   = throne(RefIndex[i]->resnam);
   align2[ai++] = throne(MobIndex[j]->resnam);

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
            the y-sequence (MobIndex)
         */
         i++;
         j++;
         while((i < nextCell.x) && (i < length1-1))
         {
            align1[ai] = throne(RefIndex[i++]->resnam);
            align2[ai++] = '-';
         }
      }
      else if(nextCell.x == i+1)
      {
         /* We are inheriting from the off-diagonal inserting a gap in
            the x-sequence (RefIndex)
         */
         i++;
         j++;
         while((j < nextCell.y) && (j < length2-1))
         {
            align1[ai] = '-';
            align2[ai++] = throne(MobIndex[j++]->resnam);
         }
      }
      else
      {
         /* Cockup!                                                     */
         fprintf(stderr,"align.c/TraceBack() internal error\n");
      }
      
      align1[ai]   = throne(RefIndex[i]->resnam);
      align2[ai++] = throne(MobIndex[j]->resnam);
   }

   /* If one sequence finished first, fill in the end with insertions   */
   if(i < length1-1)
   {
      for(j=i+1; j<length1; j++)
      {
         align1[ai]   = throne(RefIndex[j]->resnam);
         align2[ai++] = '-';
      }
   }
   else if(j < length2-1)
   {
      for(i=j+1; i<length2; i++)
      {
         align1[ai]   = '-';
         align2[ai++] = throne(MobIndex[i]->resnam);
      }
   }
   
   *align_len = ai;
   
   return(matrix[BestI][BestJ]);
}


/************************************************************************/
/*>int SearchForBestDistMat(REAL **matrix, int length1, 
                            int length2, int *BestI, int *BestJ, 
                            PDB **RefIndex, PDB **MobIndex, 
                            char *align1, char *align2)
   -------------------------------------------------------------
   Input:   REAL **matrix   N&W matrix
            int  length1    Length of first sequence
            int  length2    Length of second sequence
            int  *BestI     x position of highest score
            int  *BestJ     y position of highest score
            PDB  **RefIndex First sequence
            PDB  **MobIndex Second sequence
   Output:  char *align1    First sequence with end aligned correctly
            char *align2    Second sequence with end aligned correctly
   Returns: int             Alignment length thus far

   Searches the outside of the matrix for the best score and starts the
   alignment by putting in any starting - characters.

   15.01.01 Original based on SearchForBest()
*/
int SearchForBestDistMat(REAL **matrix, 
                         int  length1, 
                         int  length2, 
                         int  *BestI, 
                         int  *BestJ,
                         PDB  **RefIndex, 
                         PDB  **MobIndex, 
                         char *align1, 
                         char *align2)
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
         align1[ai] = throne(RefIndex[i]->resnam);
         align2[ai++] = '-';
      }
   }
   else
   {
      *BestI = 0;
      *BestJ = bestj;
      for(j=0; j<*BestJ; j++)
      {
         align1[ai] = '-';
         align2[ai++] = throne(MobIndex[j]->resnam);
      }
   }
   return(ai);
}


/************************************************************************/
/*>int CreateFitArrays(int strucnum)
   ---------------------------------
   Returns: Number of matched coordinates
            0: Failure owing to mismatch

   Creates the coordinate arrays for fitting from the currently defined 
   zones.

   30.09.92 Original
   01.10.92 Ignores undefined atoms. Fix in finding mobile atoms. 
            Corrected use of CofGs   
   09.10.92 Removed incorrect resetting of CoorCount. This broke 
            multi-zone fitting
   17.07.95 Changed screen() to printf()
   18.07.95 Added initialisation of inserts in zones
            Added calls to FormatZone()
   31.07.95 Prints numbers of residues if mismatch
            Fixed bug in counting for weights array; was only counting
            reference structure, not mobile
   13.12.95 Added printing of zone info on number of residues mismatch
   31.05.96 Added test on B-values
   18.06.96 Replaced MODE_* with ZONE_MODE_*
            Replaced FindZone() with FindZonePDB()
   06.11.96 Negative BVal cutoff interpreted as > bval
   11.11.96 Checks actual value of gUseBVal
   18.11.96 Added gIgnoreMissing handling
   03.07.97 Added a break when the mobile atom has been found. This will
            speed things up and mean that we only find the first atom
            from mobile if multiple occupancies have not been specified
            correctly (i.e. the same atoms are not together); previously
            this would core dump.
   12.01.01 gMobPDB[] now an array
   01.02.01 Added strucnum parameter
   15.02.01 Changed ValidAtom() calls to atnam_raw
   20.02.01 gMobCofG now an array
   20.02.01 -999 for start or end of structure rather than -1
   18.03.08 added call to CentreOnZone() for setting centres of 
            geometry. By: CTP
*/
int CreateFitArrays(int strucnum)
{
   PDB     *ref_start = NULL,
           *ref_stop  = NULL,
           *mob_start = NULL,
           *mob_stop  = NULL,
           *p,
           *q,
           *r,
           *m;
   ZONE    *z;
   char    ref_insert,
           mob_insert;
   int     ref_resnum,
           mob_resnum,
           ref_nres,
           mob_nres,
           CoorCount = 0,
           Found,
           i,
           natom1, 
           natom2;
   
         
   if(gRefCoor==NULL || gMobCoor[strucnum]==NULL)
   {
      printf("   Error==> A coordinate array is undefined!\n");
      return(0);
   }

   /* Allocate memory for a weights array                               */
   for(p=gRefPDB, natom1=0; p!=NULL; NEXT(p))
      natom1++;
   for(p=gMobPDB[strucnum], natom2=0; p!=NULL; NEXT(p))
      natom2++;
   
   if(gWeights != NULL)
      free(gWeights);
   if((gWeights = (REAL *)malloc(MAX(natom1, natom2) * sizeof(REAL)))
      == NULL)
   {
      printf("   Error==> No memory for weights array!\n");
      return(0);
   }
   
   /* If no zones have been specified, create a single all atoms zone   */
   if(gZoneList[strucnum] == NULL)
   {
      INIT(gZoneList[strucnum],ZONE);
      gZoneList[strucnum]->chain1       = ' ';
      gZoneList[strucnum]->start1       = -999;
      gZoneList[strucnum]->startinsert1 = ' ';
      gZoneList[strucnum]->stop1        = -999;
      gZoneList[strucnum]->stopinsert1  = ' ';
      gZoneList[strucnum]->chain2       = ' ';
      gZoneList[strucnum]->start2       = -999;
      gZoneList[strucnum]->startinsert2 = ' ';
      gZoneList[strucnum]->stop2        = -999;
      gZoneList[strucnum]->stopinsert2  = ' ';
      gZoneList[strucnum]->mode         = gCurrentMode;
   }
   
   /* Step through each zone                                            */
   for(z=gZoneList[strucnum]; z!=NULL; NEXT(z))
   {
      /* Reference structure                                            */
      if(!FindZonePDB(gRefPDB, z->start1, z->startinsert1, 
                   z->stop1, z->stopinsert1, z->chain1, z->mode, 
                   &ref_start, &ref_stop))
      {
         char zone1[64],
              zone2[64];

         /* Check ranges have been found                                */
         printf("   Error==> Reference structure zone not found.\n");

         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         printf("      %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                         :"(Sequential numbering)"));

         return(0);
      }

      /* Mobile structure                                               */
      if(!FindZonePDB(gMobPDB[strucnum], z->start2, z->startinsert2, 
                   z->stop2, z->stopinsert2, z->chain2, z->mode, 
                   &mob_start, &mob_stop))
      {
         char zone1[64],
              zone2[64];

         /* Check ranges have been found                                */
         printf("   Error==> Mobile structure zone not found.\n");

         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         printf("      %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                         :"(Sequential numbering)"));
         return(0);
      }
      
      /* Check we have the same number of residues in each zone         */
      ref_nres     = 1;
      ref_resnum   = ref_start->resnum;
      ref_insert   = ref_start->insert[0];
      for(p=ref_start; p!=ref_stop; NEXT(p))
      {
         if(p->resnum != ref_resnum || p->insert[0] != ref_insert)
         {
            ref_nres++;
            ref_resnum = p->resnum;
            ref_insert = p->insert[0];
         }
      }

      mob_nres     = 1;
      mob_resnum   = mob_start->resnum;
      mob_insert   = mob_start->insert[0];
      for(p=mob_start; p!=mob_stop; NEXT(p))
      {
         if(p->resnum != mob_resnum || p->insert[0] != mob_insert)
         {
            mob_nres++;
            mob_resnum = p->resnum;
            mob_insert = p->insert[0];
         }
      }
      
      if(ref_nres != mob_nres)
      {
         char zone1[64],
              zone2[64];

         printf("   Error==> Number of residues in zone does not \
match.\n");
         /* Added 13.12.95                                              */
         FormatZone(zone1, z->chain1, 
                    z->start1, z->startinsert1, 
                    z->stop1,  z->stopinsert1);
         
         FormatZone(zone2, z->chain2, 
                    z->start2, z->startinsert2, 
                    z->stop2,  z->stopinsert2);
         
         printf("           %-16s with %-16s %s\n",
                zone1, zone2,
                ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                                         :"(Sequential numbering)"));

         /* Added 31.07.95                                              */
         printf("            Reference: %d, Mobile: %d\n",
                ref_nres, mob_nres);

         return(0);
      }


      /* Insert the atoms from this zone into the coordinate arrays     */

/*    Removed 09.10.92.....
//    CoorCount  = 0;   
*/

      ref_nres   = 0;
      ref_resnum = -999;
      ref_insert = ' ';

      for(p=ref_start; p!=ref_stop; NEXT(p))
      {
         if(p->resnum != ref_resnum || p->insert[0] != ref_insert)
         {
            /* Start of a new reference residue                         */
            ref_nres++;
            ref_resnum = p->resnum;
            ref_insert = p->insert[0];

            mob_nres   = 0;
            mob_resnum = -999;
            mob_insert = ' ';
         
            for(q=mob_start;
                q!=mob_stop; 
                NEXT(q))
            {
               if(q->resnum != mob_resnum || q->insert[0] != mob_insert)
               {
                  /* Start of a new mobile residue                      */
                  mob_nres++;
                  mob_resnum = q->resnum;
                  mob_insert = q->insert[0];
   
                  /* Consider the equivalent residues                   */
                  if(ref_nres == mob_nres)
                  {
                     /* p points to the start of a residue in reference,
                        q points to the start of a residue in mobile

                        Step through reference set
                     */
                     for(r=p;
                         r!=NULL && r->resnum==ref_resnum 
                                 && r->insert[0]==ref_insert;
                         NEXT(r))
                     {
                        /* 15.02.01 Changed from atnam to atnam_raw     */
                        if(ValidAtom(r->atnam_raw, ATOM_FITTING))
                        {
                           if(r->x == 9999.0 &&
                              r->y == 9999.0 &&
                              r->z == 9999.0)
                           {
                              if(!gQuiet)
                              {
                                 printf("   Warning==> Undefined atom in \
reference set ignored:\n");
                                 /* 28.02.01 Changed from r->atnam      */
                                 printf("      Residue: %4s %5d%2s, \
Atom: %4s\n",
                                        r->resnam,r->resnum,
                                        r->insert,r->atnam_raw);
                              }
                              
                              continue;
                           }

                           /* 31.05.96 Test the BVal
                              06.11.96 Modified to handle -ve 
                                       specifications
                              11.11.96 Checks for 1 or 2 as value
                           */
                           if(gUseBVal==1 || gUseBVal==2)
                           {
                              if(gBValue >= (REAL)0.0)
                              {
                                 if(r->bval > gBValue)
                                    continue;
                              }
                              else
                              {
                                 if(-(r->bval) > gBValue)
                                    continue;
                              }
                           }
                           
                           /* Find this atom in the mobile set          */
                           Found = FALSE;
                           for(m=q;
                               m!=NULL && m->resnum==mob_resnum 
                                       && m->insert[0]==mob_insert;
                               NEXT(m))
                           {
                              /* 28.02.01 Changed from ->atnam          */
                              if(!strcmp(r->atnam_raw,m->atnam_raw))
                              {
                                 Found = TRUE;
                                 
                                 if(m->x == 9999.0 &&
                                    m->y == 9999.0 &&
                                    m->z == 9999.0)
                                 {
                                    if(!gQuiet)
                                    {
                                       printf("   Warning==> Undefined \
atom in mobile set ignored:\n");
                                       /* 28.02.01 Changed from r->atnam*/
                                       printf("      Residue: %4s \
%5d%2s, Atom: %4s\n",
                                              m->resnam,m->resnum,
                                              m->insert,m->atnam_raw);
                                    }
                                    
                                    continue;
                                 }


                                 /* 31.05.96 Test the BVal
                                    06.11.96 Modified to handle -ve 
                                             specifications
                                    11.11.96 Checks for 1 or 3 as value
                                 */
                                 if(gUseBVal==1 || gUseBVal==3)
                                 {
                                    if(gBValue >= (REAL)0.0)
                                    {
                                       if(m->bval > gBValue)
                                          continue;
                                    }
                                    else
                                    {
                                       if(-(m->bval) > gBValue)
                                          continue;
                                    }
                                 }
                                 
                                 /* Copy the coordinates                */
                                 gRefCoor[CoorCount].x = r->x;
                                 gRefCoor[CoorCount].y = r->y;
                                 gRefCoor[CoorCount].z = r->z;
                                 gMobCoor[strucnum][CoorCount].x = m->x;
                                 gMobCoor[strucnum][CoorCount].y = m->y;
                                 gMobCoor[strucnum][CoorCount].z = m->z;
                                 gWeights[CoorCount] = 
                                    (m->bval + r->bval)/(REAL)2.0;
                                 CoorCount++;
                                 break;             /* 03.07.97 ACRM    */
                              }
                           }  
                           if(!Found)
                           {
                              if(gIgnoreMissing)
                              {
                                 if(!gQuiet)
                                 {
                                    /* 28.02.01 Changed from r->atnam   */
                                    printf("   Warning==> Ignored \
reference atom %4s not found in mobile.\n",
                                           r->atnam_raw);
                                    printf("      Reference %4s %5d%c \
Mobile %4s %5d%c\n",
                                           p->resnam,p->resnum,
                                           p->insert[0],
                                           q->resnam,q->resnum,
                                           q->insert[0]);
                                 }
                              }
                              else
                              {
                                 printf("   Error==> Atoms do not match \
in residue:\n");
                                 printf("      Reference %4s %5d%c \
Mobile %4s %5d%c\n",
                                        p->resnam,p->resnum,p->insert[0],
                                        q->resnam,q->resnum,q->insert[0]);
                                 /* 28.02.01 Changed from r->atnam      */
                                 printf("      Unable to find \
reference atom %4s in mobile.\n",
                                        r->atnam_raw);
                                 return(0);
                              }
                           }
                        }
                     }  /* End of loop through reference residue        */
                     break;
                  }  /* End of if() equivalent residues                 */
               }  /* End of start-of-new-mobile-residue                 */
            }  /* End of loop through mobile set                        */
         }  /* End of start-of-new-reference-residue                    */
      }  /* End if loop through reference set                           */
   }  /* End of loop through zones                                      */
   
   if(CoorCount == 0)
   {
      printf("   Error==> No atoms in specified zones\n");
   }
   else
   {
     /* Calculate the centre of geometry                               */
     /* 18.03.08 Call to CentreOnZone added By: CTP                    */
     gMobCofG[strucnum].x = gMobCofG[strucnum].y 
                          = gMobCofG[strucnum].z = 0.0;
     gRefCofG.x = gRefCofG.y = gRefCofG.z = 0.0;

     if(gCZoneList[strucnum] == NULL)
     {
        for(i=0; i<CoorCount; i++)
        {
           gMobCofG[strucnum].x  += gMobCoor[strucnum][i].x;
           gMobCofG[strucnum].y  += gMobCoor[strucnum][i].y;
           gMobCofG[strucnum].z  += gMobCoor[strucnum][i].z;
           gRefCofG.x            += gRefCoor[i].x;
           gRefCofG.y            += gRefCoor[i].y;
           gRefCofG.z            += gRefCoor[i].z;
        }
        
        gMobCofG[strucnum].x /= CoorCount;
        gMobCofG[strucnum].y /= CoorCount;
        gMobCofG[strucnum].z /= CoorCount;
        gRefCofG.x           /= CoorCount;
        gRefCofG.y           /= CoorCount;
        gRefCofG.z           /= CoorCount;
     }
     else
     {
        if(!CentreOnZone(strucnum))
        {
           printf("   Error==> No centre residues set.\n");
           return(0);
        }
     }
     
#ifdef DEBUG
      fprintf(stderr,"Before fitting\n%d coordinates.\n", CoorCount);
      fprintf(stderr,"Ref CofG: %8.3f %8.3f %8.3f\n",
              gRefCofG.x, gRefCofG.y, gRefCofG.z);
      fprintf(stderr,"Mob CofG: %8.3f %8.3f %8.3f\n",
              gMobCofG[strucnum].x, 
              gMobCofG[strucnum].y, 
              gMobCofG[strucnum].z);
#endif
      
      /* Move coordinate arrays to the origin                           */
      for(i=0; i<CoorCount; i++)
      {
         gRefCoor[i].x -= gRefCofG.x;
         gRefCoor[i].y -= gRefCofG.y;
         gRefCoor[i].z -= gRefCofG.z;
         gMobCoor[strucnum][i].x -= gMobCofG[strucnum].x;
         gMobCoor[strucnum][i].y -= gMobCofG[strucnum].y;
         gMobCoor[strucnum][i].z -= gMobCofG[strucnum].z;
      }
   }
   
   return(CoorCount);
}


/************************************************************************/
/*>int CentreOnZone(int strucnum)
   ------------------------------
   Returns: Number of matched coordinates
            0: Failure owing to mismatch

   Sets the centres of geometry (gRefCofG and gRefCofG[structnum]) to 
   correspond to centres of user-defined zones.

   18.03.08 Original based on CreateFitArrays() By: CTP
*/
int CentreOnZone(int strucnum)
{
   PDB     *ref_start = NULL,
           *ref_stop  = NULL,
           *mob_start = NULL,
           *mob_stop  = NULL,
           *p,
           *q,
           *r,
           *m;
   ZONE    *z;
   char    ref_insert,
           mob_insert;
   int     ref_resnum,
           mob_resnum,
           ref_nres,
           mob_nres,
           CoorCount = 0;
   BOOL    Found;

   VEC3F ref_CofG;
   VEC3F mob_CofG;

   ref_CofG.x = ref_CofG.y = ref_CofG.z = 0.0;
   mob_CofG.x = mob_CofG.y = mob_CofG.z = 0.0;
         
   if(gRefCoor==NULL || gMobCoor[strucnum]==NULL)
   {
      printf("   Error==> A coordinate array is undefined!\n");
      return(0);
   }
   
   /* If no zones have been specified, return.                          */
   if(gCZoneList[strucnum] == NULL)
   {
      printf("   CentreOnZone: No Zones Specified.\n");
      return(0);
   }
   
   /* Step through each zone                                            */
   for(z=gCZoneList[strucnum]; z!=NULL; NEXT(z))
   {
      /* Reference structure                                            */
      if(!FindZonePDB(gRefPDB, z->start1, z->startinsert1, 
                      z->stop1, z->stopinsert1, z->chain1, z->mode, 
                      &ref_start, &ref_stop))
      {

         /* Check ranges have been found                                */
         printf("   Error==> Reference centre residue not found.\n");
         return(0);
      }
      
      /* Mobile structure                                               */
      if(!FindZonePDB(gMobPDB[strucnum], z->start2, z->startinsert2, 
                      z->stop2, z->stopinsert2, z->chain2, z->mode, 
                      &mob_start, &mob_stop))
      {
         /* Check ranges have been found                                */
         printf("   Error==> Mobile centre residue not found.\n");
         return(0);
      }
     

      /* Find atom coordinates                                          */
      ref_nres   = 0;
      ref_resnum = -999;
      ref_insert = ' ';

      for(p=ref_start; p!=ref_stop; NEXT(p))
      {
         if(p->resnum != ref_resnum || p->insert[0] != ref_insert)
         {
            /* Start of a new reference residue                         */
            ref_nres++;
            ref_resnum = p->resnum;
            ref_insert = p->insert[0];

            mob_nres   = 0;
            mob_resnum = -999;
            mob_insert = ' ';
         
            for(q=mob_start;
                q!=mob_stop; 
                NEXT(q))
            {
               if(q->resnum != mob_resnum || q->insert[0] != mob_insert)
               {
                  /* Start of a new mobile residue                      */
                  mob_nres++;
                  mob_resnum = q->resnum;
                  mob_insert = q->insert[0];
   
                  /* Consider the equivalent residues                   */
                  if(ref_nres == mob_nres)
                  {
                     /* p points to the start of a residue in reference,
                        q points to the start of a residue in mobile

                        Step through reference set
                     */
                     for(r=p;
                         r!=NULL && r->resnum==ref_resnum 
                                 && r->insert[0]==ref_insert;
                         NEXT(r))
                     {
                        /* 15.02.01 Changed from atnam to atnam_raw     */
                        if(ValidAtom(r->atnam_raw, ATOM_FITTING))
                        {
                           if(r->x == 9999.0 &&
                              r->y == 9999.0 &&
                              r->z == 9999.0)
                           {
                              if(!gQuiet)
                              {
                                 printf("   Warning==> Undefined atom in \
reference set ignored:\n");
                                 /* 28.02.01 Changed from r->atnam      */
                                 printf("      Residue: %4s %5d%2s, \
Atom: %4s\n",
                                        r->resnam,r->resnum,
                                        r->insert,r->atnam_raw);
                              }
                              
                              continue;
                           }
                           
                           /* Find this atom in the mobile set          */
                           Found = FALSE;
                           for(m=q;
                               m!=NULL && m->resnum==mob_resnum 
                                       && m->insert[0]==mob_insert;
                               NEXT(m))
                           {
                              /* 28.02.01 Changed from ->atnam          */
                              if(!strcmp(r->atnam_raw,m->atnam_raw))
                              {
                                 Found = TRUE;
                                 
                                 if(m->x == 9999.0 &&
                                    m->y == 9999.0 &&
                                    m->z == 9999.0)
                                 {
                                    if(!gQuiet)
                                    {
                                       printf("   Warning==> Undefined \
atom in mobile set ignored:\n");
                                       /* 28.02.01 Changed from r->atnam*/
                                       printf("      Residue: %4s \
%5d%2s, Atom: %4s\n",
                                              m->resnam,m->resnum,
                                              m->insert,m->atnam_raw);
                                    }
                                    
                                    continue;
                                 }

                                 /* Sum of Coordinates                  */
                                 ref_CofG.x += r->x;
                                 ref_CofG.y += r->y;
                                 ref_CofG.z += r->z;
                                 mob_CofG.x += m->x;
                                 mob_CofG.y += m->y;
                                 mob_CofG.z += m->z;
                                 CoorCount++;
                                 break;
                              }
                           }  
                           if(!Found)
                           {
                              if(gIgnoreMissing)
                              {
                                 if(!gQuiet)
                                 {
                                    /* 28.02.01 Changed from r->atnam   */
                                    printf("   Warning==> Ignored \
reference atom %4s not found in mobile.\n",
                                           r->atnam_raw);
                                    printf("      Reference %4s %5d%c \
Mobile %4s %5d%c\n",
                                           p->resnam,p->resnum,
                                           p->insert[0],
                                           q->resnam,q->resnum,
                                           q->insert[0]);
                                 }
                              }
                              else
                              {
                                 printf("   Error==> Atoms do not match \
in residue:\n");
                                 printf("      Reference %4s %5d%c \
Mobile %4s %5d%c\n",
                                        p->resnam,p->resnum,p->insert[0],
                                        q->resnam,q->resnum,q->insert[0]);
                                 /* 28.02.01 Changed from r->atnam      */
                                 printf("      Unable to find \
reference atom %4s in mobile.\n",
                                        r->atnam_raw);
                                 return(0);
                              }
                           }
                        }
                     }  /* End of loop through reference residue        */
                     break;
                  }  /* End of if() equivalent residues                 */
               }  /* End of start-of-new-mobile-residue                 */
            }  /* End of loop through mobile set                        */
         }  /* End of start-of-new-reference-residue                    */
      }  /* End if loop through reference set                           */
   }  /* End of loop through zones                                      */
   
   if(CoorCount == 0)
   {
      printf("   Error==> No atoms in specified zones\n");
   }
   else
   {
      /* Calculate the centre of geometry                               */
      gRefCofG.x           = ref_CofG.x / CoorCount;
      gRefCofG.y           = ref_CofG.y / CoorCount;
      gRefCofG.z           = ref_CofG.z / CoorCount;
      gMobCofG[strucnum].x = mob_CofG.x / CoorCount;
      gMobCofG[strucnum].y = mob_CofG.y / CoorCount;
      gMobCofG[strucnum].z = mob_CofG.z / CoorCount;
   }
   
   return(CoorCount);
}


/************************************************************************/
/*>void SetSymmetricalAtomPAirs()
   ------------------------------
   Sets the atom pairs used when auto-matching symmetrical atoms.
   (eg CD1 - CD2 and CE1 - CE2 in Tyr)

   04.06.08 Original  By: CTP
*/
void SetSymmetricalAtomPAirs(void)
{
   /* Charged                                                           */
   strcpy(gSymType[0][0],"ARG ");
   strcpy(gSymType[0][1]," NH1");
   strcpy(gSymType[0][2]," NH2");
   
   strcpy(gSymType[1][0],"ASP ");
   strcpy(gSymType[1][1]," OD1");
   strcpy(gSymType[1][2]," OD2");
   
   strcpy(gSymType[2][0],"GLU ");
   strcpy(gSymType[2][1]," OE1");
   strcpy(gSymType[2][2]," OE2");
   
   /* Aromatic                                                          */
   strcpy(gSymType[3][0],"PHE ");
   strcpy(gSymType[3][1]," CD1");
   strcpy(gSymType[3][2]," CD2");
   strcpy(gSymType[4][0],"PHE ");
   strcpy(gSymType[4][1]," CE1");
   strcpy(gSymType[4][2]," CE2");
   
   strcpy(gSymType[5][0],"TYR ");
   strcpy(gSymType[5][1]," CD1");
   strcpy(gSymType[5][2]," CD2");
   strcpy(gSymType[6][0],"TYR ");
   strcpy(gSymType[6][1]," CE1");
   strcpy(gSymType[6][2]," CE2");
   
   /* Amide Nitrogen and Oxygen                                         */
   strcpy(gSymType[7][0],"ASN ");
   strcpy(gSymType[7][1]," OD1");
   strcpy(gSymType[7][2]," ND2");
   
   strcpy(gSymType[8][0],"GLN ");
   strcpy(gSymType[8][1]," OE1");
   strcpy(gSymType[8][2]," NE2");
   
   /* Prochiral Methyls                                                 */
   strcpy(gSymType[9][0],"VAL ");
   strcpy(gSymType[9][1]," CG1");
   strcpy(gSymType[9][2]," CG2");
   
   strcpy(gSymType[10][0],"LEU ");
   strcpy(gSymType[10][1]," CD1");
   strcpy(gSymType[10][2]," CD2");
   
   /* Default Match Settings                                            */
   gSymType[0][3][0]  = TRUE;   /* Arg */
   gSymType[1][3][0]  = TRUE;   /* Asp */
   gSymType[2][3][0]  = TRUE;   /* Glu */
   gSymType[3][3][0]  = TRUE;   /* Phe */
   gSymType[4][3][0]  = TRUE;   /* Phe */
   gSymType[5][3][0]  = TRUE;   /* Tyr */
   gSymType[6][3][0]  = TRUE;   /* Tyr */
   gSymType[7][3][0]  = FALSE;  /* Asn */
   gSymType[8][3][0]  = FALSE;  /* Gln */
   gSymType[9][3][0]  = FALSE;  /* Val */
   gSymType[10][3][0] = FALSE;  /* Leu */

   return;
}

/************************************************************************/
/*>void ApplyMatrixCOOR(COOR *incoords, REAL matrix[3][3],int  ncoor)
   ------------------------------------------------------------------
   I/O:   COOR *incoords     Coordinate Array
   Input: REAL matrix[3][3]  Matrix to apply
          int  ncoor         Number of coordinates.

   Apply a rotation matrix to a coordinate array. 
   Based on ApplyMatrixPDB() by ACRM.

   18.08.08 Original By: CTP
*/
void ApplyMatrixCOOR(COOR *incoords,
                     REAL matrix[3][3],
                     int  ncoor)
{
  int   i;
  VEC3F outcoords;

  for(i=0; i<ncoor; i++)
  {
     if(incoords[i].x != 9999.0 && incoords[i].y != 9999.0 && 
        incoords[i].z != 9999.0)
     {
        MatMult3_33(incoords[i],matrix,&outcoords);
        
        incoords[i].x = outcoords.x;
        incoords[i].y = outcoords.y;
        incoords[i].z = outcoords.z;
     }
  }
  
  return;
}


/************************************************************************/
/*>void CalculateRotationMatrix(REAL RotAngle, REAL Matrix[3][3])
   --------------------------------------------------------------
   Calculate rotation matrix around Z axis for input angle, RotAngle, in
   degrees.

   18.08.08 Original By: CTP
   16.02.09 Rewritten as wrapper for bioplib function CreateRotMat()
*/
void CalculateRotationMatrix(REAL RotAngle, REAL Matrix[3][3])
{
   /* Convert to radians                                                */
   REAL psi = RotAngle * -1.0 * (PI/180.0); 
   
   /* Rotation Matrix                                                   */
   CreateRotMat('z', psi, Matrix);
   return;
}


/************************************************************************/
/*>REAL FitSingleStructure(int strucnum, BOOL single_iteration)
   ------------------------------------------------------------
   Fits single mobile structure to the reference structure and returns
   the RMSD. Called during all vs all comparisons, when setting the order
   for fitting structures in order and when fitting structures in order.

   Function returns -1.0 if error.

   20.10.08 Original based on FitStructures() By: CTP
*/
REAL FitSingleStructure(int strucnum, BOOL single_iteration)
{
   ZONE  *z1,
         *z2;
   int   atmnum,
         NCoor,
         niter;
   REAL  rmstot,
         rmsprev = (-100.0),
         deltaRMS,
         rmscurr = -1.0;
   BOOL  final = FALSE;
   
   
   gFitted = FALSE;
   
   if(!gRefFilename[0])
   {
      printf("   Error==> Reference structure undefined.\n");
      return(-1.0);
   }
   if(!gMobFilename[0][0])
   {
      printf("   Error==> Mobile structure undefined.\n");
      return(-1.0);
   }
   
   /* First copy the zones for display to match those for fitting       */
   if(gRZoneList[strucnum] != NULL)
   {
      FREELIST(gRZoneList[strucnum],ZONE);
      gRZoneList[strucnum] = NULL;
   }
   for(z1=gZoneList[strucnum]; z1!=NULL; NEXT(z1))
   {
      /* Allocate an entry in RMS zone list                             */
      if(gRZoneList[strucnum])
      {
         /* Move to end of zone list                                    */
         z2=gRZoneList[strucnum];
         LAST(z2);
         ALLOCNEXT(z2,ZONE);
      }
      else
      {
         INIT(gRZoneList[strucnum],ZONE);
         z2 = gRZoneList[strucnum];
      }
      
      if(z2==NULL)
      {
         printf("   Error==> No memory for RMS zone!\n");
      }
      else
      {
         /* Add this zone to the RMS zone list                          */
         z2->chain1       = z1->chain1;
         z2->start1       = z1->start1;
         z2->startinsert1 = z1->startinsert1;
         z2->stop1        = z1->stop1;
         z2->stopinsert1  = z1->stopinsert1;
         z2->chain2       = z1->chain2;
         z2->start2       = z1->start2;
         z2->startinsert2 = z1->startinsert2;
         z2->stop2        = z1->stop2;
         z2->stopinsert2  = z1->stopinsert2;
         z2->mode         = z1->mode;
      }
   }
   
   /* Now copy the atoms for RMS calculation                            */
   gNOTRMSAtoms = gNOTFitAtoms;
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      strcpy(gRMSAtoms[atmnum],gFitAtoms[atmnum]);
   
   if(gMultiCount > 1)
   {
      /* Keep looping counting the iterations                           */
      for(niter=0; ; niter++)
      {
         /* printf ("   Multi-structure fit iteration %d\n", niter);    */
         rmstot = (REAL)0.0;
         
         /*** Fit single structure                                    ***/
         
         /* Set up arrays for fitting                                   */
         if((NCoor=CreateFitArrays(strucnum))!=0)
         {
            /* Reset the convergence criterion                          */
            CheckForConvergence(0, strucnum);
            
            /* Perform the fit                                          */
            if(DoFitting(NCoor, strucnum))
            {
               while(gIterate && !final)
               {
                  if((NCoor = UpdateFitArrays(strucnum))!=0)
                  {
                     if(!DoFitting(NCoor, strucnum))
                        return(-1.0);
                     if(CheckForConvergence(NCoor, strucnum))
                        break;
                  }
                  else
                  {
                     break;
                  }
               }
               
               /* Find RMS - Do not update during final iteration       */
               /* rmstot += ShowRMS(FALSE,NULL,strucnum,!final,FALSE);  */
               rmscurr = ShowRMS(FALSE,NULL,strucnum,
                                 !final && !single_iteration,FALSE);
               rmstot += rmscurr;
               if(gIterate && !gQuiet)
               {
/***
//                  printf("      (Over %d equivalenced CA-atoms)\n",
//                  NCoor);
***/
               }
            }
         }
         
         deltaRMS = (rmstot - rmsprev);
         rmsprev = rmstot;
         
         /* If we've converged or done too many iterations, do a final  */
         /* iteration then break out.                                   */
         if(final) break;
         
         if((ABS(deltaRMS) < MULTI_ITER_STOP) ||
            (niter > MAXMULTIITER))
            final = TRUE;
         
         /* Break out if single iteration                               */
         if(single_iteration) break;
      }
   }
   
   return(rmscurr);
}


/************************************************************************/
/*>void FitStructuresInOrder(REAL sortorder[][2])
   ----------------------------------------------
   Fit structures in order based on 2D array, sortorder. The first element
   of sortorder is the structure number and the second is a score used for
   sorting (eg RMSD).

   20.10.08 Original based on FitStructures() By: CTP
*/
void FitStructuresInOrder(REAL sortorder[][2])
{
   ZONE  *z1,
         *z2;
   int   atmnum,
         NCoor,
         strucnum,
         niter,
         i;
   REAL  rmstot,
         rmsprev = (-100.0),
         deltaRMS;
   BOOL  final = FALSE;
   
   gFitted = FALSE;

   if(!gRefFilename[0])
   {
      printf("   Error==> Reference structure undefined.\n");
      return;
   }
   if(!gMobFilename[0][0])
   {
      printf("   Error==> Mobile structure undefined.\n");
      return;
   }
   if(gMultiCount <= 1)
   {
      printf("   Error==> ");
      printf("ORDERFIT can only be used with multiple structures.\n");
      return;
   }
   
   if(!gQuiet)
   {
      printf("   Fitting structures...\n");
   }

   /* First copy the zones for display to match those for fitting       */
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      if(gRZoneList[strucnum] != NULL)
      {
         FREELIST(gRZoneList[strucnum],ZONE);
         gRZoneList[strucnum] = NULL;
      }
      for(z1=gZoneList[strucnum]; z1!=NULL; NEXT(z1))
      {
         /* Allocate an entry in RMS zone list                          */
         if(gRZoneList[strucnum])
         {
            /* Move to end of zone list                                 */
            z2=gRZoneList[strucnum];
            LAST(z2);
            ALLOCNEXT(z2,ZONE);
         }
         else
         {
            INIT(gRZoneList[strucnum],ZONE);
            z2 = gRZoneList[strucnum];
         }
         
         if(z2==NULL)
         {
            printf("   Error==> No memory for RMS zone!\n");
         }
         else
         {
            /* Add this zone to the RMS zone list                       */
            z2->chain1       = z1->chain1;
            z2->start1       = z1->start1;
            z2->startinsert1 = z1->startinsert1;
            z2->stop1        = z1->stop1;
            z2->stopinsert1  = z1->stopinsert1;
            z2->chain2       = z1->chain2;
            z2->start2       = z1->start2;
            z2->startinsert2 = z1->startinsert2;
            z2->stop2        = z1->stop2;
            z2->stopinsert2  = z1->stopinsert2;
            z2->mode         = z1->mode;
         }
      }
   }

   /* Now copy the atoms for RMS calculation                            */
   gNOTRMSAtoms = gNOTFitAtoms;
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      strcpy(gRMSAtoms[atmnum],gFitAtoms[atmnum]);
      
   /* Keep looping counting the iterations                              */
   for(niter=0; ; niter++)
   {
      printf ("   Multi-structure fit iteration %d\n", niter);
      
      rmstot = (REAL)0.0;
      
      /* Loop through the structures we are fitting                     */
      /* for(strucnum=0; strucnum<gMultiCount; strucnum++) */
      for(i=0; i<gMultiCount; i++)
      {
         
         /* Set strucnum */
         strucnum = (int)sortorder[i][0];
         printf("   Fitting structure %d\n",strucnum + 1);
         
         /* Set up arrays for fitting                                   */
         if((NCoor=CreateFitArrays(strucnum))!=0)
         {
            /* Reset the convergence criterion                          */
            CheckForConvergence(0, strucnum);
            
            /* Perform the fit                                          */
            if(DoFitting(NCoor, strucnum))
            {
               while(gIterate && !final)
               {
                  if((NCoor = UpdateFitArrays(strucnum))!=0)
                  {
                     if(!DoFitting(NCoor, strucnum))
                        return;
                     if(CheckForConvergence(NCoor, strucnum))
                        break;
                  }
                  else
                  {
                     break;
                  }
               }
               
               /* Find RMS - Do not update during final iteration       */
               rmstot += ShowRMS(FALSE,NULL,strucnum,!final,FALSE);
               if(gIterate && !gQuiet)
               {
                  printf("      (Over %d equivalenced CA-atoms)\n",
                         NCoor);
               }
            }
         }
      }
      deltaRMS = (rmstot - rmsprev);
      rmsprev = rmstot;
      
      /* If we've converged or done too many iterations, do a final
         iteration then break out.                                 
      */
      if(final) break;
      
      if((ABS(deltaRMS) < MULTI_ITER_STOP) ||
         (niter > MAXMULTIITER))
         final = TRUE;
   }

   return;
}




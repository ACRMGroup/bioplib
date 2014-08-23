/*************************************************************************

   Program:    ProFit
   File:       NWAlign.c
   
   Version:    V3.1
   Date:       31.03.09
   Function:   Protein Fitting program. 
   
   Copyright:  SciTech Software 1992-2009
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  25.09.92 Original
   V0.5  08.10.93 Various tidying for Unix & chaned for booklib 
   V0.6  05.01.94 Modified MDMFILE for Unix getenv()
   V0.7  24.11.94 The DATAENV environment variable is now handled by code
                  in bioplib/align.c/ReadMDM; Checks the return from this
                  Fixed bug in multi-zone align
   V0.8  17.07.95 Replaced screen() stuff with printf()
                  Only allowed on single chains
   V1.0  18.07.95 Insert codes now work.
                  First official release (at last!).
   V1.1  20.07.95 Skipped
   V1.2  22.07.95 Added GAPPEN command making gap penalty global variable
   V1.3  31.07.95 Skipped
   V1.4  14.08.95 Skipped
   V1.5  21.08.95 Fixed bug in mapping alignment to zones. Also bug in
                  Bioplib align() routine
   V1.5b 15.11.95 Now also prints a score normalised by the length of the
                  shorter sequence.
   V1.6  20.11.95 Added ReadAlignment() code
   V1.6b 22.11.95 Added check in SetNWZones() for a deletion in both
                  sequences
   V1.6c 13.12.95 The check added in 1.6b wasn't working. Fixed!
   V1.6g 18.06.96 Changed MODE_* to ZONE_MODE_*
   V1.7  23.07.96 Skipped
   V1.7g 06.05.98 Rewrite of SetNWZones()
   V1.8  07.05.98 Skipped for release
   V2.0  01.03.01 Additions for iterative zone updating
   V2.1  28.03.01 Parameter for ITERATE and added CENTRE command
   V2.2  20.12.01 Skipped for release
   V2.3  01.12.04 Fixed bugs in removing double deletions from alignment
                  with multiple structures
   V2.4  03.06.05 Skipped for release
   V2.5  07.06.05 Skipped for release
   V2.6  23.04.08 Added VerifySequence(). By: CTP
   V3.0  06.11.08 Added multi-chain alignment and ability to print fitted
                  zones as an alignment.
   V3.0  25.11.08 Changed output format for printing fitted zones as 
                  alignment.
   V3.0  13.01.08 Included output of fitting zones as PIR alignment.
   V3.1  31.03.09 Skipped for release.   

*************************************************************************/
/* Includes
*/
#include "ProFit.h"

/************************************************************************/
/*>NWAlign(int strucnum)
   ---------------------
   28.09.92 Framework
   09.10.92 Original
   05.01.94 Modified to get data directory name from environment variable
            under Unix
   17.07.95 Replaced screen() with printf()
            Check only one chain in each structure.
   22.07.95 Made gap penalty a variable
   15.11.95 Also prints a score normalised by the length of the shorter 
            sequence.
   01.02.01 Added strucnum parameter
   16.07.08 Changed alignment function from align() to affinealign()
            to allow for inclusion of gap extension penalty, gGapPenExt.
            By: CTP
*/
void NWAlign(int strucnum)
{
   static   int FirstCall = TRUE;
   int      ref_len,
            mob_len,
            align_len,
            offset,
            score,
            i,    j,
            ai,   aj;
   char     *ref_align = NULL,
            *mob_align = NULL;
         
   printf("   Performing N&W alignment...\n");

   if(FirstCall)
   {
      if(!ReadMDM(MDMFILE))
      {
         printf("   Error==> Unable to read mutation data matrix\n");
         return;
      }
      
      FirstCall = FALSE;
   }
   
   /* Make checks that structures read                                  */
   if(gRefSeq==NULL || gMobSeq[strucnum]==NULL)
   {
      printf("   Error==> Structures have not been read!\n");
      return;
   }

   /* Check for numbers of chains                                       */
   if((countchar(gRefSeq,'*') > 0) || 
      (countchar(gMobSeq[strucnum],'*') > 0))
   {
      printf("   Error==> Structures must have only one chain for \
alignment\n");
      return;
   }

   /* Find sequence lengths                                             */
   ref_len = strlen(gRefSeq);
   mob_len = strlen(gMobSeq[strucnum]);
   
   /* Allocate memory for alignment sequences                           */
   if((ref_align = (char *)malloc((ref_len+mob_len)*sizeof(char)))==NULL)
   {
      printf("   Warning==> No memory for alignment!\n");
      return;
   }
   if((mob_align = (char *)malloc((ref_len+mob_len)*sizeof(char)))==NULL)
   {
      printf("   Warning==> No memory for alignment!\n");
      free(ref_align);
      return;
   }
   
   /* Perform the alignment                                             */
   /* Alignment function modified 16.07.08 */
   /*
   score = align(gRefSeq, ref_len, gMobSeq[strucnum], mob_len, FALSE, 
                 FALSE, gGapPen, ref_align, mob_align, &align_len);
   */

   score = affinealign(gRefSeq, ref_len, gMobSeq[strucnum], mob_len, 
                       FALSE, FALSE, gGapPen, gGapPenExt, 
                       ref_align, mob_align, &align_len);


   if(!score)
   {
      printf("   Error==> Unable to perform alignment!\n");
      return;
   }
   
   /* Display the fitted sequences                                      */
   offset = 0;
   printf("   ");
   for(i=0,ai=0,aj=0; ai<align_len; ai++)    /* Prints ref sequence     */
   {
      char  buffer[8];
      
      if(++i>60)  /* If printed 60 chars, print equiv section of mob seq*/
      {
         i=1;
         printf("\n   ");
         for(j=offset; j<60+offset; j++)
         {
            sprintf(buffer,"%c",mob_align[j]);
            printf(buffer);
         }
         printf("\n\n   ");
         offset += 60;
      }
      printf("%c",ref_align[ai]);
   }
   printf("\n   ");

   for(j=offset; j<align_len; j++)          /* Print remains of mob seq */
   {
      printf("%c",mob_align[j]);
   }
   printf("\n\n   ");
   
   printf("Score: %d Normalised score: %.2f\n",
          score,
          (REAL)score/(REAL)(MIN(ref_len,mob_len)));
   
   /* Clear any current fitting zones                                   */
   SetFitZone("CLEAR", strucnum);
   
   /* Now set zones based on alignment                                  */
   SetNWZones(ref_align, mob_align, align_len, NULL, NULL, strucnum);
   
   /* Free allocated memory                                             */
   free(ref_align);
   free(mob_align);
   
   return;
}



/************************************************************************/
/*>void ReadAlignment(char *alnfile)
   ---------------------------------
   Read the first two sequences out of an alignment file in PIR format
   and set up zones based on the alignment.

   20.11.95 Original    By: ACRM
   01.02.01 Modified to cope with multiple structures.
   01.12.04 Fixed problem with multiple structures - didn't actually
            handle removing double deletions properly since the
            sequences were modified in-place
   23.04.08 added check of read sequences against stored 
            sequences. By: CTP
   03.12.08 Adapted to read multi-chain sequences.
   13.01.09 Tidied code.

*/
void ReadAlignment(char *alnfile)
{
   FILE *fp;
   char *seqa[MAXCHAIN],
        *seqb[MAXCHAIN],
        *tseqa = NULL,
        *ref_string = NULL,
        *mob_string = NULL;
   BOOL punct, error;
   int  i,
        nchain,
        strucnum = 0,
        chainlength = 0;
   
   /* Open the PIR alignment file for reading                           */
   if((fp=fopen(alnfile,"r"))==NULL)
   {
      printf("   Error==> Unable to read alignment file (%s)\n", alnfile);
      return;
   }
   
   /* Read the first sequence from the file                             */
   nchain = ReadPIR(fp, TRUE, seqa, MAXCHAIN, NULL, &punct, &error);

   /* No sequence found */
   if(nchain == 0)
   {
      printf("   Error==> No sequence read from alignment file.\n");
      fclose(fp);
      return;
   }

   /* Make reference string including whole sequence */
   for(i=0;i<nchain;i++) chainlength += strlen(seqa[i]);
   ref_string = (char *)malloc((nchain + chainlength)*sizeof(char));
   strcpy(ref_string,seqa[0]);
   free(seqa[0]);
   for(i=1;i<nchain;i++) 
   {
      strcat(ref_string,"*");
      strcat(ref_string,seqa[i]);
      free(seqa[i]);
   }

   /* Verify file sequence against Reference sequence          */
   if(!VerifySequence(ref_string, gRefSeq))
   {
      printf("   Error==> Reference sequence doesn't match.\n");
      fclose(fp);
      return;
   }

   /* Eliminate chain breaks */
   ChainBreakToGap(ref_string);
   

   /* 01.12.04 Allocate a buffer to store a working copy                */
   if((tseqa=(char *)malloc((1+strlen(ref_string))*sizeof(char)))==NULL)
   {
      printf("   Error==> No memory for working copy of sequence.\n");
      fclose(fp);
      return;
   }

   /* Read the second sequence from the file                            */
   while((nchain = ReadPIR(fp, TRUE, seqb, MAXCHAIN, NULL, 
                           &punct, &error)))
   {
      /* Check sequence found */
      if(nchain == 0)
      {
         printf("   Error==> No sequence read from alignment file.\n");
         fclose(fp);
         return;
      }

      
      /* Make mobile string including whole sequence */
      for(i=0;i<nchain;i++) chainlength += strlen(seqb[i]);
      mob_string = (char *)malloc((nchain + chainlength)*sizeof(char));
      strcpy(mob_string,seqb[0]);
      free(seqb[0]);
      for(i=1;i<nchain;i++)
      {
         strcat(mob_string,"*");
         strcat(mob_string,seqb[i]);
         free(seqb[i]);
      }

      /* 23.04.08 Verify file sequence against Mobile sequence          */
      /*if(!VerifySequence(seqb[0], gMobSeq[strucnum]))*/
      if(!VerifySequence(mob_string, gMobSeq[strucnum]))
      {
         printf("   Error==> Mobile sequence: %d doesn't match.\n",
                strucnum+1);
         fclose(fp);
         return;
      }
      
      /* Eliminate chain breaks */
      ChainBreakToGap(mob_string);
      
      /* 01.12.04 Make working copy of sequence A                       */
      /*strcpy(tseqa,seqa[0]);*/
      strcpy(tseqa,ref_string);
      
      /* Clear any current fitting zones                                */
      SetFitZone("CLEAR", strucnum);
      
      /* Remove any deletions which appear in both sequences            */
      /* 01.12.04 Changed to tseqa rather than seqa[0]                  */
      /*if(!RemoveDoubleDeletions(tseqa, seqb[0]))*/
      if(!RemoveDoubleDeletions(tseqa, mob_string))
      {
         printf("   Warning==> No memory to remove double deletions.\n");
         printf("              Will try to remove them as we go...\n");
      }
      

      /* Now set zones based on alignment                               */
      /* 01.12.04 Changed to tseqa rather than seqa[0]                  */
      /*
        SetNWZones(tseqa, seqb[0], MIN(strlen(tseqa), strlen(seqb[0])),
        NULL, NULL, strucnum);
      */

      SetNWZones(tseqa, mob_string, MIN(strlen(tseqa), 
                                        strlen(mob_string)),
                 NULL, NULL, strucnum); 
      
      /*free(seqb[0]);*/
      free(mob_string);
      
      if(++strucnum > gMultiCount)
      {
         printf("   Warning==> Alignment file contains more sequences than there\n");
         printf("              are structures.\n");
         break;
      }
   }

   if(strucnum < gMultiCount)
   {
      printf("   Warning==> Insufficient sequences in alignment file.\n");
      printf("              Fitting may fail!\n");
   }
   
   /* Free allocated memory and close file                              */
   free(tseqa);
   /*free(seqa[0]);*/
   free(ref_string);
   fclose(fp);
   
   /* Convert to sequential numbering with breaks between chains */
   ConvertAllZones(ZONE_MODE_RESNUM);
   ConvertAllZones(ZONE_MODE_SEQUENTIAL);
}


/************************************************************************/
/*>BOOL RemoveDoubleDeletions(char *seqa, char *seqb)
   --------------------------------------------------
   Remove deletions which appear in both sequences when reading an
   alignment file. This often occurs when the two sequences have come
   from part of a multiple alignment.

   13.12.95 Original    By: ACRM
   26.11.04 Fixed to allocate both arrays to longer length
*/
BOOL RemoveDoubleDeletions(char *seqa, char *seqb)
{
   char *copya = NULL,
        *copyb = NULL;
   int  i, j,
        lena,
        lenb,
        maxlen;

   lena = strlen(seqa);
   lenb = strlen(seqb);
   maxlen = MAX(lena, lenb);

   /* Create temporary storage for the sequences                        */
   copya = (char *)malloc((maxlen+1) * sizeof(char));
   copyb = (char *)malloc((maxlen+1) * sizeof(char));
   if(copya==NULL || copyb==NULL)
   {
      if(copya!=NULL) free(copya);
      if(copyb!=NULL) free(copyb);
      return(FALSE);
   }

   /* Copy in the sequences skipping any double deletions               */
   for(i=0, j=0; i<MAX(lena, lenb); i++)
   {
      if((seqa[i] != '-') || (seqb[i] != '-'))
      {
         if(i<lena)
            copya[j] = seqa[i];
         if(i<lenb)
            copyb[j] = seqb[i];
         j++;
      }
   }
   copya[j] = copyb[j] = '\0';
   
   /* Copy back into the original strings                               */
   strcpy(seqa, copya);
   strcpy(seqb, copyb);
   
   /* Free up the temporary storage                                     */
   free(copya);
   free(copyb);
   
   return(TRUE);
}

/************************************************************************/
/*>SetNWZones(char *ref_align, char *mob_align, int align_len,
              PDB **RefIndex, PDB **MobIndex, int strucnum)
   -----------------------------------------------------------
   Searches through the N&W sequence alignment and creates fitting zones
   from the equivalent regions.

   09.10.92 Original
   24.11.94 Fixed bug causing it to lose first zone in multi-zone match
   17.07.95 Replaced screen() with printf()
   18.07.95 Added initialisation of inserts in zones
   21.08.95 Fixed bug in additional non-existant zone being added when
            last zone not at end of chain
   22.11.95 Added check on deletion in both sequences
   13.12.95 Wasn't doing this check when stepping through a block of 
            deleteions. Fixed.
   18.06.96 Changed MODE_* to ZONE_MODE_*

   06.05.98 Completely rewritten! New version is 27% shorter, MUCH
            simpler and fixes a bug which occurred when a zone had only
            one residue.
   15.01.01 Simplified even further by making each residue an individual
            zone. this is not as elegant, but makes the implementation
            of distance checking much easier. If RefIndex and MobIndex
            are NULL, it behaves as before. If not then the distance
            between an atom pair is checked before adding the residue
            pair to the zone. Finally calls MergeZones() to merge
            adjacent zones.
   01.02.01 Added strucnum parameter
   20.02.01 Added check on gLimit[]
*/
void SetNWZones(char *ref_align, char *mob_align, int align_len,
                PDB **RefIndex, PDB **MobIndex, int strucnum)
{
   int   i,
         start,
         stop,
         ref_resnum  = 0,
         mob_resnum  = 0;
   ZONE  *z;

   if(gZoneList[strucnum])
   {
      FREELIST(gZoneList[strucnum], ZONE);
      gZoneList[strucnum] = NULL;
   }
   
   /* Get the start and stop of the region we are going to look at from
      gLimit[] if it has been specified. Otherwise just use the whole 
      alignment length
   */
   start = ((gLimit[0] < 1)||(gLimit[1] < 1)) ? 
            0 : (gLimit[0] - 1);
   stop  = ((gLimit[0] < 1)||(gLimit[1] < 1)) ? 
            align_len : (gLimit[1]);

   if(start > align_len)
      start = align_len-1;
   if(stop > align_len)
      stop  = align_len;

   for(i=0; i<start; i++)
   {
      /* Find offsets for first zone                                    */
      if(ref_align[i] != '-') ref_resnum++;
      if(mob_align[i] != '-') mob_resnum++;
   }

   for(i=start; i<stop; i++)
   {
      /* Find the residue number in each structure                      */
      if(ref_align[i] != '-') ref_resnum++;
      if(mob_align[i] != '-') mob_resnum++;
      
      if((ref_align[i] != '-') && (mob_align[i] != '-'))
      {
         if(((RefIndex==NULL) && (MobIndex==NULL)) ||
            (DISTSQ(RefIndex[ref_resnum-1], MobIndex[mob_resnum-1]) <=
             gMaxEquivDistSq))
         {
            /* Allocate and store the zone                              */
            if(gZoneList[strucnum])
            {
               /* Move to end of zone list                              */
               z=gZoneList[strucnum];
               LAST(z);
               ALLOCNEXT(z,ZONE);
            }
            else
            {
               INIT(gZoneList[strucnum],ZONE);
               z = gZoneList[strucnum];
            }
            if(z==NULL)
            {
               printf("   Error==> No memory for N&W fitting zones!\n");
               return;
            }
            
            z->chain1       = ' ';
            z->start1       = ref_resnum;
            z->startinsert1 = ' ';
            z->stop1        = ref_resnum;
            z->stopinsert1  = ' ';
            z->chain2       = ' ';
            z->start2       = mob_resnum;
            z->startinsert2 = ' ';
            z->stop2        = mob_resnum;
            z->stopinsert2  = ' ';
            z->mode         = ZONE_MODE_SEQUENTIAL;
         }
      }
   }

   MergeZones(strucnum);
   
   /* Set fitting flags                                                 */
   gFitted      = FALSE;
   gUserFitZone = TRUE;
}


/************************************************************************/
/*>void MergeZones(int strucnum)
   -----------------------------
   Merges zones describing sequentially numbered adjacent amino acids 

   15.01.01 Original   By: ACRM
   01.02.01 Added strucnum parameter
*/
void MergeZones(int strucnum)
{
   ZONE *z  = NULL,
        *zn = NULL;
   BOOL converged = TRUE;
   
   if(gZoneList[strucnum])
   {
      do
      {
         /* Assume we have converged                                    */
         converged = TRUE;
         for(z=gZoneList[strucnum]; z!=NULL; NEXT(z))
         {
            zn = z->next;
            if(zn)
            {
               /* If both zones are in sequential mode                  */
               if((z->mode  == ZONE_MODE_SEQUENTIAL) &&
                  (zn->mode == ZONE_MODE_SEQUENTIAL))
               {
                  /* See if the two zones are sequential                */
                  if((zn->start1 == (z->stop1 + 1)) &&
                     (zn->start2 == (z->stop2 + 1)))
                  {
                     z->stop1 = zn->stop1;
                     z->stop2 = zn->stop2;
                     z->next = zn->next;
                     free(zn);
                     converged = FALSE;
                  }
               }
            }
         }
      }  while(!converged);
   }
}



/************************************************************************/
/*>BOOL VerifySequence(char *seqa, char *seqb)
   -------------------------------------------
   Compare sequence A to sequence B ignoring deletions in sequence.

   23.04.08 Original By: CTP 
*/
BOOL VerifySequence(char *seqa, char *seqb)
{
   int i=0;
   int j=0;
   
   /* Return if no sequences */
   if(!seqa || !seqb)
      return(FALSE);
   
   /* Verify */
   for(i=0;i<strlen(seqa);i++)
   {
      /* Skip '-' in Sequence A */
      if(seqa[i] == '-')
         continue;
      
      /* Skip '-' in Sequence B */
      while(j<strlen(seqb) && seqb[j] == '-')
         j++;
      
      /* Return if hit end of Sequence B */
      if(j==strlen(seqb))
         return(FALSE);      
      
      /* DEBUG */
      /*
        printf("%c -- %c \n",seqa[i], seqb[j]);
      */
      
      /* Compare Residues */
      if(seqa[i] != seqb[j])
         return(FALSE);
      
      j++;
   }
   
   /* All tests completeted */
   return(TRUE);
}


/************************************************************************/
/*>BOOL BuildAlignment(char **seqs, ZONE *zones, char **alns, 
                       char **flagseq)
   ----------------------------------------------------------
   Inputs:    char **seqs      Sequences to align
              ZONE *zones      Linked list of zones
   Outputs:   char **alns      Aligned sequences
              char **flagseq   Seq of flags to indicate aligned resides
   Returns:   BOOL             OK?

   21.08.06 Original   By: ACRM
   17.06.08 Incorporated into ProFit with minor modification to ZONE 
            datatype By: CTP
*/
BOOL BuildAlignment(char **seqs, ZONE *zones, char **alns, char **flagseq)
{
   int i;
   int len[2];
   int idx[2], diff;
   AA  *aa[2];
   ZONE *z;
   
   /* Build the sequences as linked lists                               */
   if((aa[0] = BuildAAList(seqs[0]))==NULL)
      return(FALSE);
   if((aa[1] = BuildAAList(seqs[1]))==NULL)
      return(FALSE);
   
   /* Cycle through the zones                                           */
   for(z=zones; z!=NULL; NEXT(z))
   {
      /* Find the position in the linked list of the start or the 2
         equivalent zones
      */
      if((idx[0] = FindAAListOffsetByResnum(aa[0], z->start1))==(-1))
         return(FALSE);
      if((idx[1] = FindAAListOffsetByResnum(aa[1], z->start2))==(-1))
         return(FALSE);

      /* Determine the difference in the positions and make the appopriate
         number of inserts in the right sequence
      */
      if(idx[0] > idx[1])
      {
         diff = idx[0] - idx[1];
         if((aa[1] = 
             InsertResiduesInAAListAt(aa[1], '-', diff, idx[1]-1))==NULL)
            return(FALSE);
      }
      else
      {
         diff = idx[1] - idx[0];
         if((aa[0] = 
             InsertResiduesInAAListAt(aa[0], '-', diff, idx[0]-1))==NULL)
            return(FALSE);
      }

      /* Flag the aligned residues                                      */
      for(i=z->start1; i<=z->stop1; i++)
      {
         SetAAListFlagByResnum(aa[0], i);
      }
      for(i=z->start2; i<=z->stop2; i++)
      {
         SetAAListFlagByResnum(aa[1], i);
      }
   }

   /* Add inserts at the end of the shorter sequence                    */
   len[0] = GetAAListLen(aa[0]);
   len[1] = GetAAListLen(aa[1]);
   if(len[0] > len[1])
   {
      diff = len[0] - len[1];
      if((aa[1] = 
          InsertResiduesInAAListAt(aa[1], '-', diff, len[1]))==NULL)
         return(FALSE);
   }
   else
   {
      diff = len[1] - len[0];
      if((aa[0] = 
          InsertResiduesInAAListAt(aa[0], '-', diff, len[0]))==NULL)
         return(FALSE);
   }

   /* Convert the linked lists back into sequences                      */
   if((alns[0] = BuildSeqFromAAList(aa[0]))==NULL)
      return(FALSE);
   if((alns[1] = BuildSeqFromAAList(aa[1]))==NULL)
      return(FALSE);

   /* Build a sequence version of the flags to show aligned residues    */
   if((*flagseq = BuildFlagSeqFromAAList(aa[0], '*'))==NULL)
     return(FALSE);

   /* Free the linked lists                                             */
   FREELIST(aa[0], AA);
   FREELIST(aa[1], AA);

   return(TRUE);
}


/************************************************************************/
/*>BOOL BuildAlignment(char **seqs, ZONE *zones, char **alns, 
                       char **flagseq)
   ----------------------------------------------------------
   Inputs:    char **seqs      Sequences to align
              ZONE *zones      Linked list of zones
   Outputs:   char **alns      Aligned sequences
              char **flagseq   Seq of flags to indicate aligned resides
   Returns:   BOOL             OK?

   21.08.06 Original   By: ACRM
   17.06.08 Incorporated into ProFit with minor modification to ZONE 
            datatype By: CTP
   23.06.08 Modified to build alignment from residue A to residue B.
   25.11.08 Changed output format - unaligned residues are matched with 
            gaps.
*/
BOOL BuildAlignmentAB(char **seqs_in, ZONE *zones, char **alns, 
                      char **flagseq, int *res)
{
   int i;
   int len[2];
   int idx[2], diff;
   int insert_idx = 0;
   int insert_len[2];
   int zonelength = 0;
   AA  *aa[2];
   ZONE *z;
 
   char *seqs[2];

   /* Make truncated sequences By: CTP                                  */
   seqs[0] = seqs[1] = NULL;
   seqs[0] = malloc((strlen(seqs_in[0])+1) * sizeof(char));
   seqs[1] = malloc((strlen(seqs_in[1])+1) * sizeof(char));

   TruncateSeq(seqs[0], seqs_in[0], res[0], res[1]);
   TruncateSeq(seqs[1], seqs_in[1], res[2], res[3]);
  
   /* Build the sequences as linked lists                               */
   if((aa[0] = BuildAAList(seqs[0]))==NULL)
      return(FALSE);
   if((aa[1] = BuildAAList(seqs[1]))==NULL)
      return(FALSE);

   /* Run through the zones                                             */
   for(z=zones; z!=NULL; NEXT(z))
   {
      /* Check in region  By: CTP                                       */
      if((z->start1 < res[0]) || (z->stop1 > res[1]) ||
         (z->start2 < res[2]) || (z->stop2 > res[3]))
      {
         continue;
      }
      
      /* Find the position in the linked list of the start or the 2
         equivalent zones
      */
      if((idx[0] = 
          FindAAListOffsetByResnum(aa[0], z->start1 - res[0] + 1))==(-1))
         return(FALSE);
      if((idx[1] = 
          FindAAListOffsetByResnum(aa[1], z->start2 - res[2] + 1))==(-1))
         return(FALSE);
      
      /* Determine the difference in the positions and make the appopriate
         number of inserts in the right sequences
      */
      
      /* Check Length of Zone                                           */
      if(z->stop1 - z->start1 != z->stop2 - z->start2)
         return(FALSE);
      zonelength = z->stop1 - z->start1;
      
      /* Insert into reference seq                                      */
      insert_len[0] = idx[1] - insert_idx - 1;
      if((aa[0] = InsertResiduesInAAListAt(aa[0], '-', insert_len[0], 
                                           idx[0]-1))==NULL)
         return(FALSE);
      
      /* Insert into mobile seq                                         */
      insert_len[1] = idx[0] - insert_idx - 1;
      if((aa[1] = InsertResiduesInAAListAt(aa[1], '-', insert_len[1], 
                                           insert_idx))==NULL)
         return(FALSE);
      
      /* Update insert start point                                      */
      insert_idx += insert_len[0] + insert_len[1] + zonelength + 1;
      
      /* Flag the aligned residues                                      */
      for(i=z->start1 - res[0] + 1; i<=z->stop1 - res[0] + 1; i++)
      {
         SetAAListFlagByResnum(aa[0], i);
      }
      
      for(i=z->start2 - res[2] + 1; i<=z->stop2 - res[2] + 1; i++)
      {
         SetAAListFlagByResnum(aa[1], i);
      }
   }
   
   /* Deal with end of chains                                           */
   len[0] = GetAAListLen(aa[0]);
   len[1] = GetAAListLen(aa[1]);
   
   /* Find length to end of ref and mobile chains                       */
   insert_len[0] = len[1] - insert_idx;
   insert_len[1] = len[0] - insert_idx;
   
   if((aa[0] = InsertResiduesInAAListAt(aa[0], '-', insert_len[0], 
                                        len[0]))==NULL)
      return(FALSE);
   if((aa[1] = InsertResiduesInAAListAt(aa[1], '-', insert_len[1], 
                                        insert_idx))==NULL)
      return(FALSE);
   
   /* Add inserts at the end of the shorter sequence                    */
   len[0] = GetAAListLen(aa[0]);
   len[1] = GetAAListLen(aa[1]);
   if(len[0] > len[1])
   {
      diff = len[0] - len[1];
      if((aa[1] = 
          InsertResiduesInAAListAt(aa[1], '-', diff, len[1]))==NULL)
         return(FALSE);
   }
   else
   {
      diff = len[1] - len[0];
      if((aa[0] = 
          InsertResiduesInAAListAt(aa[0], '-', diff, len[0]))==NULL)
         return(FALSE);
   }
   
   /* Convert the linked lists back into sequences                      */
   if((alns[0] = BuildSeqFromAAList(aa[0]))==NULL)
      return(FALSE);
   if((alns[1] = BuildSeqFromAAList(aa[1]))==NULL)
      return(FALSE);
   
   /* Build a sequence version of the flags to show aligned residues    */
   if((*flagseq = BuildFlagSeqFromAAList(aa[0], '*'))==NULL)
      return(FALSE);
   
   /* Free the linked lists                                             */
   FREELIST(aa[0], AA);
   FREELIST(aa[1], AA);
   free(seqs[0]);
   free(seqs[1]);
   
   return(TRUE);
}


/************************************************************************/
/*>void  PrintSequence(char *sequence)
   -----------------------------------
   Prints sequence.

   21.08.06 Original   By: ACRM
   14.01.08 added output to file. By: CTP
*/
void PrintSequence(FILE *fp, char *sequence)
{ 
  int width  = 60;
  int i      =  0;
  int j      =  0;
  char line[MAXBUFF];

  if(strlen(sequence) > 0)
  {
     for(i=0, j=0; i<strlen(sequence); i++)
     {
        line[j++] = sequence[i];
        
        if(j == width || i == strlen(sequence)-1)
        {
           line[j] = '\0';
           fprintf(fp,"   %s\n",line);
           j=0;
        }
     }
  }
  else 
  {
     fprintf(fp,"   Undefined\n");
  }

  return;
}


/************************************************************************/
/*>void PrintNiceAlignment(char *ref_align, char *mob_align)
   ---------------------------------------------------------
   Prints a pairwise alignment in user-friendly format.

   23.07.08 Original By: CTP
   14.01.09 Added output to file.
*/
void PrintNiceAlignment(FILE *fp,char *ref_align, char *mob_align)
{
   int i,j;
   char refline[61];
   char mobline[61];
   int  width = 60;
   
   int alignlength = strlen(ref_align);
   
   /* Deal with zero length sequences.                                  */
   if(!alignlength)
   {
      fprintf(fp,"   Undefined\n\n");
      return;
   }
   
   /* Print sequence                                                    */
   for(i=0, j=0; i<alignlength; i++)
   {
      refline[j] = ref_align[i];
      mobline[j] = mob_align[i];
      j++;
      
      if(j == width || i == alignlength-1)
      {
         refline[j] = '\0';
         mobline[j] = '\0';
         
         fprintf(fp,"   %s\n",  refline);
         fprintf(fp,"   %s\n\n",mobline);
         j=0;
      }
   }
   fprintf(fp,"\n");
   return;
}


/************************************************************************/
/*>int KillChainBreak(char *outstring, char *instring)
   ---------------------------------------------------
   Removes chainbreak characters '*' from sequence instring.
   
   23.07.08 Original. By: CTP
*/
int KillChainBreak(char *outstring, char *instring)
{
   int x, y;
   
   for(x=0, y=0; x<strlen(instring); x++)
   {
      if(instring[x] != '*')
      {
         outstring[y++] = instring[x];
      }
   }
   outstring[y] = '\0';
   
   return(0);
}


/************************************************************************/
/*>int ChainBreakToGap(char *outstring, char *instring)
   ----------------------------------------------------
   Converts chainbreak characters '*' to gaps '-' in sequence instring.
   
   02.12.08 Original By: CTP
*/
int ChainBreakToGap(char *instring)
{
   int x;
   
   for(x=0; x<strlen(instring); x++)
   {
      if(instring[x] == '*')
         instring[x] = '-';
   }
   return(0);
}


/************************************************************************/
/*>int TruncateSeq(char *outstring, char *instring, int start, int stop)
   ---------------------------------------------------------------------
   Writes section of sequence between positions start and stop in 
   instring to outstring.
   
   23.07.08 Original. By: CTP
*/
int TruncateSeq(char *outstring, char *instring, int start, int stop)
{
   int x, y;

   for(x=start-1, y=0; x < stop; x++)
   {
      outstring[y++] = instring[x];
   }
   outstring[y] = '\0';

   return(0);
}


/************************************************************************/
/*>int AlignmentFromZones(BOOL fasta)
   --------------------------------
   Generates an alignment based on ProFit fitting Zones.
   
   The default output is a (user-friendly) pairwise alignment with the 
   reference and mobile sequences printed as pairs of 60-character wide 
   lines.

   The pir flag sets the printout to (machine-friendly) FASTA formatting 
   for the chain names and sequences.

   23.07.08 Original based on program by ACRM  By: CTP
   10.09.08 Ensured sequentially numbered fitting zones had breaks
            between chains.
   14.01.09 Added output to file.
*/
int AlignmentFromZones(char *filename, BOOL fasta)
{
   char *seqs[2], *alns[2], ids[2][MAXBUFF];
   char *flagseq = NULL;
   ZONE *chainlist[2];
   ZONE *ref, *mob, *z;
   int  i=0;
   FILE *fp = stdout;
   
   
   /* Set Pointers to NULL                                              */
   chainlist[0] = chainlist[1] = NULL;
   alns[0] = alns[1] = NULL;
   seqs[0] = seqs[1] = NULL;
   
   /* Convert to sequential numbering with breaks between chains        */
   if(ConvertAllZones(ZONE_MODE_RESNUM) ||
      ConvertAllZones(ZONE_MODE_SEQUENTIAL))
   {
      printf("   Error: Could not find zones.\n");
      return(1);
   }
   
   /* Sort Zones                                                        */
   SortAllZones();
   
   /* Set chainlist and sequence for reference                          */
   strcpy(ids[0],gRefFilename);
   chainlist[0] = ChainList(gRefPDB);
   seqs[0] = malloc((strlen(gRefSeq)+1) * sizeof(char));
   KillChainBreak(seqs[0],gRefSeq);
   
   /* Open output file/pipe                                             */
   if(filename)
   {
      if((fp=OpenOrPipe(filename))==NULL)
      {
         printf("   Warning: Unable to open output file\n");
         fp = stdout;
      }
   }
   
   for(i=0; i<MAXSTRUC; i++)
   {
      if(!gZoneList[i]) break;
      
      fprintf(fp,"\n   Mobile Structure: %2d\n",i+1);
      
      if(!OneToOneChains(i))
      {
         fprintf(fp,"\n   Warning: Chain aligned to more than one \
chain\n\n");
      }
      
      if(!SequentialZones(i))
      {
         fprintf(fp,"\n");
         fprintf(fp,"   Error: Could not convert zones to alignment.\n");
         fprintf(fp,"          Zones must not overlap and must be in\n");
         fprintf(fp,"          sequence along chain.\n\n");
         continue;
      }
      
      /* Set chainlist and sequence for mobile                          */
      chainlist[1] = ChainList(gMobPDB[i]);
      strcpy(ids[1],gMobFilename[i]);
      seqs[1] = malloc((strlen(gMobSeq[i])+1) * sizeof(char));
      KillChainBreak(seqs[1],gMobSeq[i]);
      
      for(ref = chainlist[0]; ref!=NULL; NEXT(ref))
      {
         for(mob = chainlist[1]; mob!=NULL; NEXT(mob))
         {
            for(z=gZoneList[i]; z!=NULL; NEXT(z))             
            {
               if((z->start1>=ref->start1)&&(z->stop1<=ref->stop1)&&
                  (z->start2>=mob->start1)&&(z->stop2<=mob->stop1)&&
                  (z->mode == ZONE_MODE_SEQUENTIAL))
               {
                  int position[4];
                  char chainids[2][MAXBUFF];
                  
                  position[0] = ref->start1;
                  position[1] = ref->stop1;
                  position[2] = mob->start1;
                  position[3] = mob->stop1;
                  
                  sprintf(chainids[0],"%s Chain '%c'",ids[0],ref->chain1);
                  sprintf(chainids[1],"%s Chain '%c'",ids[1],mob->chain1);
                  
                  /* NULL Alignment Pointers                            */
                  alns[0] = alns[1] = NULL;
                  flagseq = NULL;
                  
                  if(BuildAlignmentAB(seqs, gZoneList[i], alns, &flagseq, 
                                      position))
                  {
                     if(fasta)
                     {
                        /* Print FASTA Alignment                        */
                        fprintf(fp,"   >%s\n",chainids[0]);
                        PrintSequence(fp, alns[0]);
                        fprintf(fp,"\n");
                        fprintf(fp,"   >%s\n",chainids[1]);
                        PrintSequence(fp, alns[1]);
                        fprintf(fp,"\n\n");
                     }
                     else 
                     {
                        /* Print Nice Alignment                         */
                        fprintf(fp,"   %s\n   %s\n",
                                chainids[0], chainids[1]);
                        PrintNiceAlignment(fp,alns[0],alns[1]);
                     }
                  } 
                  
                  /* Free Alignment Memory                              */
                  if(alns[0]) free(alns[0]);
                  if(alns[1]) free(alns[1]);
                  if(flagseq) free(flagseq);
                  
                  /* NULL Alignment Pointers                            */
                  alns[0] = alns[1] = NULL;
                  flagseq = NULL;
                  
                  break;
               }
            }
         }
      }
      
      /* Free mobile sequence memory                                    */
      if(seqs[1]) free(seqs[1]);
      if(chainlist[1]) FREELIST(chainlist[1],ZONE);
      seqs[1]      = NULL;
      chainlist[1] = NULL;
   }
   
   /* Close output file/pipe                                            */
   if(fp != stdout)
      CloseOrPipe(fp);
   
   if(chainlist[0]) FREELIST(chainlist[0],ZONE);
   if(chainlist[1]) FREELIST(chainlist[1],ZONE);
   if(seqs[0]) free(seqs[0]);
   if(seqs[1]) free(seqs[1]);
   if(alns[0]) free(alns[0]);
   if(alns[1]) free(alns[1]);
   if(flagseq) free(flagseq);
   
   return(0);
}


/************************************************************************/
/*>ZONE *AlignToZone(char *ref_align, char *mob_align, int   ref_start,  
                     int   mob_start)
   --------------------------------------------------------------------
   Derive a linked list of zones from an alignment.

   23.07.08 Original based on SetNWZones() By: CTP
*/
ZONE *AlignToZone(char *ref_align, char *mob_align, 
                  int   ref_start,  int   mob_start)
{
  ZONE *zonelist = NULL;
  ZONE *z, *zn;
  int i, ref, mob;
  BOOL converged;

  for(i=0, ref=0, mob=0; i<strlen(ref_align); i++)
  {
     /* Add single-residue zone                                         */
     if((ref_align[i] != '-') && (mob_align[i] != '-'))
     {
        /* Allocate memory                                              */
        if(!zonelist)
        {
           INIT(zonelist,ZONE);
           z = zonelist;
        }
        else 
        {
           z = zonelist;
           LAST(z);
           ALLOCNEXT(z,ZONE);
        }
        
        /* Set Start and Stop                                           */
        z->chain1       = ' ';
        z->start1       = ref + ref_start;
        z->startinsert1 = ' ';
        z->stop1        = ref + ref_start;
        z->stopinsert1  = ' ';
        z->chain2       = ' ';
        z->start2       = mob + mob_start;
        z->startinsert2 = ' ';
        z->stop2        = mob + mob_start;
        z->stopinsert2  = ' ';
        z->mode         = ZONE_MODE_SEQUENTIAL;
        z->next         = NULL;
     }
     
     /* Increment residue count                                         */
     if(ref_align[i] != '-') ref++;
     if(mob_align[i] != '-') mob++;
  }
  
  
  /* Merge Zones                                                        */
  if(zonelist)
  {
     do
     {
        /* Assume we have converged                                     */
        converged = TRUE;
        for(z=zonelist; z!=NULL; NEXT(z))
        {
           zn = z->next;
           if(zn)
           {
              /* See if the two zones are sequential                    */
              if((zn->start1 == (z->stop1 + 1)) &&
                 (zn->start2 == (z->stop2 + 1)))
              {
                 z->stop1 = zn->stop1;
                 z->stop2 = zn->stop2;
                 z->next = zn->next;
                 free(zn);
                 converged = FALSE;
              }
           }
        }
     }  while(!converged);
  }
  
  return(zonelist);
}


/************************************************************************/
/*>void AlignChainStandard(int strucnum)
   -------------------------------------
   Default method of doing pairwise alignments in ProFit.
   Performs chain by chain pairwise alignments.

   23.07.08 Original based on NWAlign()  By: CTP
   29.08.08 Added normalised score. 
*/
void AlignChainStandard(int strucnum)
{
   int  ref_len, mob_len, align_len, score;
   char *ref_seq_all = NULL,
        *mob_seq_all = NULL,
        *ref_seq     = NULL,
        *mob_seq     = NULL,
        *ref_align   = NULL,
        *mob_align   = NULL;
   
   ZONE *ref_chain   = NULL,
        *mob_chain   = NULL,
        *ref, *mob;
   ZONE *zonelist, *z;
  
   printf("\n   Mobile Structure: %2d\n",strucnum+1);
   
   /* Get Sequences and Chains                                          */
   ref_seq_all = malloc((strlen(gRefSeq)+1) * sizeof(char));
   ref_seq     = malloc((strlen(gRefSeq)+1) * sizeof(char));
   KillChainBreak(ref_seq_all,gRefSeq);
   ref_chain = ChainList(gRefPDB);
   
   mob_seq_all = malloc((strlen(gMobSeq[strucnum])+1) * sizeof(char));
   mob_seq     = malloc((strlen(gMobSeq[strucnum])+1) * sizeof(char));
   KillChainBreak(mob_seq_all,gMobSeq[strucnum]);
   mob_chain = ChainList(gMobPDB[strucnum]);
   
   /* Check same number of chains                                       */
   /* Print warning if number of chains doesn't match.                  */
   if(!gQuiet)
   {
      for(ref = ref_chain, mob = mob_chain; ref != NULL && mob != NULL; 
          NEXT(ref), NEXT(mob))
      {}
      
      if(ref || mob)
      {
         printf("   Warning: Number of chains does not match.\n");
      }
   }
   
   /* Free gZoneList                                                    */
   if(gZoneList[strucnum])
   {
      FREELIST(gZoneList[strucnum],ZONE);
      gZoneList[strucnum] = NULL;
      gUserFitZone        = FALSE;
   }
   
   /* Align Chain By Chain                                              */
   for(ref = ref_chain, mob = mob_chain; ref != NULL && mob != NULL; 
       NEXT(ref), NEXT(mob))
   {
      /* Set sequences                                                  */
      TruncateSeq(ref_seq,ref_seq_all,ref->start1,ref->stop1);
      TruncateSeq(mob_seq,mob_seq_all,mob->start1,mob->stop1);
      
      ref_len = strlen(ref_seq);
      mob_len = strlen(mob_seq);
      
      ref_align = (char *)malloc((ref_len+mob_len)*sizeof(char));
      mob_align = (char *)malloc((ref_len+mob_len)*sizeof(char));
      
      score = affinealign(ref_seq, ref_len, mob_seq, mob_len, FALSE, 
                          FALSE, gGapPen, gGapPenExt, ref_align, 
                          mob_align, &align_len);
      
      ref_align[align_len] = '\0';
      mob_align[align_len] = '\0';
      
      printf("   %s Chain '%c'\n",gRefFilename,ref->chain1);
      printf("   %s Chain '%c'\n",gMobFilename[strucnum],mob->chain1);
      
      /* printf("   Score: %d\n",score); */
      printf("   Score: %d Normalised score: %.2f\n",
             score, (REAL)score/(REAL)(MIN(ref_len,mob_len)));
      PrintNiceAlignment(stdout,ref_align,mob_align);
      
      /* Convert to Zones                                               */
      zonelist = AlignToZone(ref_align, mob_align, ref->start1, 
                             mob->start1);
      
      /* Append zonelist to gZoneList                                   */
      if(!gZoneList[strucnum])
      {
         /* Make zone list new gZoneList                                */
         gZoneList[strucnum] = zonelist ;
         if(zonelist) gUserFitZone = TRUE;
      }
      else
      {
         /* Append zonelist to gZoneList                                */
         z = gZoneList[strucnum];
         LAST(z); 
         z->next = zonelist;
      }
      
      /* Free Memory                                                    */
      if(ref_align) free(ref_align);
      if(mob_align) free(mob_align);
   }
   
   /* Tidy                                                              */
   if(ref_seq_all) free(ref_seq_all);
   if(mob_seq_all) free(mob_seq_all);
   if(ref_seq)     free(ref_seq);
   if(mob_seq)     free(mob_seq);
   
   return;
}


/************************************************************************/
/*>void AlignWholeSequence(int strucnum)
   -------------------------------------
   Method of doing pairwise alignments in ProFit.
   Aligns whole sequence ignoring chain breaks. 

   23.07.08 Original based on NWAlign()  By: CTP
   29.08.08 Added normalised score. 
*/
void AlignWholeSequence(int strucnum)
{
   int  ref_len, mob_len, align_len, score;
   char *ref_seq_all = NULL,
        *mob_seq_all = NULL;
   char *ref_align   = NULL,
        *mob_align   = NULL;
   
   ZONE *zonelist; 

   /* Get Sequences and Chains                                          */
   ref_seq_all = malloc((strlen(gRefSeq)+1) * sizeof(char));
   KillChainBreak(ref_seq_all,gRefSeq);
   
   mob_seq_all = malloc((strlen(gMobSeq[strucnum])+1) * sizeof(char));
   KillChainBreak(mob_seq_all,gMobSeq[strucnum]);
   
   ref_len = strlen(ref_seq_all);
   mob_len = strlen(mob_seq_all);
   
   ref_align = (char *)malloc((ref_len+mob_len)*sizeof(char));
   mob_align = (char *)malloc((ref_len+mob_len)*sizeof(char));
   
   score = affinealign(ref_seq_all, ref_len, mob_seq_all, mob_len, 
                       FALSE, FALSE, gGapPen, gGapPenExt, 
                       ref_align, mob_align, &align_len);
   
   ref_align[align_len] = '\0';
   mob_align[align_len] = '\0';
   
   printf("   Mobile Structure: %2d\n",strucnum+1);
   printf("   %s vs %s\n",gRefFilename,gMobFilename[strucnum]);
   /* printf("   Score: %d\n",score); */
   printf("   Score: %d Normalised score: %.2f\n",
          score, (REAL)score/(REAL)(MIN(ref_len,mob_len)));
   PrintNiceAlignment(stdout,ref_align,mob_align);
   
   /* Convert to Zones                                                  */
   zonelist = AlignToZone(ref_align, mob_align, 1, 1);
   
   /* Set breaks between chains by converting to residue and back       */
   ConvertZoneList(zonelist, strucnum, ZONE_MODE_RESNUM);
   ConvertZoneList(zonelist, strucnum, ZONE_MODE_SEQUENTIAL);
   
   /* Replace global fit zones                                          */
   if(gZoneList[strucnum]) FREELIST(gZoneList[strucnum],ZONE);
   gZoneList[strucnum] = zonelist ;
   gUserFitZone = TRUE;
   
   /* Free Memory                                                       */
   if(ref_align) free(ref_align);
   if(mob_align) free(mob_align);
   
   /* Tidy                                                              */
   if(ref_seq_all) free(ref_seq_all);
   if(mob_seq_all) free(mob_seq_all);
   
   return;
}


/************************************************************************/
/*>void AlignZone(ZONE *alignzone, int strucnum, BOOL append)
   -------------------------------------
   Method of doing pairwise alignments in ProFit.
   Performs pairwise alignment on zone. 

   23.07.08 Original based on NWAlign()  By: CTP
   29.08.08 Added normalised score. 
*/
void AlignZone(ZONE *alignzone, int strucnum, BOOL append)
{
   int  ref_len, mob_len, align_len, score;
   char *ref_seq_all = NULL,
        *mob_seq_all = NULL,
        *ref_seq     = NULL,
        *mob_seq     = NULL,
        *ref_align   = NULL,
        *mob_align   = NULL;

   ZONE *ref_chain   = NULL,
        *mob_chain   = NULL;
   ZONE *zonelist, *z;

   char zone1[64], zone2[64];
   
   printf("\n   Mobile Structure: %2d\n",strucnum+1);

   /* Convert to sequential                                             */
   if(alignzone->mode == ZONE_MODE_RESNUM)
   {
      if(ConvertResidueToSequential(alignzone, strucnum))
      {
         printf("   Error: Error failed to find alignment zone.\n");
         return;
      }
   }

   /* Get Sequences and Chains                                          */
   ref_seq_all = malloc((strlen(gRefSeq)+1) * sizeof(char));
   ref_seq     = malloc((strlen(gRefSeq)+1) * sizeof(char));
   KillChainBreak(ref_seq_all,gRefSeq);
   ref_chain = ChainList(gRefPDB);

   mob_seq_all = malloc((strlen(gMobSeq[strucnum])+1) * sizeof(char));
   mob_seq     = malloc((strlen(gMobSeq[strucnum])+1) * sizeof(char));
   KillChainBreak(mob_seq_all,gMobSeq[strucnum]);
   mob_chain = ChainList(gMobPDB[strucnum]);

   /*  Select zone/zones to align                                       */
   if(!alignzone)
   {
      printf("   Error: Cannot align non-existant zone.\n");
      return;
   }

   if(alignzone->next)
   {
      printf("   Error: Stucture must be single zone.\n");
      return;
   }

   TruncateSeq(ref_seq,ref_seq_all,alignzone->start1,alignzone->stop1);
   TruncateSeq(mob_seq,mob_seq_all,alignzone->start2,alignzone->stop2);

   ref_len = strlen(ref_seq);
   mob_len = strlen(mob_seq);

   ref_align = (char *)malloc((ref_len+mob_len)*sizeof(char));
   mob_align = (char *)malloc((ref_len+mob_len)*sizeof(char));

   score = affinealign(ref_seq, ref_len, mob_seq, mob_len, FALSE, FALSE, 
                       gGapPen, gGapPenExt, ref_align, mob_align, 
                       &align_len);

   ref_align[align_len] = '\0';
   mob_align[align_len] = '\0';

   /* Format Zone                                                       */
   FormatZone(zone1,' ',alignzone->start1,' ',alignzone->stop1,' ');   
   FormatZone(zone2,' ',alignzone->start2,' ',alignzone->stop2,' ');

   /* Print Zone                                                        */
   printf("   %-16s vs %-16s (Sequential numbering)\n",zone1,zone2);
   /* printf("   Score: %d\n",score); */
   printf("   Score: %d Normalised score: %.2f\n",
          score, (REAL)score/(REAL)(MIN(ref_len,mob_len)));
   PrintNiceAlignment(stdout,ref_align, mob_align);

   /* Convert to alignment to zones                                     */
   zonelist = AlignToZone(ref_align, mob_align,
                          alignzone->start1, alignzone->start2);

   /* Add zones to zone list                                            */
   if(!gZoneList[strucnum])
   {
      /* Make zone list new gZoneList                                   */
      gZoneList[strucnum] = zonelist ;
      if(zonelist) gUserFitZone = TRUE;
   }
   else
   {
      if(append)
      {
         /* Append zonelist to gZoneList                                */
         z = gZoneList[strucnum];
         LAST(z); 
         z->next = zonelist;
      }
      else
      {
         /* Replace gZoneList with zonelist                             */
         FREELIST(gZoneList[strucnum], ZONE);
         gZoneList[strucnum] = zonelist;
         if(zonelist) gUserFitZone = TRUE;
      }
   }
   
   /* Case where no zones found ....                                    */
   if(!zonelist) printf("   Warning: No matching zones found...\n");

   /* Free Memory                                                       */
   if(ref_align) free(ref_align);
   if(mob_align) free(mob_align);

   /* Tidy                                                              */
   if(ref_seq_all) free(ref_seq_all);
   if(mob_seq_all) free(mob_seq_all);
   if(ref_seq)     free(ref_seq);
   if(mob_seq)     free(mob_seq);

   return;
}


/************************************************************************/
/*>int AlignmentWrapper(int strucnum, char *command, BOOL append)
   --------------------------------------------------------------
   Wrapper function calling AlignChainStandard(), AlignWholeSequence()
   or AlignZone(). 

   23.07.08 Original By: CTP
   30.10.08 Reset gFitted to FALSE
*/
int AlignmentWrapper(int strucnum, char *command, BOOL append)
{
   int struc, start, stop, method;
   ZONE *alignzone = NULL;

   /* Read Mutation Data Matrix File                                    */
   static int FirstCall = TRUE;
   if(FirstCall)
   {
      if(!ReadMDM(MDMFILE))
      {
         printf("   Error==> Unable to read mutation data matrix\n");
         return(1);
      }
      
      FirstCall = FALSE;
   }

   /* Set structures to align                                           */
   if(strucnum == -1)
   {
      /* Align all structures                                           */
      start = 0; 
      stop  = gMultiCount;
   }
   else 
   {
      /* Align single structure                                         */
      start = strucnum;
      stop  = strucnum + 1;
   }
   
   /* Parse Command                                                     */
   if(!command || strlen(command) == 0)
   {
      /* Do statandard align                                            */
      method = 0;
   }
   else if(!upstrncmp(command,"WHOLE",5) || !strcmp(command,"*"))
   {
      /* Perform whole sequence comparison                              */
      method = 1;
   }
   else 
   {
      /* Align zone                                                     */
      int  SeqZone = 0;
      int  start1, stop1, start2, stop2;
      char chain1, startinsert1, stopinsert1,
           chain2, startinsert2, stopinsert2;
      
      SeqZone = ParseZone(command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1, 
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2,
                          strucnum);
      
      if(SeqZone == (-2))
      {
         printf("   Error==> You cannot specify zones for both the \
reference\n");
         printf("            and mobile structures when performing \
multiple\n");
         printf("            structure fitting.\n");
         return(1);
      }
      
      INIT(alignzone, ZONE);
      alignzone->chain1       = chain1;       
      alignzone->start1       = start1;       
      alignzone->startinsert1 = startinsert1; 
      alignzone->stop1        = stop1;        
      alignzone->stopinsert1  = stopinsert1;  
      alignzone->chain2       = chain2;       
      alignzone->start2       = start2;       
      alignzone->startinsert2 = startinsert2; 
      alignzone->stop2        = stop2;        
      alignzone->stopinsert2  = stopinsert2;  
      alignzone->mode         = SeqZone ? ZONE_MODE_SEQUENTIAL : 
                                          gCurrentMode;     
      alignzone->next         = NULL;
      
      method = 2;
   }
   
   
   /* Perform Alignment                                                 */
   for(struc=start; struc<stop; struc++)
   {
      switch(method)
      {
      case 0: /* Standard Chain by Chain Alignment                      */
         AlignChainStandard(struc);
         break;
         
      case 1: /* Whole Sequence Alignment                               */
         AlignWholeSequence(struc);
         break;
         
      case 2: /* Align Zone                                             */
         AlignZone(alignzone, struc, append);
         break;
         
      default:
         break;
      }
   }
   
   /* Reset fitted flag                                                 */
   gFitted = FALSE;
   gUserFitZone = TRUE;
   
   /* Tidy Up                                                           */
   if(alignzone) FREELIST(alignzone, ZONE);
   return(0);
}


/************************************************************************/
/*>BOOL CommonZones(void)
   ----------------------
   Test to check if ZONEs for each mobile structure are identical.

   13.01.09 Original by CTP.
*/
BOOL CommonZones(void)
{
   int i;
   ZONE *z1, *z2;
   
   /* Convert to sequential numbering with breaks between chains        */
   /* Note - Could do this before calling function.                     */
   if(ConvertAllZones(ZONE_MODE_RESNUM) ||
      ConvertAllZones(ZONE_MODE_SEQUENTIAL))
   {
      printf("   Error: Could not find zones.\n");
      return(FALSE);
   }
   
   /* Sort Zones                                                        */
   SortAllZones();
   
   /* Loop - Structures                                                 */
   for(i=1;i<gMultiCount;i++)
   {
      /* Loop - Zones                                                   */
      z1 = gZoneList[0];
      z2 = gZoneList[i];
      for(;z1 != NULL && z2 != NULL; NEXT(z1), NEXT(z2))
      {
         if((gZoneList[0]->chain1       != gZoneList[i]->chain1      ) ||
            (gZoneList[0]->start1       != gZoneList[i]->start1      ) ||
            (gZoneList[0]->startinsert1 != gZoneList[i]->startinsert1) ||
            (gZoneList[0]->stop1        != gZoneList[i]->stop1       ) ||
            (gZoneList[0]->stopinsert1  != gZoneList[i]->stopinsert1 ) ||
            (gZoneList[0]->mode         != gZoneList[i]->mode        ))
         {
            return(FALSE);
         }
      }
      
      if(z1 || z2) return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>int AlignmentFromZones_PIR(void)
   --------------------------------
   Generates an alignment based on ProFit fitting Zones.
   
   Output is a multiple sequence alignment with the reference and mobile 
   sequences printed as series of 60-character wide lines.

   The format is machine-readable (well... ProFit-readable) but requires
   the zones generated by ProFit to be:

   a) Non-Overlapping
   b) Sequential along the WHOLE sequence - NOT just sequential for each 
      individual chain.
      (e.g. The user can't align Chain A with B then B with A.)

   02.12.08 Original based on AlignmentFromZones(). By: CTP
   13.01.09 Removed requirement for common zones.
   14.01.09 Added output to file.
   29.01.09 Fixed bug in checking sequetial zones. (Only checked for 
            each individual chain) and added additional error checking at
            start of function.
*/
int AlignmentFromZones_PIR(char *filename)
{
   char *seqs[MAXSTRUC + 1], *alns[MAXSTRUC + 1];
   int  i;
   FILE *fp = stdout;
   
   /*** Check Zones can be output in PIR Format                       ***/
   /* Check for User-Defined Zones                                      */
   if(!gUserFitZone)
   {
      printf("   Error: No user-defined zones found.\n");
      return(1);
   }
   
   /* Convert to sequential numbering with breaks between chains        */
   if(ConvertAllZones(ZONE_MODE_RESNUM) ||
      ConvertAllZones(ZONE_MODE_SEQUENTIAL))
   {
      printf("   Error: Could not convert zones.\n");
      return(1);
   }
   
   /* Check for sequential zones across whole sequence                  */
   for(i=0; i<gMultiCount; i++)
   {
      if(!SequentialZonesWholeSeq(i))
      {
         printf("\n");
         printf("   Error: Could not convert zones to PIR alignment.\n");
         printf("          Zones must not overlap and must be in\n");
         printf("          sequence over whole sequence.\n\n");
         return(1);
      }
   }
   
   /*** Make Sequences                                                ***/
   /* Set reference sequence as Sequence zero */
   seqs[0] = malloc(strlen(gRefSeq)+1 *sizeof(char));
   KillChainBreak(seqs[0],gRefSeq);
   
   /* Set Mobile Sequences                                              */
   for(i=1;i<=gMultiCount;i++)
   {
      seqs[i] = malloc(strlen(gMobSeq[i-1])+1 *sizeof(char));
      KillChainBreak(seqs[i],gMobSeq[i-1]);
   }
   
   /* Build Multiple Alignment from Zones                               */
   if(!BuildMultiAlignment(seqs, alns))
   {
      printf("   Error: Could not create alignment from user zones.\n");
      return(1);
   }
   
   /*** Print PIR Alignment                                           ***/
   /* Open output file/pipe                                             */
   if(filename)
   {
      if((fp=OpenOrPipe(filename))==NULL)
      {
         printf("   Warning: Unable to open output file\n");
         fp = stdout;
      }
   }
   
   /* Ref Sequence                                                      */
   fprintf(fp,">P1;REFSEQ\n");
   fprintf(fp,"Reference Sequence - %s\n",gRefFilename);
   PrintSequencePIR(fp,alns[0],60);
   
   /* Mobile Sequence                                                   */
   for(i=1;i<=gMultiCount;i++)
   {
      fprintf(fp,"\n");
      fprintf(fp,">P1;M_%04d\n",i);
      fprintf(fp,"Mobile Sequence - %s\n",gMobFilename[i-1]);
      PrintSequencePIR(fp,alns[i],60);
   }
   
   if(fp == stdout) printf("\n\n");
   
   /* DEBUG: Print Multiple Alignment */
/***
   for(i=0;i<=gMultiCount;i++) printf("%s\n",alns[i]);
***/
   
   /* Close output file/pipe                                            */
   if(fp != stdout)
      CloseOrPipe(fp);
   
   /* Free Memory                                                       */
   for(i=0;i<gMultiCount+1;i++)
   {
      if(seqs[i]) free(seqs[i]);
      if(alns[i]) free(alns[i]);
   }
   
   return(0);
}


/************************************************************************/
/*>void  PrintSequencePIR(char *sequence)
   --------------------------------------
   Prints PIR sequence.

   17.12.08 Original based on PrintSequence(). By: CTP
   14.01.08 Added output to file.
*/
void PrintSequencePIR(FILE *fp, char *sequence, int width)
{ 
   int i      =  0;
   int j      =  0;
   char line[MAXBUFF];
   
   if(strlen(sequence) > 0)
   {
      for(i=0, j=0; i<strlen(sequence); i++)
      {
         line[j++] = sequence[i];
         
         if(j == width || i == strlen(sequence)-1)
         {
            line[j] = '\0';
            fprintf(fp,"%s\n",line);
            j=0;
         }
      }
   }
   else 
   {
      fprintf(fp,"Undefined\n");
   }
   
   return;
}


/************************************************************************/
/*>BOOL BuildMultiAlignment(char **seqs, char **alns)
   --------------------------------------------------
   Inputs:    char **seqs      Sequences to align
   Outputs:   char **alns      Aligned sequences
   Returns:   BOOL             OK?

   Builds multiple sequence alignment based on zones for output in PIR 
   format. 

   Mobile sequences are added one at a time and gaps for unmatched 
   residues are added to the reference sequence and previously included
   sequences. Gaps are then added to the mobile sequence. The ends of the
   alignment are dealt with and chain breaks are added. 

   The zones used to create the alignment do not have to be identical for
   each mobile but must occur sequentially along the whole structure. In
   other words, zones can't be converted into a sequence alignment unless
   they're in sequence.
   eg 
   Zones 1-5:1-5, 6-10:11-15, 11-15:6-10 cannot be converted into a 
   sequence alignment.

   13.01.09 Original By: CTP
   29.01.09 Fixed bug handling chain ends.
*/
BOOL BuildMultiAlignment(char **seqs, char **alns)
{
   ZONE *mobzone[MAXSTRUC], *chainlist[MAXSTRUC +1]; /* Zone lists      */
   ZONE *z;                /* Zone pointer                              */
   AA   *aa[MAXSTRUC + 1]; /* AA List                                   */
   AA   *ref_aa, *mob_aa;  /* AA pointers                               */
   
   int mobile;             /* Current mobile.                           */
   int i, j, len, idx[2]; 
   int insert_idx = 0;
   int insert_len[2];
   int max_len = 0;
   int length_ref = 0;
   int length_mob = 0; 
   
   
   /* Convert sequences into AA lists                                   */
   for(i=0;i<=gMultiCount;i++)
   {
      if((aa[i] = BuildAAList(seqs[i]))==NULL) 
         return(FALSE);
   }
   
   /* Set zone pointer for each sequence                                */
   for(i=0;i<gMultiCount;i++)
   {
      mobzone[i] = gZoneList[i];
      if(mobzone[i] == NULL) 
         return(FALSE);
   }
   
   /* Insert Gaps for each mobile                                       */
   /* ---------------------------                                       */
   
   /* Cycle Through Mobile Sequences                                    */
   for(mobile=1; mobile<=gMultiCount; mobile++)
   {
      /* Reset Insert Index                                             */
      insert_idx = 0; 
      
      /* Cycle Through Mobile Zones                                     */
      for(z=gZoneList[mobile - 1]; z!=NULL; NEXT(z))
      {
         /* Find offset for reference                                   */
         if((idx[0] = FindAAListOffsetByResnum(aa[0], z->start1))==(-1))
            return(FALSE);
         
         /* Find Offset for current mobile                              */
         if((idx[1] = FindAAListOffsetByResnum(aa[mobile], 
                                               z->start2))==(-1))
            return(FALSE);
         
         /* Insert Length for Reference Sequence Gap                    */
         insert_len[0] = idx[1] - insert_idx - 1;
         
         /* Insert Length for Mobile Sequence Gap                       */
         insert_len[1] = idx[0] - insert_idx - 1;
         
         /* Add Gap to Reference and Previous Mobiles                   */
         for(i=0;i<mobile;i++)
         {
            if((aa[i] = InsertResiduesInAAListAt(aa[i], '-', 
                                                 insert_len[0], 
                                                 idx[0]-1))==NULL) 
               return(FALSE);
         }
         
         /* Add Gap to Current Mobile                                   */
         if((aa[mobile] = InsertResiduesInAAListAt(aa[mobile], '-', 
                                                   insert_len[1],
                                                   insert_idx))==NULL)
         {
            return(FALSE);      
         }
         
         /* Update Insert Index                                         */
         insert_idx = FindAAListOffsetByResnum(aa[0], z->stop1);
         
         /* Insert gaps within zone.                                    */
         /* ------------------------                                    */
         
         /* Find start of zone                                          */
         idx[0] = FindAAListOffsetByResnum(aa[0], z->start1); 
         
         /* Set pointers for Ref and Mob sequences                      */
         ref_aa = aa[0];
         mob_aa = aa[mobile];
         
         for(i=1; i<idx[0]; i++)
         {
            NEXT(ref_aa);
            NEXT(mob_aa);
         }
         
         /* Loop through zone                                           */
         for(i=idx[0]; i<=insert_idx; i++)
         {
            /* Splice gaps into aa list                                 */
            if(ref_aa->res == '-')
            {
               AA *tmp_aa = NULL;
               tmp_aa = InsertResidueInAAListAt(aa[mobile], '-', i-1);
               
               /* No memory                                             */
               if(!tmp_aa) return(FALSE);
               
               /* Do we need to reset pointer?                          */
               if(tmp_aa != aa[mobile])
               {
                  /* Reset pointer                                      */
                  aa[mobile] = tmp_aa;
                  mob_aa    = aa[mobile];
                  for(j=0; j<i && mob_aa != NULL; j++) NEXT(mob_aa);
               }
            }
            else 
            {
               NEXT(mob_aa);
            }
            
            NEXT(ref_aa);
         }
         
      }  /* End of zones loop                                           */
      
      
      /* Deal with unaligned ends                                       */
      /* ------------------------                                       */
      
      /* Find length of sequences                                       */
      length_ref = GetAAListLen(aa[0]); 
      length_mob = GetAAListLen(aa[mobile]); 
      
      /* Find Final Mobile Zone                                         */
      z=gZoneList[mobile - 1];
      LAST(z);

      /* Find offset for reference zones finish                         */
      if((idx[0] = FindAAListOffsetByResnum(aa[0], z->stop1))==(-1))
         return(FALSE);
      
      /* Find Offset for current mobile zones finish                    */
      if((idx[1] = FindAAListOffsetByResnum(aa[mobile], z->stop2))==(-1))
         return(FALSE);
      
      /* Insert Length for Reference Sequence Gap                       */
      insert_len[0] = length_mob - idx[1];
      
      /* Insert Length for Mobile Sequence Gap                          */
      insert_len[1] = length_ref - idx[0];
      
      /* Add Gap to End of Reference                                    */
      if((aa[0] = InsertResiduesInAAListAt(aa[0], '-', insert_len[0], 
                                           length_ref))==NULL)
         return(FALSE);
      
      /* Add Gap to Mobile                                              */
      if((aa[mobile] = InsertResiduesInAAListAt(aa[mobile], '-', 
                                                insert_len[1],
                                                idx[1]))==NULL)
         return(FALSE);
   }  /* End of mobiles loop                                            */
   
   
   /* Add inserts at the end of the shorter sequences                   */
   for(i=0;i<gMultiCount+1;i++)
   {
      len = GetAAListLen(aa[i]);
      max_len = (len > max_len) ? len : max_len;
   }
   
   for(i=0;i<gMultiCount+1;i++)
   {
      len = GetAAListLen(aa[i]);
      
      if(max_len > len)
      {
         if((aa[i] = InsertResiduesInAAListAt(aa[i], '-', max_len-len,
                                              len))==NULL)
            return(FALSE);
      }
   }
   
   /* Add Chain Breaks                                                  */
   /* ----------------                                                  */
   
   /* Get Chain breaks                                                  */
   chainlist[0] = ChainList(gRefPDB);
   for(i=1; i<gMultiCount+1;i++)
   {
      chainlist[i] = ChainList(gMobPDB[i-1]);
   }
   
   /* Cycle Through Structures                                          */
   for(i=0; i<gMultiCount+1;i++)
   {
      ZONE *chain = NULL;
      BOOL found  = FALSE;
      
      /* Cycle Through Chains                                           */
      chain = chainlist[i];
      NEXT(chain);
      
      for(;chain != NULL; NEXT(chain))
      {
         int index = 0;
         /* Find offset of chainbreak                                   */
         if((index = FindAAListOffsetByResnum(aa[i], 
                                              chain->start1))==(-1))
            return(FALSE);
         
         /* Is this one we did earlier?                                 */
         for(j=0; j<i-1 && !found; j++)
         {
            z=chainlist[j];
            NEXT(z);
            for(;z!=NULL && !found;NEXT(z))
            {
               int curr_index = 0;
               /* Find current offset                                   */
               if((curr_index = 
                   FindAAListOffsetByResnum(aa[j], z->start1))==(-1))
                  return(FALSE);
               
               if(curr_index == index) 
                  found = TRUE;
            } 
         }
         
         /* Insert new column into alignment and add breaks             */
         if(!found)
         {
            /* Add column to alignment                                  */
            for(j=0; j<gMultiCount+1;j++)
            {
               if(i == j)
               {
                  /* Add break                                          */
                  if((aa[j] = InsertResiduesInAAListAt(aa[j],'*',1,
                                                       index-1))==NULL) 
                     return(FALSE);
               }
               else
               {
                  /* Add gap                                            */
                  if((aa[j]=InsertResiduesInAAListAt(aa[j],'-',1,
                                                     index-1))==NULL) 
                     return(FALSE);
               }
            }             
         }
         else    
         {
            AA *aa_index = aa[i];
            
            /* Set pointer                                              */
            for(j=1;j<index-1;j++)
            {
               NEXT(aa_index);
            }
            
            /* Delete gap then insert break                             */
            DELETE(aa[i],aa_index,AA);
            
            if((aa[i]=InsertResiduesInAAListAt(aa[i],'*',1,index-2))
               == NULL) 
               return(FALSE);
         }
      }
   }
   
   /* Add Chainbreak to end of sequence                                 */
   for(i=0; i<gMultiCount+1;i++)
   {
      len = GetAAListLen(aa[i]);
      if((aa[i] = InsertResiduesInAAListAt(aa[i], '*',1, len))==NULL)
         return(FALSE);
   }
   
   /* Convert the linked lists back into sequences                      */
   for(i=0;i<gMultiCount+1;i++)
   {
      if((alns[i] = BuildSeqFromAAList(aa[i]))==NULL)
         return(FALSE);
   }
   
   /* Free Linked Lists                                                 */
   for(i=0;i<gMultiCount+1;i++)
   {
      FREELIST(aa[i],        AA);
      FREELIST(chainlist[i], ZONE);
   }
   
   return(TRUE);
}

/************************************************************************/
/**

   \file       RdSeqPDB.c
   
   \version    V1.5
   \date       14.07.22
   \brief      Read sequence from SEQRES records in a PDB file
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 1996-2022
   \author     Prof. Andrew C. R. Martin
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0   14.10.96 Original   By: ACRM
-  V1.1   25.03.14 Added CHAINMATCH. By: CTP
-  V1.2   07.07.14 Use bl prefix for functions By: CTP
-  V1.3   26.02.15 Added blReadSeqresWholePDB()  By: ACRM
-  V1.4   17.11.21 Added blFixSequence()
-  V1.5   14.07.22 GAPPEN in blFixSequence() now zero

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP File IO
   #FUNCTION  blReadSeqresPDB()
   Reads the sequence from the SEQRES records of a PDB file

   #FUNCTION  blReadSeqresWholePDB()
   Reads the sequence from the SEQRES records from header data stored in
   a WHOLEPDB structure

   #FUNCTION  blFixSequence()
   Combines SEQRES and ATOM sequence data
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "general.h"
#include "seq.h"
#include "macros.h"
#include "fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))
#define MAXBUFF   160
#define MAXCHAINS 240
#define GAPPEN    0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static STRINGLIST *RdSeqRes(FILE *fp);
static STRINGLIST *RdSeqResHeader(WHOLEPDB *wpdb);
static char *sCombineSequence(char *align1, char *align2, int align_len,
                              BOOL upper);
extern char *strdup(const char *s);


/************************************************************************/
/*>char **blReadSeqresPDB(FILE *fp, int *nchains)
   ----------------------------------------------
*//**

   \param[in]     *fp       PDB file pointer
   \param[out]    *nchains  Number of chains found
   \return                  Array of sequence strings

   Reads the sequence from the SEQRES records of a PDB file. Creates
   an array of malloc()'d character arrays in which the sequence is
   stored. Can therefore cope with any size of sequence information
   from the PDB file.

   This is not normally recommended to get the sequence for a PDB file
   this way, but is useful to detect discrepancies compared with the
   sequence described by the ATOM records.

-  14.10.96 Original   By: ACRM
-  25.03.14 Added CHAINMATCH. Chain IDs handled as strings. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
*/
char **blReadSeqresPDB(FILE *fp, int *nchains)
{
   STRINGLIST *seqres = NULL, 
              *s;
   char       currchain[2] = " ",
              chain[2]     = " ",
              **seqs,
              res[13][8];
   int        chainnum = 0,
              nres     = 0,
              i;

   *nchains = 0;
   
   /* First read the SEQRES records into a linked list                  */
   if((seqres = RdSeqRes(fp))==NULL)
      return(NULL);

   /* FIRST PASS: See how many chains there are                         */
   strncpy(currchain,&(seqres->string[11]),1);
   *nchains  = 1;
   for(s=seqres; s!=NULL; NEXT(s))
   {
      strncpy(chain,&(s->string[11]),1);
      if(!CHAINMATCH(chain,currchain))
      {
         strncpy(currchain,chain,1);
         (*nchains)++;
      }
   }

   /* Allocate an array of character pointers to store this number of
      strings
   */
   if((seqs=(char **)malloc((*nchains) * sizeof(char *)))==NULL)
   {
      FREELIST(seqres, STRINGLIST);
      return(NULL);
   }

   /* SECOND PASS: Allocate space to store each chain                   */
   chainnum  = 0;
   strcpy(currchain,"");
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s%5d",chain,&nres);
      if(!CHAINMATCH(chain,currchain))
      {
         strcpy(currchain,chain);
         if((seqs[chainnum]=(char *)malloc((nres+1)*sizeof(char))) 
            == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
         chainnum++;
      }
   }

   /* THIRD PASS: Store the sequence                                    */
   chainnum  = 0;
   nres      = 0;
   strncpy(currchain,&(seqres->string[11]),1);
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s%7x%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s",
              chain,res[0],res[1],res[2],res[3],res[4],res[5],res[6],
              res[7],res[8],res[9],res[10],res[11],res[12]);
      if(!CHAINMATCH(chain,currchain))
      {
         /* Start of new chain, terminate last one                      */
         seqs[chainnum][nres] = '\0';
         strcpy(currchain,chain);
         nres      = 0;
         chainnum++;
      }
      
      /* Store these sequence data                                      */
      for(i=0; i<13; i++)
      {
         /* Break out if not all positions were filled in               */
         if(res[i][0] == ' ')
            break;
         seqs[chainnum][nres++] = blThrone(res[i]);
      }
   }
   /* Terminate last chain                                              */
   seqs[chainnum][nres] = '\0';

   FREELIST(seqres, STRINGLIST);
   
   return(seqs);
}


/************************************************************************/
/*>char **blReadSeqresWholePDB(WHOLEPDB *wpdb, int *nchains)
   ---------------------------------------------------------
*//**

   \param[in]     wpdb      WHOLEPDB structure
   \param[out]    *nchains  Number of chains found
   \return                  Array of sequence strings

   Reads the sequence from the SEQRES records from the PDB header
   stored in a WHOLEPDB structure. Creates an array of malloc()'d
   character arrays in which the sequence is stored. Can therefore
   cope with any size of sequence information from the PDB file.

   This is not normally recommended to get the sequence for a PDB file
   this way, but is useful to detect discrepancies compared with the
   sequence described by the ATOM records.

-  26.02.15 Original based on blReadSeqresPDB()   By: ACRM
*/
char **blReadSeqresWholePDB(WHOLEPDB *wpdb, int *nchains)
{
   STRINGLIST *seqres = NULL, 
              *s;
   char       currchain[2] = " ",
              chain[2]     = " ",
              **seqs,
              res[13][8];
   int        chainnum = 0,
              nres     = 0,
              i;

   *nchains = 0;
   
   /* First read the SEQRES records into a linked list                  */
   if((seqres = RdSeqResHeader(wpdb))==NULL)
      return(NULL);

   /* FIRST PASS: See how many chains there are                         */
   strncpy(currchain,&(seqres->string[11]),1);
   *nchains  = 1;
   for(s=seqres; s!=NULL; NEXT(s))
   {
      strncpy(chain,&(s->string[11]),1);
      if(!CHAINMATCH(chain,currchain))
      {
         strncpy(currchain,chain,1);
         (*nchains)++;
      }
   }

   /* Allocate an array of character pointers to store this number of
      strings
   */
   if((seqs=(char **)malloc((*nchains) * sizeof(char *)))==NULL)
   {
      FREELIST(seqres, STRINGLIST);
      return(NULL);
   }

   /* SECOND PASS: Allocate space to store each chain                   */
   chainnum  = 0;
   strcpy(currchain,"");
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s%5d",chain,&nres);
      if(!CHAINMATCH(chain,currchain))
      {
         strcpy(currchain,chain);
         if((seqs[chainnum]=(char *)malloc((nres+1)*sizeof(char))) 
            == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
         chainnum++;
      }
   }

   /* THIRD PASS: Store the sequence                                    */
   chainnum  = 0;
   nres      = 0;
   strncpy(currchain,&(seqres->string[11]),1);
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s%7x%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s",
              chain,res[0],res[1],res[2],res[3],res[4],res[5],res[6],
              res[7],res[8],res[9],res[10],res[11],res[12]);
      if(!CHAINMATCH(chain,currchain))
      {
         /* Start of new chain, terminate last one                      */
         seqs[chainnum][nres] = '\0';
         strcpy(currchain,chain);
         nres      = 0;
         chainnum++;
      }
      
      /* Store these sequence data                                      */
      for(i=0; i<13; i++)
      {
         /* Break out if not all positions were filled in               */
         if(res[i][0] == ' ')
            break;
         seqs[chainnum][nres++] = blThrone(res[i]);
      }
   }
   /* Terminate last chain                                              */
   seqs[chainnum][nres] = '\0';

   FREELIST(seqres, STRINGLIST);
   
   return(seqs);
}



/************************************************************************/
/*>static STRINGLIST *RdSeqRes(FILE *fp)
   -------------------------------------
*//**

   \param[in]     *fp      PDB File pointer
   \return                 Linked list of SEQRES records

   Used by ReadSeqresPDB() to read the SEQRES records into a linked list.

-  14.10.96 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
static STRINGLIST *RdSeqRes(FILE *fp)
{
   static STRINGLIST *seqres = NULL;
   char              buffer[MAXBUFF];
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(!strncmp(buffer,"SEQRES",6))
      {
         if((seqres = blStoreString(seqres, buffer)) == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
      }
   }
   
   return(seqres);
}

/************************************************************************/
/*>static STRINGLIST *RdSeqResHeader(WHOLEPDB *wpdb)
   -------------------------------------------------
*//**

   \param[in]     *fp      PDB File pointer
   \return                 Linked list of SEQRES records

   Used by ReadSeqresPDB() to read the SEQRES records into a linked list.

-  26.02.15 Original based on RdSeqRes()   By: ACRM
*/
static STRINGLIST *RdSeqResHeader(WHOLEPDB *wpdb)
{
   STRINGLIST *seqres = NULL,
              *s;
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string,"SEQRES",6))
      {
         if((seqres = blStoreString(seqres, s->string)) == NULL)
         {
            FREELIST(seqres, STRINGLIST);
            return(NULL);
         }
      }
   }
   
   return(seqres);
}


/************************************************************************/
/*>char *blFixSequence(char *seqresSequence, char *atomSequence,
                       char **seqresChains, char **atomChains,
                       char **outChains, BOOL IgnoreSEQRES, 
                       int nAtomChains, BOOL upper, BOOL quiet, 
                       char *label)
   -----------------------------------------------------------------------
*//**
   \param[in]  *seqresSequence   Sequence extracted from SEQRES records 
                                 with a * separating each chain
   \param[in]  *atomSequence     Sequence extracted from ATOM records 
                                 with a * separating each chain
   \param[in]  **seqresChains    Chain labels for SEQRES
   \param[in]  **atomChains      Chain labels for ATOMS
   \param[out] **outChains       Output chain labels
   \param[in]  IgnoreSEQRES      Ignore SEQRES records for completely
                                 missing chains
   \param[in]  nAatomChains      Number of ATOM chains
   \param[in]  upper             Make output sequence all uppercase
   \param[in]  quiet             No warnings about missing residues
   \param[in]  label             Label to use in warnings (or NULL)
   \return                       Corrected sequence combining ATOM and
                                 SEQRES records (malloc'd)

   Create a final output sequence by combining the information from the
   ATOM and SEQRES records.

   Originally part of pdb2pir

-  21.08.97 Original   By: ACRM
-  22.08.97 Added chain information
-  26.08.97 A couple of bug fixes in initialising allocated memory
            Warning/error messages give the label if specified
-  22.05.09 Added IgnoreSEQRES to ignore chains that are in SEQRES but
            not in ATOM records
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  27.07.21 Returns the atomSequence if SEQRES not read
*/
char *blFixSequence(char *seqresSequence, char *atomSequence,
                    char **seqresChains, char **atomChains,
                    char **outChains, BOOL IgnoreSEQRES, int nAtomChains,
                    BOOL upper, BOOL quiet, char *label)
{
   int  i, j, len, len1, len2,
        nchain[2],
        align_len,
        NOutChain = 0,
        ArraySize = 0;
   char *ptr,
        *buffer,
        *outseq = NULL,
        *combseq = NULL,
        *align1,
        *align2,
        **seqs[2];
   BOOL DoneSEQRES[MAXCHAINS],
        DoneATOM[MAXCHAINS],
        DoInit;


   if((seqresSequence == NULL) || (atomSequence == NULL))
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outChains[i], atomChains[i]);
      }
      return(strdup(atomSequence));
   }
   
   
   /* Set flags to say we haven't handled the sequences yet             */
   for(i=0; i<MAXCHAINS; i++)
   {
      DoneSEQRES[i] = FALSE;
      DoneATOM[i]   = FALSE;
   }
   
   /* If the sequences and chains are identical just copy one of them
      and return
   */
   if(!strcmp(seqresSequence,atomSequence) &&
      !strcmp(seqresChains[0],atomChains[0])) /* FIXME! */
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outChains[i], seqresChains[i]);
      }
      return(strdup(atomSequence));
   }

   /* Create a temporary buffer to store a sequence                     */
   len1 = strlen(seqresSequence);
   len2 = strlen(atomSequence);
   if((buffer = (char *)malloc((1+MAX(len1, len2))
                               * sizeof(char)))==NULL)
   {
      return(NULL);
   }

   /* See how many chains there are and create arrays of char pointers
      to store the separate chains
   */
   nchain[0] = blCountchar(seqresSequence,   '*');
   if(len1 && seqresSequence[len1-1] != '*')
      nchain[0]++;
   nchain[1] = blCountchar(atomSequence, '*');
   if(len2 && atomSequence[len2-1] != '*')
      nchain[1]++;

   for(i=0; i<2; i++)
   {
      if((seqs[i] = (char **)malloc(nchain[i] * sizeof(char *)))==NULL)
         return(NULL);
   }

   /* Transfer the individual chains into the split arrays              */
   for(i=0; i<2; i++)
   {
      strcpy(buffer,((i==0)?seqresSequence:atomSequence));
      ptr = buffer;
      
      for(j=0; j<nchain[i]; j++)
      {
         TERMAT(ptr,'*');
         len = strlen(ptr);
         if((seqs[i][j] = (char *)malloc((1+len)*sizeof(char)))==NULL)
            return(NULL);
         strcpy(seqs[i][j], ptr);
         seqs[i][j][len] = '\0';

         ptr += strlen(ptr) + 1;
      }
   }


   /* Now align the sequences of the matching chains                    */
   for(i=0; i<nchain[0]; i++)
   {
      for(j=0; j<nchain[1]; j++)
      {
         if(CHAINMATCH(seqresChains[i], atomChains[j]))
         {
            DoneSEQRES[i] = TRUE;
            DoneATOM[j]   = TRUE;
            strcpy(outChains[NOutChain++], seqresChains[i]);
            
            if(!strcmp(seqs[0][i], seqs[1][j]))
            {
               /* If they are identical, copy to the output array
                  (+2 in the array size for * and \0)
               */
               DoInit = (ArraySize) ? FALSE : TRUE;
               ArraySize += strlen(seqs[0][i]) + 2;
               if((outseq = (char *)realloc(outseq, 
                                            ArraySize*sizeof(char)))
                  ==NULL)
               {
                  return(NULL);
               }
               if(DoInit)
                  outseq[0] = '\0';

               strcat(outseq,seqs[0][i]);
               strcat(outseq,"*");
            }
            else
            {
               /* The sequences are non-identical so we align them      */
               len1 = strlen(seqs[0][i]);
               len2 = strlen(seqs[1][j]);
               if((align1=(char *)malloc((len1+len2)*sizeof(char)))==NULL)
                  return(NULL);
               if((align2=(char *)malloc((len1+len2)*sizeof(char)))==NULL)
                  return(NULL);
            
               if(!blAlign(seqs[0][i], len1,
                           seqs[1][j], len2,
                           FALSE, TRUE, GAPPEN,
                           align1, align2, &align_len))
                  return(NULL);

               if((combseq=sCombineSequence(align1,align2,align_len, upper))
                  ==NULL)
                  return(NULL);

               free(align1); 
               free(align2);

               /* Allocate memory for the output sequence
                  (+2 in the array size for * and \0)
               */
               if(ArraySize==0)
               {
                  ArraySize = align_len + 2;
                  if((outseq = (char *)malloc(ArraySize*sizeof(char)))
                     ==NULL)
                     return(NULL);
                  outseq[0] = '\0';
               }
               else
               {
                  ArraySize += align_len + 2;
                  if((outseq = (char *)realloc(outseq,
                                               ArraySize*sizeof(char)))
                     ==NULL)
                     return(NULL);
                  outseq[ArraySize-1] = '\0';
               }

               strcat(outseq,combseq);
               strcat(outseq,"*");
               free(combseq);
            }
            break;
         }
      }
   }

   /* Add any chains from the ATOM records not yet handled              */
   for(i=0; i<nchain[1]; i++)
   {
      if(!DoneATOM[i])
      {
         /* Allocate memory for the output sequence
            (+2 in the array size for * and \0)
         */
         if(ArraySize==0)
         {
            ArraySize = strlen(seqs[1][i]) + 2;
            if((outseq = (char *)malloc(ArraySize*sizeof(char)))
               ==NULL)
               return(NULL);
            outseq[0] = '\0';
         }
         else
         {
            ArraySize += strlen(seqs[1][i]) + 2;
            if((outseq = (char *)realloc(outseq,
                                         ArraySize*sizeof(char)))
               ==NULL)
               return(NULL);
            outseq[ArraySize-1] = '\0';
         }
         
         strcat(outseq,seqs[1][i]);
         strcat(outseq,"*");
         strcpy(outChains[NOutChain++], atomChains[i]);
      }
   }

   /* Add any chains from the SEQRES records not yet handled.
      We issue a warning about these ones
   */
   if(!IgnoreSEQRES) /* 22.05.09 Added this check                       */
   {
      for(i=0; i<nchain[0]; i++)
      {
         if(!DoneSEQRES[i])
         {
            /* Allocate memory for the output sequence
               (+2 in the array size for * and \0)
            */
            if(ArraySize==0)
            {
               ArraySize = strlen(seqs[0][i]) + 2;
               if((outseq = (char *)malloc(ArraySize*sizeof(char)))
                  ==NULL)
                  return(NULL);
               outseq[0] = '\0';
            }
            else
            {
               ArraySize += strlen(seqs[0][i]) + 2;
               if((outseq = (char *)realloc(outseq,
                                            ArraySize*sizeof(char)))
                  ==NULL)
                  return(NULL);
               outseq[ArraySize-1] = '\0';
            }

            /* 22.05.09 Added this                                      */
            if(!upper)
            {
               LOWER(seqs[0][i]);
            }
            
            strcat(outseq,seqs[0][i]);
            strcat(outseq,"*");
            strcpy(outChains[NOutChain++], seqresChains[i]);
            
            if(!quiet)
            {
               if(label != NULL)
               {
                  fprintf(stderr,"Warning: Chain %s from SEQRES records \
not found in ATOM records%s%s\n",
                          seqresChains[i],
                          ((label[0])?" Label: ":""),
                          ((label[0])?label:""));
               }
               else
               {
                  fprintf(stderr,"Warning: Chain %s from SEQRES records \
not found in ATOM records\n",
                          seqresChains[i]);
               }
            }
         }
      }
   }

#ifdef DEBUG
   fprintf(stderr, "\n=============================================\n\n");
#endif

   /*** NEEDS TO FREE seqs[][] too                                    ***/
   free(buffer);         /*  11.06.15                                   */
   return(outseq);
}


/************************************************************************/
/*>static char *sCombineSequence(char *align1, char *align2, int align_len)
   ----------------------------------------------------------------
*//**
   \param[in]  *align1   ATOM   sequence aligned
   \param[in]  *align2   SEQRES sequence aligned
   \param[in]  align_len Alignment length
  
   Combine the information from the two sequences. Originally in pdb2pir

-  22.08.97 Original   By: ACRM
*/
static char *sCombineSequence(char *align1, char *align2, int align_len,
                              BOOL upper)
{
   static char *outseq = NULL;
   int         i;

   if((outseq=(char *)malloc((align_len+1) + sizeof(char)))==NULL)
      return(NULL);

#ifdef DEBUG
   fprintf(stderr,"\n----------------------------------------------\n");
   fprintf(stderr,"Alignment:\n");
   align1[align_len] = '\0';
   fprintf(stderr,"%s\n", align1);
   align2[align_len] = '\0';
   fprintf(stderr,"%s\n", align2);
#endif
   
   
   
   for(i=0; i<align_len; i++)
   {
      if((align1[i] == align2[i]) || (align1[i] == '-'))
      {
         outseq[i] = safetoupper(align2[i]);
      }
      else
      {
#ifdef DEBUG
         fprintf(stderr, "AL1: %c AL2: %c\n", align1[i], align2[i]);
#endif
         if(align2[i] == '-')
         {
            if(upper)
               outseq[i] = safetoupper(align1[i]);
            else
               outseq[i] = safetolower(align1[i]);
         }
         else
         {
            outseq[i] = safetoupper(align2[i]);
         }
      }
   }
   outseq[align_len] = '\0';
#ifdef DEBUG
   fprintf(stderr, "\n(Sequence now assembled)\n");
#endif
   return(outseq);
}



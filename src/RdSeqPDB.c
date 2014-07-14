/************************************************************************/
/**

   \file       RdSeqPDB.c
   
   \version    V1.2
   \date       07.07.14
   \brief      Read sequence from SEQRES records in a PDB file
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996-2014
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  14.10.96 Original   By: ACRM
-  V1.1  25.03.14 Added CHAINMATCH. By: CTP
-  V1.2  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
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
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
char **blReadSeqresPDB(FILE *fp, int *nchains);
static STRINGLIST *RdSeqRes(FILE *fp);


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
         seqs[chainnum][nres++] = throne(res[i]);
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


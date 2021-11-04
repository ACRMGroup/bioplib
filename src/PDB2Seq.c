/************************************************************************/
/**

   \file       PDB2Seq.c
   
   \version    V1.16
   \date       01.12.15
   \brief      Conversion from PDB to sequence and other sequence
               related routines
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 1993-2021
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
-  V1.0  29.09.92 Original
-  V1.1  07.06.93 Corrected allocation
-  V1.2  18.06.93 Handles multi-chains and skips NTER and CTER residues.
                  Added SplitSeq()
-  V1.3  09.07.93 SplitSeq() cleans up properly if allocation failed
-  V1.4  11.05.94 Added TrueSeqLen()
-  V1.5  13.05.94 Fixed bug in PDB2Seq().
                  Added KnownSeqLen().
-  V1.6  07.09.94 Fixed allocation bug in SplitSeq()
-  V1.7  19.07.95 Added check for ATOM records
-  V1.8  24.01.96 Fixed bug when no ATOM records in linked list
                  Returns a blank string
-  V1.9  26.08.97 Renamed DoPDB2Seq() with handling of Asx/Glx and
                  protein-only. Added macros to recreate the
                  old PDB2Seq() interface and similar new calls
-  V1.10 02.10.00 Added NoX option
-  V1.11 30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.12 10.06.05 Fixed bug - was undercounting by 1 for CA-only chains
-  V1.13 04.02.14 Use CHAINMATCH By: CTP
-  V1.14 07.07.14 Use bl prefix for functions By: CTP
-  V1.15 01.12.15 Added blDoPDB2SeqByChain()  By: ACRM
-  V1.16 03.11.21 HETATM PCA now handled as Q

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Conversions
   #FUNCTION  blDoPDB2Seq()
   malloc()'s an array containing the 1-letter sequence corresponding to
   an input PDB linked list.

   #FUNCTION  blDoPDB2SeqByChain()
   Creates a hash indexed by chain label containing the 1-letter code
   sequence from an input PDB linked list.
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "pdb.h"
#include "seq.h"


/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char *blDoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX)
   -------------------------------------------------------------------
*//**

   \param[in]     *pdb     PDB linked list
   \param[in]     DoAsxGlx Handle Asx and Glx as B and Z rather than X
   \param[in]     ProtOnly Don't do DNA/RNA; these simply don't get
                           done rather than being handled as X
   \param[in]     NoX      Skip amino acids which would be assigned as X
   \return                 Allocated character array containing sequence

   malloc()'s an array containing the 1-letter sequence corresponding to
   an input PDB linked list. Returns NULL if given a NULL parameter or
   memory allocation fails. Puts *'s in the sequence for multi-chains.

   This routine is normally called via the macro interfaces:
   PDB2Seq(pdb), PDB2SeqX(pdb), PDBProt2Seq(pdb), PDBProt2SeqX(pdb)
   Those with Prot in their names handle protein only; those with
   X handle Asx/Glx as B/Z rather than as X
   
-  29.09.92 Original    By: ACRM
-  07.06.93 Corrected allocation.
-  18.06.93 Handles multi-chains and skips NTER and CTER residues
-  13.05.94 Check for chain change *before* copy residue (!)
            (Bug reported by Bob MacCullum)
-  19.07.95 Added check for ATOM records
-  24.01.96 Returns blank string (rather than core dumping!) if the
            linked list contained no ATOM records
-  26.08.97 Changed to doPDB2Seq with extra parameters (DoAsxGlx & 
            ProtOnly). The old calling forms have now become macros
-  02.10.00 Added NoX
-  10.06.05 Changed the initialization of rescount, resnum, etc. so
            it correctly points to the first residue. This solves a
            bug with CA-only chains where it was undercounting by 1
-  04.02.14 Use CHAINMATCH By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  03.11.21 HETATM/PCA -> Q  By: ACRM
*/
char *blDoPDB2Seq(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, BOOL NoX)
{
   int   resnum,
         rescount,
         NBreak    = 0;
   char  insert,
         chain[8],
         *sequence = NULL;
   PDB   *p        = NULL;
   
   /* Sanity check                                                      */
   if(pdb==NULL) return(NULL);

   /* First step through the pdb linked list to see how many residues
      and chains.
      10.06.05 Fixed bug - was undercounting by one for CA-only chains
   */
   rescount = 1;
   resnum   = pdb->resnum;
   insert   = pdb->insert[0];
   strcpy(chain,pdb->chain);
   
   for(p=pdb->next; p!=NULL; NEXT(p))
   {
      if(p->resnum != resnum || p->insert[0] != insert)
      {
         if(strncmp(p->resnam,"NTER",4) && 
            strncmp(p->resnam,"CTER",4))
         {
            if(!strncmp(p->record_type,"ATOM  ",6) ||
               (!strncmp(p->record_type,"HETATM",6) &&
                !strncmp(p->resnam,"PCA",3)))
               rescount++;
         }
         
         resnum = p->resnum;
         insert = p->insert[0];
         
         /* Check for chain change                                      */
         if(!CHAINMATCH(chain,p->chain))
         {
            NBreak++;
            strcpy(chain,p->chain);
         }
      }
   }
   
   if(NBreak) rescount += NBreak;
   
   /* Allocate memory for the sequence array                            */
   sequence = malloc((rescount + 1) * sizeof(char));
   if(sequence == NULL) return(NULL);
   
   /* Step through the pdb linked list again, setting sequence array    */
   p = pdb;

   /* Skip an NTER residue                                              */
   /* 24.01.96 Added NULL check; occurs when no ATOM records present    */
/*
   while(p!=NULL && 
         (!strncmp(p->resnam,"NTER",4) || 
          strncmp(p->record_type,"ATOM  ",6))) 
      NEXT(p);
*/

   /* Skip non-residues at the start                                    */
   while(p!=NULL)
   {
      /* Jump out if this is an ATOM, but not NTER                      */
      if(!strncmp(p->record_type,"ATOM  ",6) &&
         strncmp(p->resnam,"NTER", 4))
         break;
      /* Jump out if this is a HETATM/PCA                               */
      if(!strncmp(p->record_type,"HETATM",6) &&
         !strncmp(p->resnam,"PCA",3))
         break;
      NEXT(p);
   }
   
   if(p==NULL)
   {
      sequence[0] = '\0';
      return(sequence);
   }
   
   sequence[0] = ((DoAsxGlx)?blThronex(p->resnam):blThrone(p->resnam));
   if((!ProtOnly) || (!gBioplibSeqNucleicAcid))
      rescount = 1;
   else
      rescount = 0;

   /* 02.10.00 Reset count if it's an X character and we are ignoring
      them
   */
   if(NoX && sequence[0] == 'X')   
      rescount = 0;
   
   resnum      = p->resnum;
   insert      = p->insert[0];
   strcpy(chain,p->chain);

   for(p=p->next; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->record_type,"ATOM  ",6) ||
         (!strncmp(p->record_type,"HETATM",6) &&
          !strncmp(p->resnam,"PCA ",3)))
      {
         if(p->resnum != resnum || p->insert[0] != insert)
         {
            /* Check for chain change                                   */
            if(!CHAINMATCH(chain,p->chain))
            {
               sequence[rescount++] = '*';
               strcpy(chain,p->chain);
            }
            
            /* 06.02.03 Fixed bug - was incrementing recount even when
               it was NTER/CTER
            */
            if(strncmp(p->resnam,"NTER",4) && strncmp(p->resnam,"CTER",4))
            {
               sequence[rescount] = ((DoAsxGlx) ? 
                                     blThronex(p->resnam):
                                     blThrone(p->resnam));
               if((!ProtOnly) || (!gBioplibSeqNucleicAcid))
                  rescount++;

               /* 02.10.00 Reset count if it's an X character and we are 
                  ignoring them
               */
               if(NoX && sequence[rescount-1] == 'X')   
                  rescount--;
            }

            resnum = p->resnum;
            insert = p->insert[0];
         }
      }
   }
   
   sequence[rescount] = '\0';
   
   return(sequence);
}


/************************************************************************/
/*>HASHTABLE *blDoPDB2SeqByChain(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, 
                                 BOOL NoX)
   ---------------------------------------------------------------------
*//**

   \param[in]     *pdb     PDB linked list
   \param[in]     DoAsxGlx Handle Asx and Glx as B and Z rather than X
   \param[in]     ProtOnly Don't do DNA/RNA; these simply don't get
                           done rather than being handled as X
   \param[in]     NoX      Skip amino acids which would be assigned as X
   \return                 A hash of 1-letter code sequences indexed by
                           chain label

   Reads sequence from ATOM records in 1-letter code, storing the
   results in a hash indexed by chain label.

   This routine is normally called via the macro interfaces:
   PDB2SeqByCHain(pdb), PDB2SeqXByCHain(pdb), PDBProt2SeqByChain(pdb),
   PDBProt2SeqXByChain(pdb)
   Those with Prot in their names handle protein only; those with
   X handle Asx/Glx as B/Z rather than as X
   
-  30.11.15 Original based on blDoPDB2Seq()    By: ACRM
-  03.11.21 HETATM/PCA -> Q
*/
HASHTABLE *blDoPDB2SeqByChain(PDB *pdb, BOOL DoAsxGlx, BOOL ProtOnly, 
                              BOOL NoX)
{
   int   lastresnum,
         nres       = 0,
         ArraySize  = ALLOCSIZE;
   char  lastinsert[blMAXCHAINLABEL],
         lastchain[blMAXCHAINLABEL],
         *sequence  = NULL;
   PDB   *p         = NULL;
   HASHTABLE *hash  = NULL;

   /* Ensure fist residue will be recognized as different               */
   lastresnum = (-1000);
   strcpy(lastinsert, "z");
   
   /* Sanity check                                                      */
   if(pdb==NULL) return(NULL);

   /* Initialize hash with 11 bins                                      */
   if((hash = blInitializeHash(11))==NULL)
      return(NULL);

   /* Initialize string to store the sequence                           */
   if((sequence=(char *)malloc(ArraySize * sizeof(char)))==NULL)
      return(NULL);
   sequence[0]  = '\0';

   /* Initialize chain label                                            */
   strncpy(lastchain, pdb->chain, blMAXCHAINLABEL);

   /* Step through the PDB linked list                                  */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* Only interested in ATOM records and HETATM/PCA                 */
      if(!strncmp(p->record_type, "ATOM  ", 6) ||
         (!strncmp(p->record_type, "HETATM", 6) &&
          !strncmp(p->resnam, "PCA ", 4)))
      {
         /* If chain has changed                                        */
         if(!CHAINMATCH(p->chain, lastchain))
         {
            if(blHashKeyDefined(hash, lastchain))
            {
               /*** TODO ***/
            }
            else
            {
               blSetHashValueString(hash, lastchain, sequence);
            }
            nres = 0;
         }

         /* If residue has changed                                      */
         if((p->resnum != lastresnum) ||
            !INSERTMATCH(p->insert, lastinsert) ||
            !CHAINMATCH(p->chain, lastchain))
         {
            if(strncmp(p->resnam,"NTER",4) && 
               strncmp(p->resnam,"CTER",4))
            {
               sequence[nres]=((DoAsxGlx) ? blThronex((p)->resnam) :
                               blThrone((p)->resnam));

               /* Increments count if it's not protein only or it's not 
                  a nucleic acid AND we aren't skipping Xs or it's not 
                  an X
               */
               if(((!ProtOnly) || (!gBioplibSeqNucleicAcid)) &&        
                  (!NoX || (sequence[nres] != 'X'))) 
                  nres++;

               /* Increase the array size if needed                     */
               if(nres >= ArraySize) 
               {
                  ArraySize += ALLOCSIZE;
                  if((sequence=(char *)realloc(sequence, ArraySize))
                     == NULL) 
                  {
                     blFreeHash(hash);
                     return(NULL);
                  }
               }
               
               /* Terminate the string                                  */
               sequence[nres] = '\0';
            }
            
            /* Updated information on last residue                      */
            lastresnum      = p->resnum;
            strcpy(lastinsert, p->insert);
            strcpy(lastchain,  p->chain);
         }
      }
   }
   
   if(nres)
      blSetHashValueString(hash, lastchain, sequence);

   if(sequence != NULL)
      free(sequence);
   
   return(hash);
}

#ifdef TEST
#include <stdio.h>
int main(int argc, char **argv)
{
   PDB *pdb;
   int natoms;
   char *seq;
   FILE *fp;
   HASHTABLE *seq2;
   
   if((fp = fopen(argv[1], "r"))!=NULL)
   {
      pdb=blReadPDB(fp, &natoms);
      fclose(fp);

      if((seq = blDoPDB2Seq(pdb, FALSE, FALSE, FALSE))!=NULL)
         printf("blDoPDB2Seq(): %s\n", seq);

      if((seq2 = blDoPDB2SeqByChain(pdb, FALSE, FALSE, FALSE))!=NULL)
      {
         char **chains = NULL;
         
         printf("blDoPDB2SeqByChain()\n");
         if((chains = blGetHashKeyList(seq2))!=NULL)
         {
            int i;
            
            for(i=0; chains[i]!=NULL; i++)
            {
               printf("%s : %s\n", chains[i],
                      blGetHashValueString(seq2, chains[i]));
            }
            
            blFreeHashKeyList(chains);
         }
      }
   }
   return(0);
}
#endif


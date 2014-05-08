/*************************************************************************

   Program:    
   File:       FindZonePDB.c
   
   Version:    V1.6
   Date:       07.05.14
   Function:   Routines for handling zones in PDB linked lists
   
   Copyright:  (c) SciTech Software 1993-2014
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
   V1.0  30.09.92 Original
   V1.1  16.06.93 Tidied for book. Mode now a char.
   V1.2  18.06.96 Added InPDBZone() from QTree program
   V1.3  19.09.96 Added InPDBZoneSpec()
   V1.5  20.03.14 Added blFindZonePDB() and deprecated InPDBZone() By: CTP
   V1.6  07.05.14 Moved FindZonePDB() to deprecated.h By: CTP

*************************************************************************/

/* Includes
*/
#include "SysDefs.h"
#include "MathType.h"

#include "pdb.h"
#include "macros.h"


/************************************************************************/
/*>BOOL blFindZonePDB(PDB *pdb, int start, char *startinsert, int stop, 
                      char *stopinsert, char *chain, int mode, 
                      PDB **pdb_start, PDB **pdb_stop)
   -------------------------------------------------------------
   Input:   PDB   *pdb         PDB linked list
            int   start        Resnum of start of zone
            char  *startinsert Insert code for start of zone
            int   stop         Resnum of end of zone
            char  *stopinsert  Insert code for end of zone
            char  *chain       Chain name
            int   mode         ZONE_MODE_RESNUM:     Use PDB residue 
                                                     numbers/chain
                               ZONE_MODE_SEQUENTIAL: Use sequential 
                                                     numbering
   Output:  PDB   **pdb_start  Start of zone
            PDB   **pdb_stop   End of zone
   Returns: BOOL               OK?

   Finds pointers to the start and end of a zone in a PDB linked list. The
   end is the atom *after* the specified zone

   20.03.14 Function based on FindZonePDB() but takes string for chain 
            label and insert instead of single char. 
            Function finds end of chain properly if stop set to -999.
            By: CTP
*/
BOOL blFindZonePDB(PDB   *pdb,
                   int   start,
                   char  *startinsert,
                   int   stop,
                   char  *stopinsert,
                   char  *chain,
                   int   mode,
                   PDB   **pdb_start,
                   PDB   **pdb_stop)
{
   PDB   *p;
   int   rescount,
         resnum;
   BOOL  InStop     = FALSE,
         FoundChain = FALSE;
   char  insert[8];
   
   /* To start, we don't know where either are                          */
   *pdb_start = NULL;
   *pdb_stop  = NULL;
   
   /* If both start and stop are -999, then the whole structure (or a 
      whole chain) is being specified
   */
   if((start == (-999)) && (stop == (-999)))
   {
      if(CHAINMATCH(chain," "))           /* Whole structure            */
      {
         *pdb_start = pdb;
         *pdb_stop  = NULL;
         return(TRUE);
      }
      else                                /* An individual chain        */
      {
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(CHAINMATCH(p->chain,chain))
            {
               if(*pdb_start == NULL)
               {
                  *pdb_start = p;
               }
            }
            else if(*pdb_start != NULL)   /* We've aleady got the start */
            {
               *pdb_stop = p;
               return(TRUE);
            }
         }
         if(*pdb_start==NULL) 
            return(FALSE);                /* Chain not found            */
         else
            return(TRUE);
      }
   }
   
   /* Handle one end of a zone being set to -999                        */
   if(start == -999) *pdb_start = pdb;
   if(stop  == -999) *pdb_stop  = NULL;
   
   /* If either end is still undefined                                  */
   if(*pdb_start == NULL || *pdb_stop == NULL)
   {
      /* Search reference structure for start and end of zone           */
      rescount = 1;
      resnum   = pdb->resnum;
      strcpy(insert,pdb->insert);
      InStop   = FALSE;
      FoundChain = FALSE;
   
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(mode == ZONE_MODE_RESNUM)
         {
            if(CHAINMATCH(chain," ") || CHAINMATCH(chain,p->chain))
            {  /* We are in the correct chain                           */
               FoundChain = TRUE;

               /* If start undefined, see if residue matches            */
               if(*pdb_start == NULL)
               {
                  if((p->resnum == start) && 
                     !strcmp(p->insert,startinsert))
                     *pdb_start = p;
               }
               
               /* If stop undefined, then find the following residue    */
               if(*pdb_stop == NULL)
               {
                  /* See if we have just moved out of the stop residue.
                     If so, set the stop position and return
                  */
                  if(InStop && 
                     (p->resnum != stop || strcmp(p->insert,stopinsert)))
                  {
                     *pdb_stop = p;
                     return((*pdb_start==NULL)?FALSE:TRUE);
                  }

                  /* Residue matches, so set flag to say we're in the
                     last residue of the zone.
                  */
                  if((p->resnum == stop) &&
                     !strcmp(p->insert,stopinsert))
                     InStop = TRUE;
               }
               if(*pdb_start != NULL && *pdb_stop != NULL)  /*Found both*/
                  break;
            }
            else if(InStop)
            {
               /* We will get here if InStop has been set without having 
                  found the start of the next residue. This will occur
                  if the last residue of a zone was also the last 
                  residue of a chain, since the chain name will now have
                  changed.
                  We just set *pdb_stop to this pointer and return.
               */
               *pdb_stop = p;
               return((*pdb_start==NULL)?FALSE:TRUE);
            }
            else if(stop == -999 && FoundChain == TRUE)
            {
               /* We get here if stop is set for the end of a chain and
                  a the start of new chain has been found.      By: CTP */
               *pdb_stop = p;
               return((*pdb_start==NULL)?FALSE:TRUE);               
            }
         }  /* End of ZONE_MODE_RESNUM                                  */
         else if(mode == ZONE_MODE_SEQUENTIAL)
         {
            /* Correct the residue count                                */
            if(p->resnum != resnum || strcmp(p->insert,insert))
            {
               rescount++;
               resnum = p->resnum;
               strcpy(insert,p->insert);
            }
            
            if(*pdb_start == NULL)           /* Identify zone start     */
               if(rescount == start) *pdb_start = p;
            if(*pdb_stop == NULL)            /* Identify zone stop      */
            {
               if(InStop && rescount != stop) *pdb_stop = p;
               
               if(rescount == stop) InStop = TRUE;
            }
            if(*pdb_start != NULL && *pdb_stop != NULL)   /* Found both */
               break;
         }
         else
         {
            printf("   Error==> CreateFitAtoms(): Internal confusion!\n");
         }
      }  /* End of loop through PDB linked list                         */
   }  /* End of if() one pointer undefined                              */

   /* Check start of range has been found and return                    */
   return((*pdb_start==NULL)?FALSE:TRUE);
}


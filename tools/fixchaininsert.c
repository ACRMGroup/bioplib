/************************************************************************/
/**

   \file       fixchaininsert.c
   
   \version    V1.2
   \date       22.07.14
   \brief      Look for examples where the chain name is in the insert
               column (1mfa)
   
   \copyright  (c) UCL, Dr. Andrew C. R. Martin 1997-2014
   \author     Dr. Andrew C. R. Martin
   \par
               Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
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
-  V1.0  25.07.97 Original    By: ACRM
-  V1.1  04.11.08 Changed pdb.junk to pdb.record_type
-  V1.2  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define CADISTSQ 16.0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
void CheckChainInsert(PDB *pdb);
void CheckChain(PDB *start, PDB *end);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for neighbour counting

-  05.07.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   int  natoms;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in,&natoms)) != NULL)
         {
            CheckChainInsert(pdb);
            blWritePDB(out, pdb);
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \return                     Success?

   Parse the command line
   
-  25.07.97 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void CheckChainInsert(PDB *pdb)
   -------------------------------
*//**

   Splits the PDB linked list into chains and calls CheckChain() on each
   to look for the chain name in the insert column

-  25.07.97 Original   By: ACRM
-  03.11.08 Fixed pdb.junk to pdb.record_type
*/
void CheckChainInsert(PDB *pdb)
{
   PDB *p,
       *start,          /* Start of chain                               */
       *rstart,         /* Start of current residue                     */
       *prev,
       *CAPrev = NULL;
   
   for(rstart=start=p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    != rstart->resnum)    ||
         (p->insert[0] != rstart->insert[0]) ||
         (p->chain[0]  != rstart->chain[0]))
      {
         rstart = p;
      }
      
      /* Check for move from ATOM to HETATM: treat as new chain         */
      if((prev!=NULL) && strncmp(p->record_type,prev->record_type,6))
      {
         rstart = p;
         CheckChain(start,rstart);
         start=rstart;
         CAPrev = NULL;
      }
      else if(!strncmp(p->record_type,"ATOM  ",6) && 
              !strncmp(p->atnam,"CA  ",4))
      {
         if(CAPrev!=NULL)
         {
            if(DISTSQ(p,CAPrev) > CADISTSQ)
            {
               CheckChain(start,rstart);
               start  = rstart;
            }
         }
         CAPrev = p;
      }
      prev = p;
   }
   /* Do the final chain                                                */
   CheckChain(start,NULL);
}


/************************************************************************/
/*>void CheckChain(PDB *start, PDB *end)
   -------------------------------------
*//**

   Does the actual work of looking for the chain name in the insert column

-  25.07.97 Original   By: ACRM
*/
void CheckChain(PDB *start, PDB *end)
{
   PDB *p;
   BOOL Bad = TRUE;
   
   for(p=start; p!=end; NEXT(p))
   {
      if((p->insert[0] == ' ') || (p->chain[0] != ' '))
      {
         Bad = FALSE;
         break;
      }
   }
   
   if(Bad)
   {
      for(p=start; p!=end; NEXT(p))
      {
         p->chain[0] = p->insert[0];
         p->insert[0] = ' ';
      }
   }
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  25.07.97 Original    By: ACRM
-  22.07.14 V1.2 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nFixChainInsert V1.2 (c) 1997-2014, Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: fixchaininsert [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nRuns through a PDB file and sees whether all atoms \
in a chain have an \n");
   fprintf(stderr,"insert code. If they do and the chain name is \
blank, moves the\n");
   fprintf(stderr,"insert code to the chain name field.\n");

   fprintf(stderr,"\nI/O through stdin/stdout if files not \
specified.\n\n");
}

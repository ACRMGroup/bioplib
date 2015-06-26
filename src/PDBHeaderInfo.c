/************************************************************************/
/**

   \file       PDBHeaderInfo.c
   
   \version    V1.4
   \date       11.06.15

   \brief      Get misc header info from PDB header
   
   \copyright  (c) UCL / Dr. Andrew C.R. Martin, 2015
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

   See documentation for details

**************************************************************************

   Revision History:
   =================
-  V1.0  26.03.15 Original
-  V1.1  13.05.15 blGetTitleWholePDB() returns NULL if no title line in
                  header. Fix for blGetSpeciesWholePDBChain()
                  Added blGetCompoundWholePDBMolID() and 
                  blGetSpeciesWholePDBMolID() By: CTP
-  V1.2  04.06.15 Fixed bug in dealing with compounds where the referenced
                  chains span more than one line
-  V1.3  09.06.15 Merged changes from CTP and ACRM. 
                  Updated blGetTitleWholePDB()  By: CTP
-  V1.4  11.06.15 Added blGetSeqresAsStringWholePDB(),
                  blGetModresWholePDB() and blFindOriginalResType() 
                  By: ACRM
-  V1.5  26.06.15 Added blGetBiomoleculeWholePDB() and blFreeBiomolecule()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP File IO

   #FUNCTION blGetHeaderWholePDB()
   Obtains the data in the HEADER record from WHOLEPDB info

   #FUNCTION blGetTitleWholePDB()
   Obtains the title information from WHOLEPDB info

   #FUNCTION blGetCompoundWholePDBChain()
   Obtains the compound data for a specified chain from WHOLEPDB info

   #FUNCTION blFindMolID()
   Finds the MOL_ID for a specified chain

   #FUNCTION blGetSpeciesWholePDBChain()
   Obtains the species data for a specified chain from WHOLEPDB info
   
   #FUNCTION blGetCompoundWholePDBCMolID()
   Obtains the compound data for a specified MOL_ID from WHOLEPDB info
   
   #FUNCTION blGetSpeciesWholePDBMolID()
   Obtains the species data for a specified MOL_ID from WHOLEPDB info

   #FUNCTION blGetSeqresAsStringWholePDB()
   Obtain the sequence from the SEQRES records

   #FUNCTION blGetModresWholePDB()
   Obtain the MODRES data

   #FUNCTION blFindOriginalResType()
   Find the original residue type for a modified residue from MODRES
   data

   #FUNCTION blGetBiomoleculeWholePDB()
   Obtain the biomolecule data

   #FUNCTION blFreeBiomolecule()
   Free the biomolecule data
*/

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "pdb.h"
#include "seq.h"
#include "macros.h"
#include "fsscanf.h"
#include "general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXWORD     8
#define MAXBUFF   160
#define ALLOCSIZE  80

#define CLEAR_BIOMOL(b)                         \
   do{                                          \
   (b)->details         = NULL;                 \
   (b)->numBiomolecules = 0;                    \
   (b)->biomolNumber    = 0;                    \
   (b)->authorUnit[0]   = '\0';                 \
   (b)->softwareUnit[0] = '\0';                 \
   (b)->chains          = NULL;                 \
   (b)->biomt           = NULL;                 \
   }  while(0)
   
/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void checkRemark300(BIOMOLECULE *biomolecule, STRINGLIST *s, 
                           int *SkipStandardRemark);
static void checkRemark350(BIOMOLECULE *biomolecule, STRINGLIST *s);

/************************************************************************/
/*>BOOL blGetHeaderWholePDB(WHOLEPDB *wpdb, 
                            char *header,  int maxheader,
                            char *date,    int maxdate,
                            char *pdbcode, int maxcode)
   ---------------------------------------------------
*//**

   \param[in]     *wpdb        WHOLEPDB structure pointer
   \param[out]    *header      String containing header text
   \param[in]     maxheader    Max length for storing header
   \param[out]    *date        Date string
   \param[in]     maxdate      Max length for storing date
   \param[out]    *pdbcode     PDB code
   \param[in]     maxcode      Max length for storing PDB code
   \return                     TRUE:  Found HEADER
                               FALSE: Didn't find HEADER

   Obtains information from the PDB HEADER record

-  26.03.15  Original   By: ACRM
*/
BOOL blGetHeaderWholePDB(WHOLEPDB *wpdb, 
                         char *header,  int maxheader,
                         char *date,    int maxdate,
                         char *pdbcode, int maxcode)
{
   STRINGLIST *s;
   int        i;
   BOOL       retval = FALSE;

   /* Blank all the strings                                             */
   for(i=0; i<maxheader; i++) header[i]  = '\0';
   for(i=0; i<maxdate;   i++) date[i]    = '\0';
   for(i=0; i<maxcode;   i++) pdbcode[i] = '\0';

   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "HEADER", 6))
      {
         retval = TRUE;
         strncpy(header,  s->string+10, MIN(40, maxheader));
         KILLTRAILSPACES(header);
         strncpy(date,    s->string+50, MIN( 9, maxdate));
         strncpy(pdbcode, s->string+62, MIN( 4, maxcode));
         break;
      }
   }

   return(retval);
}


/************************************************************************/
/*>char *blGetTitleWholePDB(WHOLEPDB *wpdb)
   ----------------------------------------
*//**
   \param[in]    *wpdb    WHOLEPDB structure
   \return                Tit;le from PDB file (malloc()'d)

   Extracts the title from a PDB file malloc()ing a string in which to
   store the data. This must be freed by user code

-  28.04.15 Original   By: ACRM
-  11.05.15 Return NULL if TITLE line absent. By: CTP
-  09.06.15 Add columns 11 to 80 to title string for both start and 
            continuation lines. By: CTP
*/
char *blGetTitleWholePDB(WHOLEPDB *wpdb)
{
   char       *title = NULL,
              *cleanTitle = NULL;
   STRINGLIST *s;
   BOOL       inTitle = FALSE;

   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "TITLE ", 6))
      {
         char buffer[MAXPDBANNOTATION];
         strcpy(buffer, s->string);
         TERMINATE(buffer);

         /* append cols 11-80 to title string                           */
         title = blStrcatalloc(title, buffer+10);

         if(title == NULL)
            return(NULL);
      }
      else if(inTitle)
      {
         break;
      }
   }

   /* title line not found                                              */
   if(title == NULL)
      return(NULL);

   cleanTitle = blCollapseSpaces(title);
   free(title);
   KILLTRAILSPACES(cleanTitle);
   
   return(cleanTitle);
}

/************************************************************************/
/*>static STRINGLIST *FindNextMolIDRecord(STRINGLIST *start, char *type)
   ---------------------------------------------------------------------
*//**
   \param[in]   *start    Start of stringlist containing header
   \param[in]   *type     Type of header record - COMPND or SOURCE
   \return                Pointer to start of the next molecule ID in
                          the appropriate header records

   Find the next MOL_ID within the specified header record type (COMPND
   or SOURCE)

   28.04.15  Original   By: ACRM
*/
static STRINGLIST *FindNextMolIDRecord(STRINGLIST *start, char *type)
{
   STRINGLIST *s;
   
   if(start==NULL)
      return(NULL);
   
   for(s=start->next; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, type, 6))
      {
         if(strstr(s->string, "MOL_ID:"))
            return(s);
      }
   }
   return(NULL);
}


/************************************************************************/
/*>static BOOL ExtractField(STRINGLIST *molidStart, 
                            STRINGLIST *molidStop, char *data,
                            char *type, char *field)
   ------------------------------------------------------------
*//**
   \param[in]   *molidStart    Start of a set of header records 
   \param[in]   *molidStop     Start of next set of headers (or NULL)
   \param[out]  *data          Storage for extracted string
   \param[in]   *type          Record type (COMPND or SOURCE)
   \param[in]   *field         Sub-record field of interest
   \return                     Success

   Extracts data for a field from a COMPND or SOURCE record. The field 
   data after the field specfication and is terminated by a ;

   Returns FALSE if field not found.

   28.04.15  Original   By: ACRM
*/
static BOOL ExtractField(STRINGLIST *molidStart, STRINGLIST *molidStop,
                         char *data, char *type, char *field)
{
   STRINGLIST *s;
   BOOL       GotField = FALSE;
   char       *chp,
              buffer[MAXPDBANNOTATION];

   data[0] = '\0';

   for(s=molidStart; s!=molidStop; NEXT(s))
   {
      if(strncmp(s->string, type, 6))
         break;

      chp = NULL;

      if(GotField && isdigit(s->string[9]))
      {
         /* We have found the field already on previous line and this is
            marked as a continuation line            
         */
         chp = s->string+10;
      }
      else
      {
         /* Look for this field                                         */
         if((chp=strstr(s->string, field))!=NULL)
         {
            GotField = TRUE;
            /* Step over the field name                                 */
            chp += strlen(field);
            if(*chp == ' ')
               chp++;
         }
      }
      
      if(GotField && (chp != NULL))
      {
         /* Copy into the buffer                                        */
         strncpy(buffer, chp, MAXPDBANNOTATION);
         /* Remove spaces                                               */
         TERMINATE(buffer);
         KILLTRAILSPACES(buffer);
         /* Add to output data                                          */
         blStrncat(data, buffer, MAXPDBANNOTATION);

         /* Exit if the string contains a ;                             */
         if((chp=strchr(data, ';'))!=NULL)
         {
            *chp = '\0';
            return(TRUE);
         }
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>BOOL blGetCompoundWholePDBChain(WHOLEPDB *wpdb, char *chain, 
                                   COMPND *compnd)
   ------------------------------------------------------------
*//**
   \param[in]    *wpdb    WHOLEPDB structure
   \param[in]    *chain   Chain label of interest
   \param[out]   *compnd  Data from the COMPND records
   \return       BOOL     Success

   Extracts the COMPND data for a specified chain. Returns FALSE if the
   chain isn't found

   28.04.15 Original   By: ACRM
   04.06.15  Modified to use the ExtractField() routine instead of
             duplicating code here. Fixes a bug in dealing with compounds
             where the referenced chains span more than one line
*/
BOOL blGetCompoundWholePDBChain(WHOLEPDB *wpdb, char *chain, 
                                COMPND *compnd)
{
   STRINGLIST *molidFirst,
              *molidStart,
              *molidStop;
   int        molid;

   compnd->molid         = 0;
   compnd->molecule[0]   = '\0';
   compnd->chain[0]      = '\0';
   compnd->fragment[0]   = '\0';
   compnd->synonym[0]    = '\0';
   compnd->ec[0]         = '\0';
   compnd->engineered[0] = '\0';
   compnd->mutation[0]   = '\0';
   compnd->other[0]      = '\0';

#ifdef DEBUG
   molid = blFindMolID(wpdb, chain);
   fprintf(stderr,"DEBUG: Chain %s molid %d\n", chain, molid);
   if(molid == 0)
      return(FALSE);
#else
   if((molid = blFindMolID(wpdb, chain)) == 0)
      return(FALSE);
#endif
   
   molidFirst = FindNextMolIDRecord(wpdb->header, "COMPND");

   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      char buffer[MAXPDBANNOTATION];
      int  thisMolid = 0;

      molidStop  = FindNextMolIDRecord(molidStart, "COMPND");

      ExtractField(molidStart, molidStop,
                   buffer,             "COMPND", "MOL_ID:");
      sscanf(buffer,"%d", &thisMolid);

      if(thisMolid == molid)
      {
         ExtractField(molidStart, molidStop,
                      compnd->molecule,   "COMPND","MOLECULE:");
         ExtractField(molidStart, molidStop,
                      compnd->chain,      "COMPND", "CHAIN:");
         ExtractField(molidStart, molidStop,
                      compnd->fragment,   "COMPND", "FRAGMENT:");
         ExtractField(molidStart, molidStop,
                      compnd->synonym,    "COMPND", "SYNONYM:");
         ExtractField(molidStart, molidStop,
                      compnd->ec,         "COMPND", "EC:");
         ExtractField(molidStart, molidStop,
                      compnd->engineered, "COMPND", "ENGINEERED:");
         ExtractField(molidStart, molidStop,
                      compnd->mutation,   "COMPND", "MUTATION:");
         ExtractField(molidStart, molidStop,
                      compnd->other,      "COMPND", "OTHER:");
         ExtractField(molidStart, molidStop,
                      buffer,             "COMPND", "MOL_ID:");
         sscanf(buffer,"%d", &(compnd->molid));
         return(TRUE);
      }
   }

   return(FALSE);
}


/************************************************************************/
/*>int blFindMolID(WHOLEPDB *wpdb, char *chain)
   --------------------------------------------
*//**
   \param[in]    *wpdb   WHOLEPDB structure
   \param[in]    *chain  Chain label
   \return               MOL_ID or 0 if chain not found

   Finds the MOL_ID for a specified chain

   28.04.15  Original   By: ACRM
   04.06.15  Modified to use the ExtractField() routine instead of
             duplicating code here. Fixes a bug in dealing with compounds
             where the referenced chains span more than one line
*/
int blFindMolID(WHOLEPDB *wpdb, char *chain)
{
   STRINGLIST *molidFirst,
              *molidStart,
              *molidStop;

   molidFirst = FindNextMolIDRecord(wpdb->header, "COMPND");

   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      char buffer[MAXPDBANNOTATION],
           *chp,
           word[MAXWORD];
      
      molidStop  = FindNextMolIDRecord(molidStart, "COMPND");
      ExtractField(molidStart, molidStop, buffer, "COMPND", "CHAIN:");
      
      /* Check the chains to see if our chain is there            */
      chp = buffer;
      
      do {
         int molid = 0;
               
         chp=blGetWord(chp, word, MAXWORD);
         if(!strcmp(word, chain))
         {
            ExtractField(molidStart, molidStop,
                         buffer,             "COMPND", "MOL_ID:");
            sscanf(buffer,"%d", &molid);
            return(molid);
         }
      }  while(chp!=NULL);
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL blGetSpeciesWholePDBChain(WHOLEPDB *wpdb, char *chain,
                                  PDBSOURCE *source)
   -----------------------------------------------------------
*//**
   \param[in]    *wpdb    WHOLEPDB structure
   \param[in]    *chain   Chain label
   \param[out]   *source  SOURCE information for chain
   \return                Success (chain found?)

   Extracts the SOURCE data for a specified chain
*/
BOOL blGetSpeciesWholePDBChain(WHOLEPDB *wpdb, char *chain,
                               PDBSOURCE *source)
{
   STRINGLIST *s,
              *molidFirst = NULL,
              *molidStart = NULL,
              *molidStop  = NULL;
   int        molid    = 0;

   source->scientificName[0] = '\0';
   source->commonName[0]     = '\0';
   source->strain[0]         = '\0';
   source->taxid             = 0;

   if((molid = blFindMolID(wpdb, chain)) == 0)
      return(FALSE);
   
   molidFirst = FindNextMolIDRecord(wpdb->header, "SOURCE");

   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      molidStop  = FindNextMolIDRecord(molidStart, "SOURCE");
      for(s=molidStart; s!=molidStop; NEXT(s))
      {
         char buffer[MAXPDBANNOTATION];
         int  thisMolid = 0;

         ExtractField(molidStart, molidStop,
                      buffer,             "SOURCE", "MOL_ID:");
         sscanf(buffer,"%d", &thisMolid);

         if(thisMolid == molid)
         {
            ExtractField(molidStart, molidStop,
                         source->scientificName, "SOURCE","ORGANISM_SCIENTIFIC:");
            ExtractField(molidStart, molidStop,
                         source->commonName,     "SOURCE", "ORGANISM_COMMON:");
            ExtractField(molidStart, molidStop,
                         source->strain,         "SOURCE", "STRAIN:");
            ExtractField(molidStart, molidStop,
                         buffer,                 "SOURCE", "ORGANISM_TAXID:");
            sscanf(buffer,"%d",&source->taxid);
            return(TRUE);
         }
      }
   }

   return(FALSE);
}

/************************************************************************/
/*>BOOL blGetCompoundWholePDBMolID(WHOLEPDB *wpdb, int molid, 
                                   COMPND *compnd)
   ------------------------------------------------------------
*//**
   \param[in]    *wpdb    WHOLEPDB structure
   \param[in]    molid    MOL_ID from PDB-format COMPND record.
   \param[out]   *compnd  Data from the COMPND records
   \return       BOOL     Success

   Extracts the COMPND data for a specified MOL_ID. Returns FALSE if the
   MOL_ID isn't found

-  13.05.15 Original based on blGetCompoundWholePDBChain().  By: CTP
*/
BOOL blGetCompoundWholePDBMolID(WHOLEPDB *wpdb, int molid, 
                                 COMPND *compnd)
{
   STRINGLIST *molidFirst,
              *molidStart,
              *molidStop,
              *s;

   /* reset compnd                                                      */
   compnd->molid         = 0;
   compnd->molecule[0]   = '\0';
   compnd->chain[0]      = '\0';
   compnd->fragment[0]   = '\0';
   compnd->synonym[0]    = '\0';
   compnd->ec[0]         = '\0';
   compnd->engineered[0] = '\0';
   compnd->mutation[0]   = '\0';
   compnd->other[0]      = '\0';

   /* find start of compnd records                                      */
   molidFirst = FindNextMolIDRecord(wpdb->header, "COMPND");

   /* get compound record                                               */
   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      molidStop  = FindNextMolIDRecord(molidStart, "COMPND");
      for(s=molidStart; s!=molidStop; NEXT(s))
      {
         char buffer[MAXPDBANNOTATION];
         int  thisMolid = 0;

         ExtractField(molidStart, molidStop,
                      buffer,             "COMPND", "MOL_ID:");
         sscanf(buffer,"%d", &thisMolid);

         if(thisMolid == molid)
         {
            ExtractField(molidStart, molidStop,
                         compnd->molecule,   "COMPND","MOLECULE:");
            ExtractField(molidStart, molidStop,
                         compnd->chain,      "COMPND", "CHAIN:");
            ExtractField(molidStart, molidStop,
                         compnd->fragment,   "COMPND", "FRAGMENT:");
            ExtractField(molidStart, molidStop,
                         compnd->synonym,    "COMPND", "SYNONYM:");
            ExtractField(molidStart, molidStop,
                         compnd->ec,         "COMPND", "EC:");
            ExtractField(molidStart, molidStop,
                         compnd->engineered, "COMPND", "ENGINEERED:");
            ExtractField(molidStart, molidStop,
                         compnd->mutation,   "COMPND", "MUTATION:");
            ExtractField(molidStart, molidStop,
                         compnd->other,      "COMPND", "OTHER_DETAILS:");
            ExtractField(molidStart, molidStop,
                         buffer,             "COMPND", "MOL_ID:");
            sscanf(buffer,"%d", &(compnd->molid));
            return(TRUE);
         }
      }
   }

   return(FALSE);
}

/************************************************************************/
/*>BOOL blGetSpeciesWholePDBMolID(WHOLEPDB *wpdb, int molid,
                                  PDBSOURCE *source)
   -----------------------------------------------------------
*//**
   \param[in]    *wpdb    WHOLEPDB structure
   \param[in]    molid    MOL_ID from PDB-format SOURCE record.
   \param[out]   *source  SOURCE information for chain
   \return                Success (chain found?)

   Extracts the SOURCE data for a specified MOL_ID. Returns FALSE if not 
   found.
   
   12.05.15 Original based on blGetSpeciesWholePDBChain().  By: CTP
*/
BOOL blGetSpeciesWholePDBMolID(WHOLEPDB *wpdb, int molid,
                                PDBSOURCE *source)
{
   STRINGLIST *s,
              *molidFirst = NULL,
              *molidStart = NULL,
              *molidStop  = NULL;

   /* reset source */
   source->scientificName[0] = '\0';
   source->commonName[0]     = '\0';
   source->strain[0]         = '\0';
   source->taxid             = 0;

   /* find start of source records */
   molidFirst = FindNextMolIDRecord(wpdb->header, "SOURCE");


   /* get source record */
   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      molidStop  = FindNextMolIDRecord(molidStart, "SOURCE");
      for(s=molidStart; s!=molidStop; NEXT(s))
      {
         char buffer[MAXPDBANNOTATION];
         int  thisMolid = 0;

         ExtractField(molidStart, molidStop,
                      buffer,             "SOURCE", "MOL_ID:");
         sscanf(buffer,"%d", &thisMolid);

         if(thisMolid == molid)
         {
            ExtractField(molidStart, molidStop,
                         source->scientificName, "SOURCE","ORGANISM_SCIENTIFIC:");
            ExtractField(molidStart, molidStop,
                         source->commonName,     "SOURCE", "ORGANISM_COMMON:");
            ExtractField(molidStart, molidStop,
                         source->strain,         "SOURCE", "STRAIN:");
            ExtractField(molidStart, molidStop,
                         buffer,                 "SOURCE", "ORGANISM_TAXID:");
            sscanf(buffer,"%d",&source->taxid);
            return(TRUE);
         }
      }
   }

   return(FALSE);
}


/************************************************************************/
/*>char *blGetSeqresAsStringWholePDB(WHOLEPDB *wpdb, char **chains, 
                                     MODRES *modres, BOOL doNucleic)
   -----------------------------------------------------------------
*//**
   \param[in]  *wpdb      Pointer to whole PDB structure
   \param[out] **chains   Chain labels for the chains - may also be
                          set to NULL if you don't want to record
                          chain labels
   \param[in]  *modres    Linked list of MODRES information. May be 
                          NULL if you don't want to translate 
                          non-standard amino acids.
   \param[in]  doNucleic  Read sequence for nucleic acid chains
   \return                malloc()'d Sequence from SEQRES, chains 
                          separated by a *

   Reads sequence from SEQRES records into a character string in 1-letter
   code. Chains are terminated by * characters.

-  21.08.97 Original   by: ACRM
-  22.08.97 Added chains parameter
-  26.08.97 No longer reads DNA/RNA
-  07.03.07 Added code to check for modified amino acids
            Now reads from wpdb rather than from the file
-  07.11.14 Initialize lastchain
-  11.06.15 Moved to bioplib - doNucleic is now a paramater instead of
            a global; chains is now an array of strings
-  12.06.15 Frees memory and returns NULL if no SEQRES found
*/
char *blGetSeqresAsStringWholePDB(WHOLEPDB *wpdb, char **chains, 
                                  MODRES *modres, BOOL doNucleic)
{
   static char *sequence = NULL;
   char        buffer[MAXBUFF],
               chain[blMAXCHAINLABEL],
               lastchain[blMAXCHAINLABEL],
               seq3[13][4];
   int         i,
               nchain    = 0,
               nres      = 0,
               ArraySize = ALLOCSIZE;
   BOOL        AddStar   = FALSE;
   STRINGLIST  *s;
   
   lastchain[0] = '\0';

   if((sequence=(char *)malloc(ArraySize * sizeof(char)))==NULL)
   {
      return(NULL);
   }
   sequence[0] = '\0';
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      strncpy(buffer, s->string, MAXBUFF);
      TERMINATE(buffer);
      if(!strncmp(buffer,"SEQRES",6))
      {
         fsscanf(buffer, 
                 "%11x%1s%7x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x\
%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s",
                 chain,
                 seq3[0],  seq3[1],  seq3[2],  seq3[3],  seq3[4], 
                 seq3[5],  seq3[6],  seq3[7],  seq3[8],  seq3[9],
                 seq3[10], seq3[11], seq3[12]);

         if((nres == 0) && !AddStar)
         {
            /* This is the first line so we set the lastchain           */
            strcpy(lastchain, chain);
            if(chains!=NULL)
               strncpy(chains[nchain++], chain, blMAXCHAINLABEL);
         }
         else if(nres+15 >= ArraySize)
         {
            /* Allocate more space if needed                            */
            ArraySize += ALLOCSIZE;
            if((sequence=(char *)realloc((void *)sequence, 
                                         ArraySize*sizeof(char)))
               == NULL)
            {
               return(NULL);
            }
         }

         if(!CHAINMATCH(chain, lastchain))
         {
            sequence[nres++] = '*';
            strcpy(lastchain, chain);
            if(chains!=NULL)
               strncpy(chains[nchain++], chain, blMAXCHAINLABEL);
         }

         for(i=0; i<13; i++)
         {
            AddStar = TRUE;
            if(!strncmp(seq3[i],"   ",3))
               break;
            sequence[nres] = blThronex(seq3[i]);

            /* 07.03.07 Added code to check for modified amino acids    */
            if(sequence[nres] == 'X')
            {
               char tmpthree[8];
               if(modres != NULL)   /* 11.06.15                         */
               {
                  blFindOriginalResType(seq3[i], tmpthree, modres);
                  sequence[nres] = blThronex(tmpthree);
               }
            }
               
            if(!gBioplibSeqNucleicAcid || doNucleic)
               nres++;
         }
      }
   }

   /* If no SEQRES records found, then free the memory and return NULL  */
   if(!strlen(sequence))
   {
      free(sequence);
      sequence = NULL;
   }
   else
   {
      if(AddStar)
      {
         sequence[nres++] = '*';
      }
      sequence[nres++] = '\0';
      if(chains!=NULL)
         chains[nchain][0] = '\0';
   }

   return(sequence);
}



/************************************************************************/
/*>MODRES *blGetModresWholePDB(WHOLEPDB *wpdb)
   -------------------------------------------
*//**
   \param[in]
   \return

   Reads MODRES records from a Whole PDB structure and returns a linked
   list containing the information

-  07.03.07  Original   By: ACRM
-  11.06.15  Moved to Bioplib
*/
MODRES *blGetModresWholePDB(WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   char *ch;
   MODRES *modres = NULL,
          *m = NULL;
   
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "MODRES", 6))
      {
         if(m==NULL)
         {
            INIT(modres, MODRES);
            m = modres;
         }
         else
         {
            ALLOCNEXT(m, MODRES);
         }
         if(m==NULL)
         {
            fprintf(stderr,"pdb2pir: Error! No memory for modres\n");
            exit(1);
         }
         
         ch = s->string+12;
         strncpy(m->modres, ch, 3);
         PADCHARMINTERM(m->modres, ' ', 4);
         
         ch = s->string+24;
         strncpy(m->origres, ch, 3);
         PADCHARMINTERM(m->origres, ' ', 4);
         if(m->origres[0] == ' ')
         {
            strncpy(m->origres, "XXX ", 4);
         }
      }
   }
   return(modres);
}


/************************************************************************/
/*>void blFindOriginalResType(char *modAA, char *stdAA, MODRES *modres)
   --------------------------------------------------------------------
*//**
   \param[in]  *modAA    Non-standard (modified) amino acid name
   \param[out] *stdAA    Standard amino acid from which it was derived
   \param[in]  *modres   MODRES linked list

   Uses the MODRES information to identify the original (standard)
   amino acid from which a modified amino acid was derived

-  07.03.07  Original   By: ACRM
-  11.06.15  Moved to bioplib, renamed routine and parameters
*/
void blFindOriginalResType(char *modAA, char *stdAA, MODRES *modres)
{
   MODRES *m;
   for(m=modres; m!=NULL; NEXT(m))
   {
      if(!strncmp(modAA, m->modres, 3))
      {
         strncpy(stdAA, m->origres, 3);
         PADCHARMINTERM(stdAA, ' ', 4);
         return;
      }
   }
}

/************************************************************************/
/*>static void checkRemark300(BIOMOLECULE *biomolecule, STRINGLIST *s, 
                              int *SkipStandardRemark)
   -------------------------------------------------------------------
*//**
   \param[in,out] *biomolecule        Pointer to a BIOMOLECULE structure
                                      that we populate
   \param[in]     *s                  The current header line
   \param[in,out] *SkipStandardRemark Counter of lines into the REMARK 300
                                      standard comments

   Checks header lines to see if they are REMARK 300 and, if so, extracts
   the maximum number of biomolecules and any 'details' comments that
   appear after the standard REMARK 300 comments

   biomolecule is a linked list that only has one item for REMARK 300.
   The 'details' that this code finds are just placed in the initial
   entry of the linked list.

-  26.06.15  Original   By: ACRM
*/
static void checkRemark300(BIOMOLECULE *biomolecule, STRINGLIST *s, 
                           int *SkipStandardRemark)
{
   if(!strncmp(s->string, "REMARK 300", 10))
   {
      if(!strncmp(s->string, "REMARK 300 BIOMOLECULE:", 23))
      {
         char buffer[80],
              *chp;
         
         *SkipStandardRemark = 1;
         
         /* Copy the actual data and remove trailing spaces             */
         strncpy(buffer, s->string+24, 80);
         TERMINATE(buffer);
         KILLTRAILSPACES(buffer);
         
         /* Now move to the last space                                  */
         if((chp = strrchr(buffer, ' '))==NULL)
         {
            chp = buffer;
         }
         else
         {
            chp++;
         }
         
         sscanf(chp, "%d", &(biomolecule->numBiomolecules));
      }
      
      if(*SkipStandardRemark)
      {
         if((*SkipStandardRemark)++ > 5)
         {
            char buffer[80];
            
            strncpy(buffer, s->string+11, 80);
            TERMINATE(buffer);
            KILLTRAILSPACES(buffer);
            if(strlen(buffer))
            {
               biomolecule->details = 
                  blStoreString(biomolecule->details, buffer);
            }
         }
      }
   }
}


/************************************************************************/
/*>static void checkRemark350(BIOMOLECULE *biomolecule, STRINGLIST *s)
   -------------------------------------------------------------------
*//**
   \param[in,out]  *biomolecule  Pointer to the BIOMOLECULE structure that
                                 we are populating
   \param[in]      *s            The current header record

   Checks a header line for REMARK 350. If more than one biomolecule is
   found, the biomolecule linked list is extended for each new one. Only
   the first one in the list will have the REMARK 300 details and the 
   number of biomolecules.

   This routine stores the author and software determined assembly size
   and the list of chains associated with a given biomolecule. It also
   contains a linked list of BIOMT structures which have the 
   transformation matrices needed to recreate the biomolecule.

-  26.06.15  Original   By: ACRM
*/
static void checkRemark350(BIOMOLECULE *biomolecule, STRINGLIST *s)
{
   BIOMOLECULE *bm         = NULL;
   BIOMT       *biomt      = NULL;
   static BOOL firstRecord = TRUE;
   
   bm=biomolecule;
   LAST(bm);

   if(!strncmp(s->string, "REMARK 350", 10))
   {
      if(!strncmp(s->string, "REMARK 350 BIOMOLECULE:", 23))
      {
         if(!firstRecord)
         {
            /* Allocate space for new biomolecule                          */
            ALLOCNEXT(bm, BIOMOLECULE);
            if(bm == NULL)
               return;
            CLEAR_BIOMOL(bm);
         }
         sscanf(s->string+23, "%d", &(bm->biomolNumber));
         
         firstRecord = FALSE;
      }
      else if(!strncmp(s->string, "REMARK 350 AUTHOR DETERMINED", 28))
      {
         strncpy(bm->authorUnit, s->string+46, 40);
         TERMINATE(bm->authorUnit);
         KILLTRAILSPACES(bm->authorUnit);
      }
      else if(!strncmp(s->string, "REMARK 350 SOFTWARE DETERMINED", 30))
      {
         strncpy(bm->softwareUnit, s->string+53, 40);
         TERMINATE(bm->softwareUnit);
         KILLTRAILSPACES(bm->softwareUnit);
      }
      else if(!strncmp(s->string, 
                       "REMARK 350 APPLY THE FOLLOWING TO CHAINS:", 41) ||
              !strncmp(s->string, 
                       "REMARK 350                    AND CHAINS:", 41))
      {
         char buffer[80],
              *chp;
         int  pos;
         
         for(chp=s->string+42, pos=0; 
             ((*chp!='\n') && (*chp!='\0'));
             chp++)
         {
            if(*chp != ' ')
            {
               buffer[pos++] = *chp;
            }
         }
         buffer[pos] = '\0';

         /* Remove any chain information already stored                 */
         if(bm->chains != NULL)
         {
            free(bm->chains);
            bm->chains = NULL;
         }
         
         bm->chains = blStrcatalloc(bm->chains, buffer);
      }
      else if(!strncmp(s->string, "REMARK 350   BIOMT", 18))
      {
         char buffer[80];
         int  line,
              entry;
         REAL val[4];
         
         strncpy(buffer, s->string+18, 80);
         TERMINATE(buffer);
         if(sscanf(buffer, "%d %d %lf %lf %lf %lf",
                   &line, &entry, &val[0], &val[1], &val[2], &val[3]))
         {
            /* Nothing defined yet so create entry and set entry number */
            if(bm->biomt == NULL)
            {
               INIT(bm->biomt, BIOMT);
               biomt = bm->biomt;
               if(biomt!=NULL)
                  biomt->biomtNum = entry;
            }
            else /* We already have something so go to end of list      */
            {
               biomt = bm->biomt;
               LAST(biomt);
            }

            /* If this entry is a different entry number, allocate a 
               new item
            */
            if(entry != biomt->biomtNum)
               ALLOCNEXT(biomt, BIOMT);
            
            /* Copy in the data                                         */
            if(biomt != NULL)
            {
               int i;
               
               biomt->biomtNum = entry;
               for(i=0; i<3; i++)
                  biomt->rotMatrix[line-1][i] = val[i];
               biomt->transMatrix[line-1] = val[3];
            }
         }
      }
   }
}



/************************************************************************/
/*>BIOMOLECULE *blGetBiomoleculeWholePDB(WHOLEPDB *wpdb)
   -----------------------------------------------------
*//**
   \param[in]   *wpdb   Pointer to WHOLEPDB structure
   \return              Pointer to malloc()'d BIOMOLECULE linked list

   Reads the biomolecule assembly data from REMARK 300 and REMARK 350
   headers.

   The returned pointer is a linked list of BIOMOLECULE structures.
   The first item in the list contains the information from REMARK 300:
   the total number of biomolecule assemblies and any additional
   information provided in REMARK 300 beyond the standard comments.
   The additional information (biomolecule->details) is a STRINGLIST
   linked list.

   The list of chains will appear in biomolcule->chains as a comma
   separated list and the transformation to be applied to them
   will appear in the biomolecule->biomt structure. This is also 
   a linked list so if multiple transformations are required these
   will appear in further entries to this list.

   Additional entries will be found in the BIOMOLECULE linked list
   for each different biomolecule.

   See the test code in PDBHeaderInfo.c for example usage.

TODO - return NULL if no data found!

-  26.06.15  Original   By: ACRM
*/
BIOMOLECULE *blGetBiomoleculeWholePDB(WHOLEPDB *wpdb)
{
   STRINGLIST  *s;
   BIOMOLECULE *biomolecule = NULL;
   int         SkipStandardRemark = 0;
   
   INIT(biomolecule, BIOMOLECULE);
   if(biomolecule == NULL)
      return(NULL);
   CLEAR_BIOMOL(biomolecule);

   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      checkRemark300(biomolecule, s, &SkipStandardRemark);
      checkRemark350(biomolecule, s); 
   }

   return(biomolecule);
}


/************************************************************************/
/*>void blFreeBiomolecule(BIOMOLECULE *biomolecule)
   ------------------------------------------------
   \param[in]    *biomolecule    Pointer to BIOMOLECULE linked list

   Frees the data associated with a BIOMOLECULE linked list

-  26.06.15 Original   By: ACRM
*/
void blFreeBiomolecule(BIOMOLECULE *biomolecule)
{
   BIOMOLECULE *bm;
   
   if(biomolecule == NULL)
      return;

   /* Free the REMARK 300 data stored in the first item in the list     */
   if(biomolecule->details != NULL) 
      FREESTRINGLIST(biomolecule->details); 
   
   /* Free chain and BIOMT data associated with each biomolecule        */
   for(bm=biomolecule; bm!=NULL; NEXT(bm))
   {
      if(bm->chains != NULL)
         free(bm->chains);

      if(bm->biomt)
         FREELIST((bm->biomt), BIOMT);
   }

   /* Free the biomolecule linked list                                  */
   FREELIST(biomolecule, BIOMOLECULE);
}




/************************************************************************/
#ifdef TEST
int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   FILE     *in;
   char     header[80],
            date[16],
            pdbcode[8],
            *title,
            **chainLabels;
   int      nChains,
            i;
   PDBSOURCE species;
   BIOMOLECULE *biomolecule = NULL,
               *bm = NULL;
   STRINGLIST  *s;
   
   if((in=fopen(argv[1], "r"))!=NULL)
   {
      if((wpdb = ReadWholePDB(in))!=NULL)
      {
         if(blGetHeaderWholePDB(wpdb, 
                                header, 80,
                                date,   16,
                                pdbcode, 8))
         {
            printf("Header:   '%s'\n", header);
            printf("Date:     '%s'\n", date);
            printf("PDB code: '%s'\n", pdbcode);
         }

         if((title = blGetTitleWholePDB(wpdb))!=NULL)
         {
            printf("Title:    '%s'\n", title);
         }

         if((biomolecule = blGetBiomoleculeWholePDB(wpdb))!=NULL)
         {
            printf("Number of Biomolecules: %d\n", biomolecule->numBiomolecules);
            for(s=biomolecule->details; s!=NULL; NEXT(s))
            {
               printf("REMARK 300 Details: %s\n", s->string);
            }

            for(bm=biomolecule; bm!=NULL; NEXT(bm))
            {
               BIOMT *bmt;
               
               printf("Biomolecule: %d\n", bm->biomolNumber);
               printf("   Author Unit:   %s\n", bm->authorUnit);
               printf("   Software Unit: %s\n", bm->softwareUnit);
               printf("   Chains:        %s\n", (bm->chains?bm->chains:""));

               for(bmt=bm->biomt; bmt!=NULL; NEXT(bmt))
               {
                  int i;
                  
                  printf("   Matrix %d\n", bmt->biomtNum);
                  for(i=0; i<3; i++)
                  {
                     printf("      %8.6f %8.6f %8.6f   %8.6f\n",
                            bmt->rotMatrix[i][0],
                            bmt->rotMatrix[i][1],
                            bmt->rotMatrix[i][2],
                            bmt->transMatrix[i]);
                  }
               }
            }

            blFreeBiomolecule(biomolecule);
         }

         chainLabels = blGetPDBChainLabels(wpdb->pdb, &nChains);
         for(i=0; i<nChains; i++)
         {
            COMPND compound;

            printf("\n\n>>>Chain: %s\n", chainLabels[i]);

            blGetCompoundWholePDBChain(wpdb, chainLabels[i], &compound);

            printf("molid:      %d\n", compound.molid);
            printf("molecule:   %s\n", compound.molecule);
            printf("chain:      %s\n", compound.chain);
            printf("fragment:   %s\n", compound.fragment);
            printf("synonym:    %s\n", compound.synonym);
            printf("ec:         %s\n", compound.ec);
            printf("engineered: %s\n", compound.engineered);
            printf("mutation:   %s\n", compound.mutation);
            printf("other:      %s\n", compound.other);

            if(blGetSpeciesWholePDBChain(wpdb, chainLabels[i], &species))
            {
               printf("Scientific name: %s\n", species.scientificName);
               printf("Common name:     %s\n", species.commonName);
               printf("Strain:          %s\n", species.strain);
               printf("Tax ID:          %d\n", species.taxid);
            }
            
            free(chainLabels[i]);
         }
         free(chainLabels);
      }
   }
   
   return(0);
}
#endif


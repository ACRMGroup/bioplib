/************************************************************************/
/**

   \file       PDBHeaderInfo.c
   
   \version    V1.0
   \date       26.03.15
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

   #FUNCTION blGetSpeciesWholePDBChain()
   Obtains the species data for a specified chain from WHOLEPDB info

*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "pdb.h"
#include "macros.h"

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
char *collapseSpaces(char *inText);


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
         char buffer[MAXBUFF];
         strcpy(buffer, s->string);
         TERMINATE(buffer);
         
         if(buffer[9] == ' ')
            title = blStrcatalloc(title, buffer+10);
         else
            title = blStrcatalloc(title, buffer+11);

         if(title == NULL)
            return(NULL);
      }
      else if(inTitle)
      {
         break;
      }
   }

   cleanTitle = collapseSpaces(title);
   free(title);
   KILLTRAILSPACES(cleanTitle);
   
   return(cleanTitle);
}

/************************************************************************/
char *collapseSpaces(char *inText)
{
   int  nchar = 0;
   char *ch, *chp, *chq,
        *outText = NULL;

   if(inText==NULL)
      return(NULL);
   
   /* Count characters skipping repeated spaces                         */
   chp=NULL;
   for(ch=inText; *ch!='\0'; ch++)
   {
      if(((chp==NULL) || (*chp!=' ')) || (*ch!=' '))
      {
         nchar++;
      }
      chp=ch;
   }
   /* Increment count for terminator                                    */
   nchar++;

   /* Allocate new space                                                */
   if((outText=(char *)malloc(nchar * sizeof(char)))==NULL)
      return(NULL);

   /* Copy characters skipping repeated spaces                          */
   chp=NULL;
   chq=outText;
   for(ch=inText; *ch!='\0'; ch++)
   {
      if(((chp==NULL) || (*chp!=' ')) || (*ch!=' '))
      {
         *chq = *ch;
         chq++;
      }
      chp=ch;
   }
   *chq = '\0';

   return(outText);
}


typedef struct _compnd
{
   char  molecule[MAXBUFF],
         chain[MAXBUFF],
         fragment[MAXBUFF],
         synonym[MAXBUFF],
         ec[MAXBUFF],
         engineered[MAXBUFF],
         mutation[MAXBUFF],
         other[MAXBUFF];
}  COMPND;

typedef struct _pdbsource
{
   char scientificName[160],
        commonName[160],
        strain[160];
   int  taxid;
}  PDBSOURCE;

   
/************************************************************************/
STRINGLIST *FindNextMolIDRecord(STRINGLIST *start)
{
   STRINGLIST *s;
   
   if(start==NULL)
      return(NULL);
   
   for(s=start->next; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "COMPND", 6))
      {
         if(strstr(s->string, "MOL_ID:"))
            return(s);
      }
   }
   return(NULL);
}

#define MAXWORD 8


/************************************************************************/
BOOL PopulateCompndField(STRINGLIST *molidStart, STRINGLIST *molidStop,
                         char *data, char *field)

{
   STRINGLIST *s;
   BOOL       GotField = FALSE;
   char       *chp,
              buffer[MAXBUFF];

   data[0] = '\0';

   for(s=molidStart; s!=molidStop; NEXT(s))
   {
      if(strncmp(s->string, "COMPND", 6))
         break;

      chp = NULL;

      if(GotField)
      {
         /* We have found the field already on previous line            */
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
         strncpy(buffer, chp, MAXBUFF);
         /* Remove spaces                                               */
         TERMINATE(buffer);
         KILLTRAILSPACES(buffer);
         /* Add to output data                                          */
         blStrncat(data, buffer, MAXBUFF);

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
BOOL blGetCompoundWholePDBChain(WHOLEPDB *wpdb, char *chain, 
                                COMPND *compnd)
{
   STRINGLIST *molidFirst,
              *molidStart,
              *molidStop,
              *s;
   BOOL       foundChain = FALSE;


   compnd->molecule[0]   = '\0';
   compnd->chain[0]      = '\0';
   compnd->fragment[0]   = '\0';
   compnd->synonym[0]    = '\0';
   compnd->ec[0]         = '\0';
   compnd->engineered[0] = '\0';
   compnd->mutation[0]   = '\0';
   compnd->other[0]      = '\0';

   molidFirst = FindNextMolIDRecord(wpdb->header);

   for(molidStart=molidFirst; molidStart!=NULL; molidStart=molidStop)
   {
      molidStop  = FindNextMolIDRecord(molidStart);
      for(s=molidStart; s!=molidStop; NEXT(s))
      {
         char *chp;
         char buffer[MAXBUFF],
              word[MAXWORD];
         
         if((chp=strstr(s->string, "CHAIN:"))!=NULL)
         {
            strncpy(buffer, chp+6, MAXBUFF);
            TERMAT(buffer, ';');
            KILLLEADSPACES(chp,buffer);

            /* Check the chains to see if our chain is there            */
            do {
               chp=blGetWord(chp, word, MAXWORD);
               if(!strcmp(word, chain))
               {
                  foundChain = TRUE;
                  goto found;
               }
            }  while(chp!=NULL);
         }
      }
   }
   
found:
   if(foundChain)
   {
      PopulateCompndField(molidStart, molidStop,
                          compnd->molecule,   "MOLECULE:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->chain,      "CHAIN:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->fragment,   "FRAGMENT:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->synonym,    "SYNONYM:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->ec,         "EC:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->engineered, "ENGINEERED:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->mutation,   "MUTATION:");
      PopulateCompndField(molidStart, molidStop,
                          compnd->other,      "OTHER:");
      return(TRUE);
   }
   
   return(FALSE);
}

/************************************************************************/
BOOL blGetSpeciesWholePDBChain(WHOLEPDB *wpdb, char *chain,
                               PDBSOURCE *source)
{
   source->scientificName[0] = '\0';
   source->commonName[0]     = '\0';
   source->strain[0]         = '\0';
   source->taxid             = 0;
   
   return(FALSE);
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
   
   if((in=fopen("test.pdb", "r"))!=NULL)
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

         chainLabels = blGetPDBChainLabels(wpdb->pdb, &nChains);
         for(i=0; i<nChains; i++)
         {
            COMPND compound;

            printf("\n\n>>>Chain: %s\n", chainLabels[i]);

            blGetCompoundWholePDBChain(wpdb, chainLabels[i], &compound);

            printf("molecule: %s\n", compound.molecule);
            printf("chain: %s\n", compound.chain);
            printf("fragment: %s\n", compound.fragment);
            printf("synonym: %s\n", compound.synonym);
            printf("ec: %s\n", compound.ec);
            printf("engineered: %s\n", compound.engineered);
            printf("mutation: %s\n", compound.mutation);
            printf("other: %s\n", compound.other);

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


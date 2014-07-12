/************************************************************************/
/**

   \file       WholePDB.c
   
   \version    V1.7
   \date       07.07.14
   \brief      
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2014
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
-  V1.0  30.05.02 Original
-  V1.1  12.06.08 CTP Added include for port.h
-  V1.2  13.06.08 popen() and pclose() prototypes skipped for Mac OS X.
-  V1.3  17.03.09 popen() prototype skipped for Windows. By: CTP
-  V1.4  22.04.14 Added handling of PDBML files. By CTP
-  V1.5  21.06.14 Updated writing of PDBML files. By: CTP
-  V1.6  06.07.14 Defined _XOPEN_SOURCE and __USE_XOPEN. Required for 
                  time.h on some linux systems. By: CTP
-  V1.7  07.07.14 Use renamed functions with bl prefix. By: CTP

*************************************************************************/
/* Includes
*/
#include "port.h"    /* Required before stdio.h                         */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "macros.h"
#include "general.h"
#include "pdb.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#define _XOPEN_SOURCE
#define __USE_XOPEN
#include <time.h>


/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF      160
#define XML_BUFFER  1024

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static WHOLEPDB *bldoReadWholePDB(FILE *fpin, BOOL atomsonly);
static STRINGLIST *blParseHeaderPDBML(FILE *fpin);

#if !defined(__APPLE__) && !defined(MS_WINDOWS)
FILE *popen(char *, char *);
#endif
#ifndef __APPLE__
int  pclose(FILE *);
#endif

/************************************************************************/
/*>void blFreeWholePDB(WHOLEPDB *wpdb)
   ---------------------------------
*//**

   \param[in]     *wpdb    WHOLEPDB structure to be freed

   Frees the header, trailer and atom content from a WHOLEPDB structure

-  30.05.02  Original   By: ACRM
-  07.07.14  Renamed to blFreeWholePDB() By: CTP
*/
void blFreeWholePDB(WHOLEPDB *wpdb)
{
   FreeStringList(wpdb->header);
   FreeStringList(wpdb->trailer);
   FREELIST(wpdb->pdb, PDB);
   free(wpdb);
}

/************************************************************************/
/*>BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
   ----------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes a PDB file including header and trailer information.
   Output in PDBML-format if flags set.

-  21.06.14  Original   By: CTP
*/
BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
{
   if((gPDBXMLForce == FORCEXML_XML) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == TRUE))
   {
      /* Write PDBML file (omitting header and footer data) */
      blWriteAsPDBML(fp, wpdb->pdb);
   }
   else
   {
      /* Check format */
      if(blFormatCheckWritePDB(wpdb->pdb) == FALSE)
      {
         return FALSE;
      }

      /* Write whole PDB File */
      blWriteWholePDBHeader(fp, wpdb);
      blWriteAsPDB(fp, wpdb->pdb);
      blWriteWholePDBTrailer(fp, wpdb);
   }
   
   return TRUE;
}


/************************************************************************/
/*>void blWriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
   ----------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes the header of a PDB file 

-  30.05.02  Original   By: ACRM
-  21.06.14  Renamed to blWriteWholePDBHeader() By: CTP
*/
void blWriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      fputs(s->string, fp);
   }
}


/************************************************************************/
/*>void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
   -----------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes the trailer of a PDB file 

-  30.05.02  Original   By: ACRM
-  21.06.14  Renamed to blWriteWholePDBTrailer() By: CTP
*/
void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   
   for(s=wpdb->trailer; s!=NULL; NEXT(s))
   {
      fputs(s->string, fp);
   }
}


/************************************************************************/
/*>WHOLEPDB *blReadWholePDB(FILE *fpin)
   ------------------------------------
*//**

   \param[in]     *fpin     File pointer
   \return                  Whole PDB structure containing linked
                            list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

-  07.03.07 Made into a wrapper to doReadWholePDB()
-  07.07.14 Use bldoReadWholePDB() Renamed to blReadWholePDB() By: CTP
*/
WHOLEPDB *blReadWholePDB(FILE *fpin)
{
   return(bldoReadWholePDB(fpin, FALSE));
}

/************************************************************************/
/*>WHOLEPDB *blReadWholePDBAtoms(FILE *fpin)
   -----------------------------------------
*//**

   \param[in]     *fpin     File pointer
   \return                  Whole PDB structure containing linked
                            list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

-  07.03.07 Made into a wrapper to doReadWholePDB()
-  07.07.14 Use bldoReadWholePDB() Renamed to blReadWholePDBAtoms() 
            By: CTP
*/
WHOLEPDB *blReadWholePDBAtoms(FILE *fpin)
{
   return(bldoReadWholePDB(fpin, TRUE));
}


/************************************************************************/
/*>static WHOLEPDB *bldoReadWholePDB(FILE *fpin, BOOL atomsonly)
   -------------------------------------------------------------
*//**

   \param[in]     *fpin       File pointer
   \param[in]     atomsonly   TRUE:  Read ATOM records only
                              FALSE: Read ATOM & HETATM records
   \return                    Whole PDB structure containing linked
                              list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

-  30.05.02 Original   By: ACRM
-  07.03.07 Made into a doXXX routine to add a atomsonly parameter
-  05.06.07 Added support for Unix compress'd files
-  22.04.14 Handles PDBML format. By: CTP
-  07.07.14 Use Renamed ReadPDB functions. Use blParseHeaderPDBML().
            Renamed to bldoReadWholePDB() By: CTP

   TODO FIXME!!!!! Move all this into doReadPDB so that we don't worry 
   about rewinding any more
*/
static WHOLEPDB *bldoReadWholePDB(FILE *fpin, BOOL atomsonly)
{
   WHOLEPDB *wpdb;
   char     buffer[MAXBUFF];
   FILE     *fp = fpin;
   BOOL     pdbml_format = FALSE;
   
#ifdef GUNZIP_SUPPORT
   int      signature[3],
            i,
            ch;
   char     cmd[80];
#endif

   if((wpdb=(WHOLEPDB *)malloc(sizeof(WHOLEPDB)))==NULL)
      return(NULL);

   wpdb->pdb     = NULL;
   wpdb->header  = NULL;
   wpdb->trailer = NULL;
   
#ifdef GUNZIP_SUPPORT
   cmd[0] = '\0';
   
   /* See whether this is a gzipped file                                */
   for(i=0; i<3; i++)
      signature[i] = fgetc(fpin);
   for(i=2; i>=0; i--)
      ungetc(signature[i], fpin);
   if(((signature[0] == (int)0x1F) &&    /* gzip                        */
       (signature[1] == (int)0x8B) &&
       (signature[2] == (int)0x08)) ||
      ((signature[0] == (int)0x1F) &&    /* 05.06.07 compress           */
       (signature[1] == (int)0x9D) &&
       (signature[2] == (int)0x90)))
   {
      /* It is gzipped so we'll open gunzip as a pipe and send the data
         through that into a temporary file
      */
      cmd[0] = '\0';
      sprintf(cmd,"gunzip >/tmp/readpdb_%d",(int)getpid());
      if((fp = (FILE *)popen(cmd,"w"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
      while((ch=fgetc(fpin))!=EOF)
         fputc(ch, fp);
      pclose(fp);

      /* We now reopen the temporary file as our PDB input file         */
      sprintf(cmd,"/tmp/readpdb_%d",(int)getpid());
      if((fp = fopen(cmd,"r"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
   }
#endif   


   /* Check file format */
   pdbml_format = blCheckFileFormatPDBML(fp);

   /* Read the header from the PDB file                                 */
   if(!pdbml_format)
   {
      while(fgets(buffer,MAXBUFF,fp))
      {
         if(!strncmp(buffer, "ATOM  ", 6) ||
            !strncmp(buffer, "HETATM", 6) ||
            !strncmp(buffer, "MODEL ", 6))
      {
         break;
      }
      if((wpdb->header = StoreString(wpdb->header, buffer))==NULL)
         return(NULL);
      }
   }
   else
   {
      wpdb->header = blParseHeaderPDBML(fp);
   }
   
   /* Read the coordinates                                              */
   rewind(fp);
   if(atomsonly)
   {
      wpdb->pdb = blReadPDBAtoms(fp, &(wpdb->natoms));
   }
   else
   {
      wpdb->pdb = blReadPDB(fp, &(wpdb->natoms));
   }

   /* Read the trailer                                                  */
   rewind(fp);
   if(!pdbml_format)
   {
      while(fgets(buffer,MAXBUFF,fp))
      {
         if(!strncmp(buffer, "CONECT", 6) ||
            !strncmp(buffer, "MASTER", 6) ||
            !strncmp(buffer, "END   ", 6))
         {
            wpdb->trailer = StoreString(wpdb->trailer, buffer);
         }
      }
   }
   else
   {
      wpdb->trailer = StoreString(wpdb->trailer, "END   \n");
   }
   
   return(wpdb);
}

/************************************************************************/
/*>static STRINGLIST *blParseHeaderPDBML(FILE *fpin)
   -------------------------------------------------
*//**

   \param[in]     *fpin   File pointer
   \return                STRINGLIST with basic header information.

   Parses a PDBML file and creates HEADER and TITLE lines.

-  22.04.14 Original. By: CTP
-  07.07.14 Renamed to blParseHeaderPDBML() By: CTP

*/
static STRINGLIST *blParseHeaderPDBML(FILE *fpin)
{
   xmlParserCtxtPtr ctxt;
   xmlDoc  *document;
   xmlNode *root_node = NULL, 
           *node      = NULL,
           *subnode   = NULL,
           *n         = NULL;
   int     size_t;
   char    xml_buffer[XML_BUFFER];
   xmlChar *content, *attribute;
   
   STRINGLIST *wpdb_header = NULL,
              *title_lines = NULL;

   char header_line[82]  = "",
        title_line[82]   = "",
        header_field[41] = "",
        pdb_field[5]     = "",
        date_field[10]   = "",
        title_field[71]  = "";
        
   int cut_from = 0, 
       cut_to   = 0,
       nlines   = 0,
       i        = 0;

   struct tm date = {0};

   
   /* Generate Document From Filehandle                                 */
   size_t = fread(xml_buffer, 1, XML_BUFFER, fpin);
   ctxt = xmlCreatePushParserCtxt(NULL, NULL, xml_buffer, size_t, "file");
   while ((size_t = fread(xml_buffer, 1, XML_BUFFER, fpin)) > 0) 
   {
      xmlParseChunk(ctxt, xml_buffer, size_t, 0);
   }
   xmlParseChunk(ctxt, xml_buffer, 0, 1);
   document = ctxt->myDoc;
   xmlFreeParserCtxt(ctxt);

   if(document == NULL){ return(NULL); } /*        failed to parse file */


   /* Parse Document Tree */
   root_node = xmlDocGetRootElement(document);
   for(node = root_node->children; node; node = node->next)
   {
      if(node->type != XML_ELEMENT_NODE){ continue; }

      /* get header                                                     */
      if(!strcmp("struct_keywordsCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("pdbx_keywords",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);
               strncpy(header_field,(char *) content,40);
               xmlFree(content);
            }
         }
      }

      /* get date                                                       */
      if(!strcmp("database_PDB_revCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("date_original",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);
               
               /* convert date format */
               strptime((char *) content,"%Y-%m-%d",&date);
               strftime(date_field, sizeof(date_field), "%d-%b-%y", &date);
               UPPER(date_field);
               
               xmlFree(content);
            }
         }
      }

      /* get pdb code                                                   */
      if(!strcmp("entryCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            if(strcmp("entry",(char *) subnode->name)){ continue; }
            attribute = xmlGetProp(subnode,(xmlChar *) "id");
            strncpy(pdb_field,(char *)attribute,4);
            pdb_field[4] = '\0';
            xmlFree(attribute);
         }
      }

      /* get title                                                      */
      if(!strcmp("structCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("title",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);

               /* Get title lines as STRINGLIST                         */
               for(i=0; i<strlen((char *) content); i++)
               {
                  if(content[i] == ' ' ) cut_to = i;
                  if(i == strlen((char *) content) - 1) cut_to = i+1;

                  /* split and store title line                         */
                  if( (i && !((i - cut_from)%70)) || 
                      i == strlen((char *) content)-1 )
                  {
                     nlines++;
                     cut_to = (cut_from == cut_to) ? i : cut_to;
                     strncpy(title_field,
                             (char *) &content[cut_from],
                             cut_to - cut_from);
                     title_field[cut_to - cut_from] = '\0';
                     PADMINTERM(title_field,70);
                     cut_from = cut_to;
                     i        = cut_to;

                     if(nlines == 1)
                     {
                        sprintf(title_line, "TITLE     %s\n",
                                title_field);
                        title_lines = StoreString(NULL,title_line);
                     }
                     else
                     {
                        sprintf(title_line, "TITLE   %2d%s\n", nlines,
                                title_field);
                        StoreString(title_lines,title_line);
                     }
                  }
               }
               xmlFree(content);
            }
         }
      }
   }

   /* Free document */
   xmlFreeDoc(document);

   /* Cleanup xml parser                                                */
   xmlCleanupParser();

   /* Create Header Line                                                */
   if(!strlen(header_field))
   {
      strcpy(header_field,"Converted from PDBML");
   }
   sprintf(header_line, "HEADER    %-40s%9s   %4s\n",
           header_field, date_field, pdb_field);
   
   /* Make Stringlist                                                   */
   wpdb_header = StoreString(wpdb_header, header_line);
   wpdb_header->next = title_lines;
   
   return(wpdb_header);
}

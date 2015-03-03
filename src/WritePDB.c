/************************************************************************/
/**

   \file       WritePDB.c
   
   \version    V1.21
   \date       02.03.15
   \brief      Write a PDB file from a linked list
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2015
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

   This routine will write a .PDB file of any size from a linked list of 
   the protein structure. This list is contained in a linked set of 
   structures of type pdb_entry. The structure is set up by including the 
   file "pdb.h". For details of the structure, see this file.

**************************************************************************

   Usage:
   ======

\code
   blWritePDB(fp, pdb)
\endcode
   
   \param[in]     *fp      A pointer to the file to write
   \param[in]     *pdb     The start of the PDB linked list.

**************************************************************************

   Revision History:
   =================
-  V1.0  08.03.89 Original
-  V1.2  28.03.90 Modified to match the correct column definition of 
                  ReadPDB V1.2   (N.B. There was no V1.1)
-  V1.3  01.06.92 Corrected header, to match standard. Autodoc'd, 
                  ANSIed. Added FPU check.
-  V1.4  10.06.93 Changed to use NEXT() macro. void types
-  V1.5  22.02.94 Added TER card at end of file
-  V1.6  15.02.01 Writes using atnam_raw so atom name is unchanged from
                  input
-  V1.7  30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.8  03.06.05 'atnam_raw' no longer includes the alternate indicator
                  which is now in 'altpos'
-  V1.9  22.09.06 Added WritePDBRecordAtnam()
-  V1.10 04.02.14 Use CHAINMATCH macro. By: CTP
-  V1.11 01.06.14 Added WritePDBML() By: CTP
-  V1.12 21.06.14 Added blWritePDB(), blFormatCheckWritePDB() and 
                  blWriteAsPDB(). Renamed WritePDBML() to blWriteAsPDBML()
                  and deprecated WritePDB(). Defined WRITEPDB_MAIN.
                  By: CTP
-  V1.13 07.07.14 Renamed functions to use bl prefix. Moved WritePDB() to 
                  deprecated.h By: CTP
-  V1.14 17.07.14 Added blSetElementSymbolFromAtomName() By: CTP
-  V1.15 16.08.14 Added writing element and charge. By: CTP
-  V1.16 18.08.14 Added XML_SUPPORT option allowing compilation without 
                  support for PDBML format. By: CTP
-  V1.17 17.02.15 Handles segid and the element and formal charge in 
                  blWritePDBRecordAtnam()   By: ACRM
-  V1.18 23.02.15 Modified blWriteAsPDB() to do proper TER cards - now
                  returns an int - the number of TER cards written
-  V1.19 24.02.15 Renamed blWriteAsPDB() to blWritePDBAsPDBorGromos()
                  and integrated Gromos support into this routine.
                  Now takes a new BOOL flag
-  V1.20 25.02.15 blWritePDBAsPDBML() now returns BOOL and checks all
                  memory allocations
-  V1.21 02.03.15 Corrected counting of ORIGX, SCALE and MTRIX records
                  in WriteMaster().
                  Added space padding to TER in blWriteTerCard()
                  Added space padding to END and to CONECT in
                  blWriteWholePDBTrailer()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP File IO

   #KEYFUNCTION blWritePDB()
   Main entry point to write a PDB linked list to a file

   #KEYFUNCTION blWritePDBRecord()
   Writes a single PDB record in PDB format

   #KEYFUNCTION  blWriteWholePDB()
   Writes a PDB file including header and trailer information.
   Output in PDBML-format if flags set.

   #FUNCTION  blWriteWholePDBHeader()
   Writes the header of a PDB file 

   #FUNCTION blWriteWholePDBHeaderNoRes()
   Writes the header of a PDB file, but skips any records associated
   with residue numbers

   #FUNCTION  blWriteWholePDBTrailer()
   Writes the trailer of a PDB file 

   #FUNCTION blFormatCheckWritePDB()
   Checks that a PDB linked list can be written as a standard PDB file
   (i.e. chain labels are no more than one character)

   #FUNCTION blWritePDBAsPDBorGROMOS()
   Writes a PDB linked list to a file in PDB format or the modified
   GROMOS version

   #FUNCTION blWritePDBRecordAtnam()
   Writes a single PDB record in PDB format using atom data from the
   atnam field rather than the atnam_raw field

   #FUNCTION blWritePDBAsPDBML()
   Write a PDB linked list to a file in PDBML XML format

   #FUNCTION blSetElementSymbolFromAtomName()
   Sets the element field based on the content of the atom name stored 
   in atnam_raw

   #FUNCTION  blWriteGromosPDB()
   Write a PDB linked list by calling blWritePDBAsPDBorGromos()

   #FUNCTION  blWriteGromosPDBRecord()
   Write a GROMOS PDB record

*/
/************************************************************************/
/* Defines required for includes
*/
#define WRITEPDB_MAIN

/************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef XML_SUPPORT /* Required to read PDBML files                      */
#include <libxml/tree.h>
#include <ctype.h>
#define XMLDIE(x) do {if((x)!=NULL) { xmlFreeDoc((x));                   \
                                      xmlCleanupParser(); }              \
                      return(FALSE);} while(FALSE)
#endif

#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/* Prototypes
*/
static void WriteMaster(FILE *fp, WHOLEPDB *wpdb, int numConect,
                        int numTer);

/************************************************************************/
/*>int blWritePDB(FILE *fp, PDB *pdb)
   ----------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write
   \return        int   Number of TER cards written (0 indicates error)

   Write a PDB linked list...

-  21.06.14 Original By: CTP
-  18.08.14 Added XML_SUPPORT option. Return error if attempting to write 
            PDBML format. By: CTP
-  23.02.15 Now returns an int
-  24.02.15 Changed to call blWritePDBAsPDBorGromos() and blWritePDBAsPDBML()
*/
int blWritePDB(FILE *fp,
               PDB  *pdb)
{
   int numTer;

   if((gPDBXMLForce == FORCEXML_XML) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == TRUE))
   {

#ifdef XML_SUPPORT
      /* Write PDBML file                                               */
      blWritePDBAsPDBML(fp, pdb);
      return(1);
#else
      /* PDBML not supported                                            */
      return FALSE;
#endif

   }
   else
   {
      /* Check format                                                   */
      if(blFormatCheckWritePDB(pdb) == FALSE)
      {
         return(0);
      }

      /* Write whole PDB File                                           */
      numTer = blWritePDBAsPDBorGromos(fp, pdb, FALSE);
   }
   
   return(numTer);
}


/************************************************************************/
/*>BOOL blFormatCheckWritePDB(PDB *pdb)
   ------------------------------------
*//**

   \param[in]     *pdb  PDB linked list to write

   Checks PDB linked list is compatible with PDB-formatted text file.

-  21.06.14 Original By: CTP
*/
BOOL blFormatCheckWritePDB(PDB *pdb)
{
   PDB *p;
   for(p = pdb ; p ; NEXT(p))
   {
      /* Check chain id is single letter */
      if(strlen(p->chain) > 1)
      {
         return FALSE;
      }
   }
   return TRUE;
}


/************************************************************************/
/*>int blWritePDBAsPDBorGromos(FILE *fp, PDB *pdb, BOOL doGromos)
   --------------------------------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write
   \return              Number of TER cards written (0=error)

   Write a PDB linked list by calls to WritePDBRecord()

-  08.03.89 Original
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 Uses NEXT macro; void type
-  08.07.93 Added insertion of TER cards
-  22.02.94 And a TER card at the end of the file
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  17.06.14 Renamed to blWriteAsPDB() By: CTP
-  07.07.14 Use blWritePDBRecord() By: CTP
-  23.02.15 Write correct format TER cards. Now returns int    By: ACRM
-  24.02.15 Renamed to blWritePDBAsPDBorGromos() and added doGromos
            flag

*/
int blWritePDBAsPDBorGromos(FILE *fp, PDB  *pdb, BOOL doGromos)
{
   PDB   *p,
         *prev = NULL;
   int   numTer = 0;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If previous was non-null and was an ATOM                       */
      if((prev!=NULL) && !strncmp(prev->record_type, "ATOM  ", 6))
      {
         /* If the chain has changed or this is the start of HETATMs    */
         if(!CHAINMATCH(p->chain, prev->chain) || 
            !strncmp(p->record_type, "HETATM", 6))
         {
            blWriteTerCard(fp, prev);
            numTer++;
         }
      }
      if(doGromos)
      {
         blWriteGromosPDBRecord(fp,p);
      }
      else
      {
         blWritePDBRecord(fp,p);
      }
      prev=p;
   }
   if(!strncmp(prev->record_type, "ATOM  ", 6))
   {
      blWriteTerCard(fp, prev);
      numTer++;
   }
   return(numTer);
}


/************************************************************************/
/*>void blWriteTerCard(FILE *fp, PDB *p)
   -------------------------------------
*//**
   \param[in]    *fp    File pointer
   \param[in]    *p     PDB record pointer

   Prints a TER card in the new PDB format - i.e. with the residue 
   information for the previous ATOM/HETATM, rather than just printing
   TER

-  23.02.15  Original   By: ACRM
-  02.03.15  Added space padding
*/
void blWriteTerCard(FILE *fp, PDB *p)
{
   if(p!=NULL)
   {
      fprintf(fp,"TER   %5d      %-4s%1s%4d%1s%s\n",
              p->atnum+1, p->resnam, p->chain, p->resnum, p->insert,
              "                                                     ");
   }
}


/************************************************************************/
/*>void blWritePDBRecord(FILE *fp, PDB *pdb)
   -----------------------------------------
*//**

   \param[in]     *fp     PDB file pointer to be written
   \param[in]     *pdb    PDB linked list record to write

   Write a PDB record

-  08.03.89 Original
-  28.03.90 Changed to match ReadPDB() V1.2 for column widths
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 void type
-  22.06.93 Changed to %lf. Ljust strings
-  11.03.94 %lf back to %f (!)
-  15.02.01 Modified to use atnam_raw
-  03.06.05 Modified to use altpos
-  07.07.14 Renamed to blWritePDBRecord() By: CTP
-  16.08.14 Write element and formal charge.  By: CTP
-  17.02.15 Added segid support   By: ACRM
*/
void blWritePDBRecord(FILE *fp,
                      PDB  *pdb)
{
   char charge = ' ',
        sign   = ' ';

   if(pdb->formal_charge && ABS(pdb->formal_charge <= 8))
   {
      charge = (char)('0' + ABS(pdb->formal_charge));
      sign   = pdb->formal_charge > 0 ? '+':'-';
   }

   fprintf(fp,"%-6s%5d %-4s%c%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%c%c\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam_raw,
           pdb->altpos,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval,
           pdb->segid,
           pdb->element,
           charge,
           sign);
}


/************************************************************************/
/*>void blWritePDBRecordAtnam(FILE *fp, PDB *pdb)
   ----------------------------------------------
*//**

   \param[in]     *fp     PDB file pointer to be written
   \param[in]     *pdb    PDB linked list record to write

   Write a PDB record using the data in atnam rather than atnam_raw

-  08.03.89 Original
-  28.03.90 Changed to match ReadPDB() V1.2 for column widths
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 void type
-  22.06.93 Changed to %lf. Ljust strings
-  11.03.94 %lf back to %f (!)
-  15.02.01 Modified to use atnam_raw
-  03.06.05 Modified to use altpos
-  22.09.05 This is like the old version which used atnam rather
            than atnam_raw
-  07.07.14 Renamed to blWritePDBRecordAtnam() By: CTP
-  17.02.15 Added element, formalcharge and segid support   By: ACRM
*/
void blWritePDBRecordAtnam(FILE *fp,
                           PDB  *pdb)
{
   char charge = ' ',
        sign   = ' ';

   if(pdb->formal_charge && ABS(pdb->formal_charge <= 8))
   {
      charge = (char)('0' + ABS(pdb->formal_charge));
      sign   = pdb->formal_charge > 0 ? '+':'-';
   }

   fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%c%c\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval,
           pdb->segid,
           pdb->element,
           charge,
           sign);
}


/************************************************************************/
/*>BOOL blWritePDBAsPDBML(FILE *fp, PDB *pdb)
   ------------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write
   \return              Success

   Write a PDB linked list in PDBML format.

-  02.06.14 Original. By: CTP
-  21.06.14 Renamed blWriteAsPDBML() and updated symbol handling. By: CTP
-  17.07.14 Use blSetElementSymbolFromAtomName() By: CTP
-  16.08.14 Use element and charge data. By: CTP
-  17.02.15 Added segid support  By: ACRM
-  24.02.15 Changed name to blWritePDBAsPDBML()
-  25.02.15 Changed to type BOOL and checks all memory allocations

*/
BOOL blWritePDBAsPDBML(FILE *fp, PDB  *pdb)
{
#ifndef XML_SUPPORT

   /* PDBML format not supported.                                       */
   return;

#else 

   /* PDBML format supported                                            */
   PDB         *p;
   xmlDocPtr   doc         = NULL;
   xmlNodePtr  root_node   = NULL, 
               sites_node  = NULL, 
               atom_node   = NULL, 
               node        = NULL;
   xmlNsPtr    pdbx        = NULL,
               xsi         = NULL;
   char        buffer[16], 
               *buffer_ptr;
   
   /* Create document                                                   */
   if((doc = xmlNewDoc((xmlChar *)"1.0"))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */
   if((doc->encoding = xmlStrdup((xmlChar *) "UTF-8"))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */
   
   /* Root node                                                         */
   if((root_node=xmlNewNode(NULL, (xmlChar *)"datablock"))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */
   xmlDocSetRootElement(doc, root_node);
   if((pdbx=xmlNewNs(root_node, (xmlChar *)"null", 
                     (xmlChar *)"PDBx"))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */
   if((xsi =xmlNewNs(root_node, (xmlChar *)"null",
                     (xmlChar *)"xsi"))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */
   xmlSetNs(root_node,pdbx);
   
   
   /* Atom_sites node                                                   */
   if((sites_node = xmlNewChild(root_node, NULL, 
                                (xmlChar *)"atom_siteCategory", 
                                NULL))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */

   
   /* Atom nodes                                                        */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* skip TER                                                       */
      if(!strncmp("TER",p->resnam,3))
      {
         continue;
      }

      /* Add atom node                                                  */
      if((atom_node = xmlNewChild(sites_node, NULL,
                                  (xmlChar *)"atom_site", NULL))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      sprintf(buffer, "%d", p->atnum);
      xmlNewProp(atom_node, (xmlChar *)"id", (xmlChar *)buffer);
      
      /*** Add atom data nodes                                        ***/

      /* B value                                                        */
      sprintf(buffer,"%.2f", p->bval);
      if((node = xmlNewChild(atom_node, NULL, 
                             (xmlChar *)"B_iso_or_equiv",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      /* coordinates                                                    */
      sprintf(buffer,"%.3f", p->x);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"Cartn_x",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      sprintf(buffer,"%.3f", p->y);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"Cartn_y",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */


      sprintf(buffer,"%.3f", p->z);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"Cartn_z",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      /* author atom site labels                                        */
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"auth_asym_id",
                             (xmlChar *)p->chain))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */


      strcpy(buffer,p->atnam);
      KILLTRAILSPACES(buffer);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"auth_atom_id",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      strcpy(buffer,p->resnam);
      KILLTRAILSPACES(buffer);
      KILLLEADSPACES(buffer_ptr,buffer);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"auth_comp_id",
                             (xmlChar *)buffer_ptr))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      
      sprintf(buffer,"%d", p->resnum);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"auth_seq_id",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */


      /* record type atom/hetatm                                        */
      strcpy(buffer,p->record_type);
      KILLTRAILSPACES(buffer);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"group_PDB",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */


      /* atom site labels                                               */
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"label_alt_id",
                             NULL))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      if(p->altpos == ' ')
      {
         xmlNewNsProp(node, xsi, (xmlChar *)"nil", (xmlChar *)"true");
      }
      else
      {
         buffer[0] = p->altpos;
         buffer[1] = '\0';
         xmlNodeSetContent(node, (xmlChar *)buffer);
      }
      
      if((node = xmlNewChild(atom_node, NULL, 
                             (xmlChar *)"label_asym_id",
                             (xmlChar *)p->chain))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      strcpy(buffer,p->atnam);
      KILLTRAILSPACES(buffer);
      if((node = xmlNewChild(atom_node, NULL, 
                             (xmlChar *)"label_atom_id",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      strcpy(buffer,p->resnam);
      KILLTRAILSPACES(buffer);
      KILLLEADSPACES(buffer_ptr,buffer);
      if((node = xmlNewChild(atom_node, NULL, 
                             (xmlChar *)"label_comp_id",
                             (xmlChar *)buffer_ptr))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      /* Note: Entity ID is not stored in PDB data structure. 
               Value set to 1 
      */
      if((node = xmlNewChild(atom_node, NULL,
                             (xmlChar *)"label_entity_id",
                             (xmlChar *)"1"))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */
      
      sprintf(buffer,"%d", p->resnum);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"label_seq_id",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      /* occupancy                                                      */
      sprintf(buffer,"%.2f", p->occ);
      if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"occupancy",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */
                         
      /* insertion code
         Note: Insertion code node only included for residues with 
               insertion codes 
      */
      if(strcmp(p->insert," "))
      {
         sprintf(buffer,"%s", p->insert);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_PDB_ins_code",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);                   /* 25.02.15                  */
      }

      /* model number
         Note: Model number is not stored in PDB data structure.
               Value set to 1
      */
      if((node = xmlNewChild(atom_node, NULL,
                             (xmlChar *)"pdbx_PDB_model_num",
                             (xmlChar *)"1"))==NULL)
         XMLDIE(doc);                      /* 25.02.15                  */

      /* formal charge
         Note: Formal charge node not included for neutral atoms 
      */
      if(p->formal_charge != 0)
      {
         sprintf(buffer,"%d", p->formal_charge);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_formal_charge",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);                   /* 25.02.15                  */
      }

      /* atom symbol
         Note: If the atomic symbol is not set in PDB data structure then
               the value set is based on columns 13-14 of pdb-formated
               text file.  
      */
      sprintf(buffer,"%s", p->element);
      KILLLEADSPACES(buffer_ptr,buffer);
      if(strlen(buffer_ptr))
      {
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"type_symbol",
                                (xmlChar *)buffer_ptr))==NULL)
            XMLDIE(doc);                   /* 25.02.15                  */

      }
      else
      {
         blSetElementSymbolFromAtomName(buffer,p->atnam_raw);
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"type_symbol",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);                   /* 25.02.15                  */
      }

      /* Segment ID 
         Note: Segment ID is not included if blank 
      */
      if(strncmp(p->segid, "    ", 4))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"seg_id",
                                (xmlChar *)p->segid))==NULL)
            XMLDIE(doc);                   /* 25.02.15                  */
      }
   }

   /* Write to doc file pointer                                         */
   xmlDocFormatDump(fp,doc,1);

   /* Free Memory                                                       */
   xmlFreeDoc(doc);
   xmlCleanupParser();

   return(TRUE);

#endif
}


/************************************************************************/
/*>void blSetElementSymbolFromAtomName(char *element, char *atom_name)
   -------------------------------------------------------------------
*//**

   \param[out]    *element       Element symbol
   \param[in]     *atom_name     Atom name

   Set the element symbol (columns 77-78 of a pdb file) based on the 
   atom name (columns 13-16 of a pdb file). 

   The atom name is stored in the PDB data stucture as atnam_raw.

-  17.07.14 Original. By: CTP

*/
void blSetElementSymbolFromAtomName(char *element, char *atom_name)
{
   char  buffer[] = "  ",
         *buffer_ptr;

   /* copy first two letters of atom name to buffer                     */
   strncpy(buffer,atom_name,2);
   buffer[2] = '\0';
   KILLLEADSPACES(buffer_ptr,buffer);

   /* remove digits                                                     */
   if(strlen(buffer_ptr) == 2 && isdigit(buffer[1]))
   {
      buffer[1] = '\0';
   }
   else if(strlen(buffer_ptr) == 2 && !isalpha(buffer[0]))
   {
      buffer_ptr += 1;
   }

   /* fix hydrogens and carbons                                         */
   if(strlen(buffer_ptr) == 2 && atom_name[3] != ' ' &&
      (atom_name[0] == 'H' || atom_name[0] == 'C' || 
       atom_name[0] == 'N' || atom_name[0] == 'O' || 
       atom_name[0] == 'P'))
   {
         if(!isalpha(atom_name[2]) || !isalpha(atom_name[3]))
         {
            buffer[1] = '\0';
         }
   }

   /* copy atom symbol to element                                       */
   strcpy(element, buffer_ptr);
}


/************************************************************************/
/*>void blWriteGromosPDB(FILE *fp, PDB *pdb)
   -----------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list by calls to WritePDBRecord()

-  08.03.89 Original
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 Uses NEXT macro; void type
-  08.07.93 Added insertion of TER cards
-  22.02.94 And a TER card at the end of the file
-  15.02.01 This is the old WritePDB()
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
-  23.02.15 Proper TER card support  By: ACRM
-  24.02.15 Now just uses blWritePDBAsPDBorGromos()
*/
void blWriteGromosPDB(FILE *fp,
                      PDB  *pdb)
{
   blWritePDBAsPDBorGromos(fp, pdb, TRUE);
}


/************************************************************************/
/*>void blWriteGromosPDBRecord(FILE *fp, PDB *pdb)
   -----------------------------------------------
*//**

   \param[in]     *fp     PDB file pointer to be written
   \param[in]     *pdb    PDB linked list record to write

   Write a PDB record

-  08.03.89 Original
-  28.03.90 Changed to match ReadPDB() V1.2 for column widths
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 void type
-  22.06.93 Changed to %lf. Ljust strings
-  11.03.94 %lf back to %f (!)
-  12.02.01 This is the old WritePDBRecord()
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blWriteGromosPDBRecord(FILE *fp,
                            PDB  *pdb)
{
   fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
           pdb->record_type,
           pdb->atnum,
           pdb->atnam,
           pdb->resnam,
           pdb->chain,
           pdb->resnum,
           pdb->insert,
           pdb->x,
           pdb->y,
           pdb->z,
           pdb->occ,
           pdb->bval);
}


/************************************************************************/
/*>BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
   ----------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes a PDB file including header and trailer information.
   Output in PDBML-format if flags set.

-  21.06.14 Original   By: CTP
-  18.08.14 Added XML_SUPPORT option. Return error if attempting to write 
            PDBML format. By: CTP
-  12.02.15 Now a wrapper to blDoWriteWhilePDB()
-  24.02.15 blWriteAsPDBML() changed to blWritePDBAsPDBML()
            blWriteAsPDB() changed to blWritePDBAsPDBorGromos()
-  25.02.15 No longer a wrapper
*/
BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
{
   int nter;

   if((gPDBXMLForce == FORCEXML_XML) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == TRUE))
   {
#ifdef XML_SUPPORT
      /* Write PDBML file (omitting header and footer data)             */
      blWritePDBAsPDBML(fp, wpdb->pdb);
#else
      /* PDBML not supported                                            */
      return(FALSE);
#endif
   }
   else
   {
      /* Check format                                                   */
      if(blFormatCheckWritePDB(wpdb->pdb) == FALSE)
      {
         return(FALSE);
      }

      /* Write whole PDB File                                           */
      blWriteWholePDBHeader(fp, wpdb);
      nter = blWritePDBAsPDBorGromos(fp, wpdb->pdb, FALSE);
      blWriteWholePDBTrailer(fp, wpdb, nter);
   }
   
   return(TRUE);
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
-  12.02.15  Added XML check   By: ACRM
*/
void blWriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;

   if((gPDBXMLForce != FORCEXML_XML) && (gPDBXML == FALSE))
   {
      for(s=wpdb->header; s!=NULL; NEXT(s))
      {
         fputs(s->string, fp);
      }
   }
}


/************************************************************************/
/*>void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb, int numTer)
   -----------------------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer
   \param[in]     numTer     Number of TER cards

   Writes the trailer of a PDB file 

-  30.05.02  Original   By: ACRM
-  21.06.14  Renamed to blWriteWholePDBTrailer() By: CTP
-  12.02.15  Added XML check   By: ACRM
-  18.02.15  Complete rewrite to use the parsed CONECT data rather than
             simply rewriting what was read in
-  23.02.15  Added numTer parameter
-  02.03.15  Padded END and CONECT
*/
void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb, int numTer)
{
   int nConect = 0;
   
   if((gPDBXMLForce != FORCEXML_XML) && (gPDBXML == FALSE))
   {
      /* Write the CONECT records                                       */
      PDB *p;
      for(p=wpdb->pdb; p!=NULL; NEXT(p))
      {
         if(p->nConect)
         {
            BOOL conectPrinted = FALSE;
            int  i, 
                 width=0;
            char format[8];

            for(i=0; i<p->nConect; i++)
            {
               if(!(i%4))
               {
                  if(conectPrinted)
                     fprintf(fp, "\n");
                  fprintf(fp,"CONECT%5d", p->atnum);
                  nConect++;
                  width = 11;
                  conectPrinted = 1;
               }
               fprintf(fp,"%5d", p->conect[i]->atnum);
               width += 5;
            }
            sprintf(format, "%%%ds\n", 80-width);
            fprintf(fp, format, " ");
         }
      }

      /* Now write the MASTER record                                    */
      WriteMaster(fp, wpdb, nConect, numTer);
      
      fprintf(fp, "END%77s\n", " ");
   }
}


/************************************************************************/
/*>static void WriteMaster(FILE *fp, WHOLEPDB *wpdb, int numConect,
                           int numTer)
   ----------------------------------------------------------------
*//**
   \param[in]  *fp          File pointer for output
   \param[in]  *wpdb        WHOLEPDB structure
   \param[in]  numConect    Number of CONECT records
   \param[in]  numTer       Number of TER cards

   Writes the PDB MASTER record. Everything is calculated in this code
   apart from the number of CONECT records and the number of TER cards
   which must be passed in

-  22.02.15 Original   By: ACRM
-  02.03.15 Corrected counting of ORIGX, SCALE and MTRIX
*/
static void WriteMaster(FILE *fp, WHOLEPDB *wpdb, int numConect,
                        int numTer)
{
   int numRemark = 0, /* Number of REMARK records                       */
       numHet    = 0, /* Number of HET records                          */
       numHelix  = 0, /* Number of HELIX records                        */
       numSheet  = 0, /* Number of SHEET records                        */
       numTurn   = 0, /* Number of TURN records                         */
       numSite   = 0, /* Number of SITE records                         */
       numXform  = 0, /* Number of coordinate transformation records 
                         (ORIGX+SCALE+MTRIX)                            */
       numCoord  = 0, /* Number of atomic coordinate records 
                         (ATOM+HETATM)                                  */
       numSeq    = 0; /* Number of SEQRES records                       */
   STRINGLIST *s;
   PDB *p;
   

   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "REMARK", 6))
         numRemark++;
      if(!strncmp(s->string, "HET   ", 6))
         numHet++;
      if(!strncmp(s->string, "HELIX ", 6))
         numHelix++;
      if(!strncmp(s->string, "SHEET ", 6))
         numSheet++;
      if(!strncmp(s->string, "TURN  ", 6))
         numTurn++;
      if(!strncmp(s->string, "SITE  ", 6))
         numSite++;
      if(!strncmp(s->string, "ORIGX",  5))
         numXform++;
      if(!strncmp(s->string, "SCALE",  5))
         numXform++;
      if(!strncmp(s->string, "MTRIX",  5))
         numXform++;
      if(!strncmp(s->string, "SEQRES", 6))
         numSeq++;
   }

   for(p=wpdb->pdb; p!=NULL; NEXT(p))
      numCoord++;
   
   fprintf(fp,"MASTER    %5d    0%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",
           numRemark,
           numHet,
           numHelix,
           numSheet,
           numTurn,
           numSite,
           numXform,
           numCoord,
           numTer,
           numConect,
           numSeq);
}


/************************************************************************/
/*>void blWriteWholePDBHeaderNoRes(FILE *fp, WHOLEPDB *wpdb)
   ---------------------------------------------------------
*//**

   \param[in]     *fp        File pointer
   \param[in]     *wpdb      Whole PDB structure pointer

   Writes the header of a PDB file, but skips any records that include
   residue numbers for cases where these may have changed.

   Skips: 
      REMARK 500
      DBREF
      HELIX
      SHEET
      SSBOND
      CISPEP

-  02.03.15  Original   By: ACRM
*/
void blWriteWholePDBHeaderNoRes(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   int i;
   char *recordLabels[] = {"REMARK 500", 
                           "DBREF ", 
                           "HELIX ", 
                           "SHEET ", 
                           "SSBOND", 
                           "CISPEP",
                           NULL};
   

   if((gPDBXMLForce != FORCEXML_XML) && (gPDBXML == FALSE))
   {
      for(s=wpdb->header; s!=NULL; NEXT(s))
      {
         BOOL printIt = TRUE;
         for(i=0; recordLabels[i] != NULL; i++)
         {
            if(!strncmp(s->string, recordLabels[i], strlen(recordLabels[i])))
            {
               printIt = FALSE;
               break;
            }
         }
         if(printIt)
            fputs(s->string, fp);
      }
   }
}



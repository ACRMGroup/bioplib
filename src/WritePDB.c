/************************************************************************/
/**

   \file       WritePDB.c
   
   \version    V1.30
   \date       09.02.18
   \brief      Write a PDB file from a linked list
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2018
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
                  Added blWriteWholePDBHeaderNoRes()
                  Added space padding to TER in blWriteTerCard()
                  Added space padding to END and to CONECT in
                  blWriteWholePDBTrailer()
-  V1.22 09.03.15 blWriteWholePDBHeaderNoRes() now skips more header lines
-  V1.23 10.03.15 blWriteWholePDBHeaderNoRes() now skips more header lines
-  V1.24 02.04.15 WriteMaster() Padded MASTER record to 80 cols.  By: CTP
-  V1.25 13.05.15 Updated blWritePDBAsPDBML() to write CONECT records, 
                  Added output of COMPND and SOURCE records.  By: CTP
-  V1.26 08.06.15 Updated XML check for output format.  By: CTP
-  V1.27 21.06.15 Added blMapChainsToEntity(). blDoWritePDBAsPDBML() sets
                  entity_id based on COMPND records if entty_id not 
                  already set.  By: CTP
-  V1.28 10.07.15 Added return value for blWritePDBAsPDBML() and 
                  blDoWritePDBAsPDBML() when XML_SUPPORT not defined 
                  By: ACRM
-  V1.29 29.07.15 Added output of PDBML seqres records from wpdb->header.
                  Added ReadSeqresChainLabelWholePDB() and 
                  ReadSeqresResidueListWholePDB().  By: CTP
-  V1.30 09.02.18 Corrected calls to ABS() in formal charges By: ACRM

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

   #FUNCTION blWritePDBAsPDBorGromos()
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

   #FUNCTION  blWriteTerCard()
   Prints a complete new PDB format TER card
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
#include "hash.h"
#include "fsscanf.h"

/************************************************************************/
/* Prototypes
*/
static void WriteMaster(FILE *fp, WHOLEPDB *wpdb, int numConect,
                        int numTer);
static BOOL blDoWritePDBAsPDBML(FILE *fp, WHOLEPDB  *wpdb, BOOL doWhole);
static BOOL blSetPDBMLDateField(char *pdbml_date, char *pdb_date);
static HASHTABLE *blMapChainsToEntity(WHOLEPDB *wpdb);
static char **ReadSeqresChainLabelWholePDB(WHOLEPDB *wpdb, int *nchains);
static STRINGLIST **ReadSeqresResidueListWholePDB(WHOLEPDB *wpdb, 
                                                  int *nchains);
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
   BOOL  doneTer = FALSE;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If previous was non-null and was an ATOM                       */
      if((prev!=NULL) && !strncmp(prev->record_type, "ATOM  ", 6))
      {
         /* If the chain has changed write a TER card                   */
         if(!doneTer && !CHAINMATCH(p->chain, prev->chain))
         {
            blWriteTerCard(fp, prev);
            numTer++;
            doneTer = TRUE;
         }
         /* If we've moved into a HETATM, then see if this HETATOM group
            is bonded to the previous residue. If it isn't, print a TER
            card
         */
         if(!doneTer && !strncmp(p->record_type, "HETATM", 6))
         {
            PDB *prevStart = blFindResidue(pdb, prev->chain, prev->resnum,
                                           prev->insert);
            if(!blAreResiduePointersBonded(prevStart, p, (REAL)0.2))
            {
               blWriteTerCard(fp, prev);
               numTer++;
               doneTer = TRUE;
            }
         }
      }

      if(doGromos)
      {
         blWriteGromosPDBRecord(fp,p);
         doneTer = FALSE;
      }
      else
      {
         blWritePDBRecord(fp,p);
         doneTer = FALSE;
      }
      prev=p;
   }

   if((!doneTer)   && 
      (prev!=NULL) && 
      !strncmp(prev->record_type, "ATOM  ", 6))
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
-  09.02.18 Corrected ABS() call
*/
void blWritePDBRecord(FILE *fp,
                      PDB  *pdb)
{
   char charge = ' ',
        sign   = ' ';

   if(pdb->formal_charge && ABS(pdb->formal_charge) <= 8)
   {
      charge = (char)('0' + ABS(pdb->formal_charge));
      sign   = (char)(pdb->formal_charge > 0 ? '+':'-');
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
-  09.02.17 Corrected ABS() call
*/
void blWritePDBRecordAtnam(FILE *fp,
                           PDB  *pdb)
{
   char charge = ' ',
        sign   = ' ';

   if(pdb->formal_charge && ABS(pdb->formal_charge) <= 8)
   {
      charge = (char)('0' + ABS(pdb->formal_charge));
      sign   = (char)(pdb->formal_charge > 0 ? '+':'-');
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
-  29.04.15 Updated to write CONECT records.  By: CTP
-  11.05.15 Made function into wrapper for blDoWritePDBAsPDBML(). By: CTP
-  10.07.15 Added return value for no XML_SUPPORT  By: ACRM

*/
BOOL blWritePDBAsPDBML(FILE *fp, PDB  *pdb)
{
#ifndef XML_SUPPORT

   /* PDBML format not supported.                                       */
   return(FALSE);                         /* 10.07.15                   */

#else 

   /* PDBML format supported.                                           */
   WHOLEPDB wpdb;
   wpdb.header  = NULL;
   wpdb.trailer = NULL;
   wpdb.natoms  =    0;
   wpdb.pdb     =  pdb;
   return(blDoWritePDBAsPDBML(fp, &wpdb, FALSE));

#endif
}

/************************************************************************/
/*>BOOL blDoWritePDBAsPDBML(FILE *fp, WHOLEPDB  *wpdb, BOOL doWhole)
   -----------------------------------------------------------------
*//**

   \param[in]     *fp      PDB file pointer to be written
   \param[in]     *wpdb    WHOLEPDB to write
   \param[in]     doWhole  Write whole pdb including header and conect 
                           records or just coordinate records.
   \return                 Success

   Write a PDB linked list in PDBML format.

-  02.06.14 Original. By: CTP
-  21.06.14 Renamed blWriteAsPDBML() and updated symbol handling. By: CTP
-  17.07.14 Use blSetElementSymbolFromAtomName() By: CTP
-  16.08.14 Use element and charge data. By: CTP
-  17.02.15 Added segid support  By: ACRM
-  24.02.15 Changed name to blWritePDBAsPDBML()
-  25.02.15 Changed to type BOOL and checks all memory allocations
-  29.04.15 Updated to write CONECT records.  By: CTP
-  11.05.15 Renamed from blWritePDBAsPDBML() to blDoWritePDBAsPDBML().
            Changed to a static function.
            blWritePDBAsPDBML() is now a wrapper for this function.
            This function takes WHOLEPDB as input instead of PDB and 
            writes wpdb->header and wpdb->trailer info if doWhole param is
            TRUE.  By: CTP
-  21.06.15 Write entity_id. Use COMPND record to set entity_id if not set
            in PDB. Set compound type to polymer. By:  CTP
-  10.07.15 Added return value for no XML_SUPPORT  By: ACRM
-  29.07.15 Added output of SEQRES records from wpdb->header.  By: CTP
*/
static BOOL blDoWritePDBAsPDBML(FILE *fp, WHOLEPDB  *wpdb, BOOL doWhole)
{
#ifndef XML_SUPPORT

   /* PDBML format not supported.                                       */
   return(FALSE);                         /* 10.07.15                   */

#else 

   /* PDBML format supported                                            */
   PDB         *p,
               *q;
   xmlDocPtr   doc         = NULL;
   xmlNodePtr  root_node   = NULL, 
               sites_node  = NULL, 
               atom_node   = NULL, 
               node        = NULL;
   xmlNsPtr    pdbx        = NULL,
               xsi         = NULL;
   char        buffer[16], 
               *buffer_ptr;
   int         conect_id   = 0,
               i, j;

   char header[82]     =   "",
        date[82]       =   "",
        pdbcode[82]    =   "",
        pdbml_date[11] =   "",
        *title         = NULL;

   COMPND    compound;
   PDBSOURCE species;
   int molid = 0;
   HASHTABLE *chain_to_entity = NULL;
   int        seqres_nchains    =    0;
   char       **seqres_chain    = NULL;
   STRINGLIST **seqres_residues = NULL,
              *s = NULL;

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
   
   /* Write Coordinate Data                                             */
   
   /* map chain to entity from compnd records                           */
   chain_to_entity = blMapChainsToEntity(wpdb);

   /* Atom_sites node                                                   */
   if((sites_node = xmlNewChild(root_node, NULL, 
                                (xmlChar *)"atom_siteCategory", 
                                NULL))==NULL)
      XMLDIE(doc);                         /* 25.02.15                  */

   
   /* Atom nodes                                                        */
   for(p=wpdb->pdb; p!=NULL; NEXT(p))
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

      /* Note: Entity ID is not set for PDB format.
               If entity_id is not set in PDB list then set from COMPND
               or default to 1.
      */
      if(p->entity_id)
      {
         /* set from PDB list */
         sprintf(buffer,"%d", p->entity_id);
      }
      else if(chain_to_entity != NULL && 
              blHashKeyDefined(chain_to_entity, p->chain))
      {
         /* set from COMPND record */
         sprintf(buffer,"%d", 
                 blGetHashValueInt(chain_to_entity, p->chain));
      }
      else
      {
         /* default */
         strcpy(buffer,"1");
      }
      if((node = xmlNewChild(atom_node, NULL,
                             (xmlChar *)"label_entity_id",
                             (xmlChar *)buffer))==NULL)
         XMLDIE(doc);


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

   /* Finished Coordinate Data                                          */
   /* Clean up and return if doWhole == FALSE                           */
   if(doWhole == FALSE)
   {
      /* Write to doc file pointer                                      */
      xmlDocFormatDump(fp,doc,1);

      /* Free Memory                                                    */
      xmlFreeDoc(doc);
      xmlCleanupParser();
      blFreeHash(chain_to_entity);

      return(TRUE);
   }


   /*** Write Header and trailer data in PDBML-format                 ***/

   /* struct_conn node                                                  */
   for(p=wpdb->pdb; p!=NULL; NEXT(p))
   {
      if(p->nConect)
      {
         if((sites_node = xmlNewChild(root_node, NULL, 
                                      (xmlChar *)"struct_connCategory", 
                                      NULL))==NULL)
            XMLDIE(doc);

         break;
      }
   }
   
   /* Conect nodes                                                      */
   for(p=wpdb->pdb; p!=NULL; NEXT(p))
   {
      /* skip TER                                                       */
      if(!strncmp("TER",p->resnam,3))
      {
         continue;
      }

      /* Add conect nodes                                               */
      for(i=0; i < p->nConect; i++)
      {
         if((atom_node = xmlNewChild(sites_node, NULL,
                                     (xmlChar *)"struct_conn", NULL))==NULL)
            XMLDIE(doc);

         /* set conect atoms */
         q = p->conect[i];
         conect_id++;         

         /* conect id */
         sprintf(buffer, "%s%d", "conect", conect_id);
         xmlNewProp(atom_node, (xmlChar *)"id", (xmlChar *)buffer);

         /* connection type */
         sprintf(buffer, "%s", "covale"); /* set type to covale         */
         if((node = xmlNewChild(atom_node, NULL,
                                (xmlChar *)"conn_type_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         /* bond length                                                 */
         sprintf(buffer, "%.3f", DIST(p,q));
         if((node = xmlNewChild(atom_node, NULL,
                                (xmlChar *)"pdbx_dist_value",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         /* atom one data                                                */
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr1_auth_asym_id",
                                (xmlChar *)p->chain))==NULL)
            XMLDIE(doc);

         strcpy(buffer,p->resnam);
         KILLTRAILSPACES(buffer);
         KILLLEADSPACES(buffer_ptr,buffer);
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr1_auth_comp_id",
                                (xmlChar *)buffer_ptr))==NULL)
            XMLDIE(doc);

         sprintf(buffer,"%d", p->resnum);
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr1_auth_seq_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         /* include alt_id if present                                   */
         if(p->altpos != ' ')
         {
            if((node = xmlNewChild(atom_node, NULL,
                                   (xmlChar *)"ptnr1_label_alt_id",
                                   NULL))==NULL)
               XMLDIE(doc);

            buffer[0] = p->altpos;
            buffer[1] = '\0';
            xmlNodeSetContent(node, (xmlChar *)buffer);
         }
         
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr1_label_asym_id",
                                (xmlChar *)p->chain))==NULL)
            XMLDIE(doc);

         strcpy(buffer,p->atnam);
         KILLTRAILSPACES(buffer);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr1_label_atom_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         strcpy(buffer,p->resnam);
         KILLTRAILSPACES(buffer);
         KILLLEADSPACES(buffer_ptr,buffer);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr1_label_comp_id",
                                (xmlChar *)buffer_ptr))==NULL)
            XMLDIE(doc);
      
         sprintf(buffer,"%d", p->resnum);
         if((node = xmlNewChild(atom_node, NULL,
                                (xmlChar *)"ptnr1_label_seq_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);
                   
         /* insertion code
            Note: Insertion code node only included for residues with 
                  insertion codes 
         */
         if(strcmp(p->insert," "))
         {
            sprintf(buffer,"%s", p->insert);
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"ptnr1_PDB_ins_code",
                                   (xmlChar *)buffer))==NULL)
               XMLDIE(doc);
         }


         /* atom two data                                               */
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr2_auth_asym_id",
                                (xmlChar *)q->chain))==NULL)
            XMLDIE(doc);

         strcpy(buffer,q->resnam);
         KILLTRAILSPACES(buffer);
         KILLLEADSPACES(buffer_ptr,buffer);
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr2_auth_comp_id",
                                (xmlChar *)buffer_ptr))==NULL)
            XMLDIE(doc);

         sprintf(buffer,"%d", q->resnum);
         if((node = xmlNewChild(atom_node, NULL, (xmlChar *)"ptnr2_auth_seq_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         /* include alt_id if present                                   */
         if(q->altpos != ' ')
         {
            if((node = xmlNewChild(atom_node, NULL,
                                   (xmlChar *)"ptnr2_label_alt_id",
                                   NULL))==NULL)
               XMLDIE(doc);

            buffer[0] = q->altpos;
            buffer[1] = '\0';
            xmlNodeSetContent(node, (xmlChar *)buffer);
         }

         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr2_label_asym_id",
                                (xmlChar *)q->chain))==NULL)
            XMLDIE(doc);

         strcpy(buffer,q->atnam);
         KILLTRAILSPACES(buffer);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr2_label_atom_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

         strcpy(buffer,q->resnam);
         KILLTRAILSPACES(buffer);
         KILLLEADSPACES(buffer_ptr,buffer);
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"ptnr2_label_comp_id",
                                (xmlChar *)buffer_ptr))==NULL)
            XMLDIE(doc);
      
         sprintf(buffer,"%d", q->resnum);
         if((node = xmlNewChild(atom_node, NULL,
                                (xmlChar *)"ptnr2_label_seq_id",
                                (xmlChar *)buffer))==NULL)
            XMLDIE(doc);

                         
         /* insertion code
            Note: Insertion code node only included for residues with 
                  insertion codes 
         */
         if(strcmp(q->insert," "))
         {
            sprintf(buffer,"%s", q->insert);
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"ptnr2_PDB_ins_code",
                                   (xmlChar *)buffer))==NULL)
               XMLDIE(doc);
         }
      }      
   }

   /* get header data                                                   */
   blGetHeaderWholePDB(wpdb, header, 82, date, 82, pdbcode, 82);
   KILLTRAILSPACES(pdbcode);
   
   /* get title                                                         */
   title = blGetTitleWholePDB(wpdb);

   /* add date node                                                     */
   if(blSetPDBMLDateField(pdbml_date, date))
   {
      if((sites_node = xmlNewChild(root_node, NULL, 
                                   (xmlChar *)"database_PDB_revCategory", 
                                   NULL))==NULL)
         XMLDIE(doc);
            
      if((atom_node = xmlNewChild(sites_node, NULL, 
                                  (xmlChar *)"database_PDB_rev",
                                  NULL))==NULL)
         XMLDIE(doc);
      xmlNewProp(atom_node, (xmlChar *)"num", (xmlChar *)"1");

      if((node = xmlNewChild(atom_node, NULL, 
                            (xmlChar *)"date",
                            (xmlChar *)pdbml_date))==NULL)
         XMLDIE(doc);
 
       if((node = xmlNewChild(atom_node, NULL, 
                            (xmlChar *)"date_original",
                            (xmlChar *)pdbml_date))==NULL)
         XMLDIE(doc);

       if((node = xmlNewChild(atom_node, NULL, 
                            (xmlChar *)"mod_type",
                            (xmlChar *)"0"))==NULL)
         XMLDIE(doc);
     
      if(strlen(pdbcode) == 4)
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"replaces",
                                (xmlChar *)pdbcode))==NULL)
            XMLDIE(doc);
      }
   }

   /* add compnd nodes                                                  */
   sites_node = NULL;
   for(i=1; blGetCompoundWholePDBMolID(wpdb, i, &compound); i++)
   {
      molid = i;

      if(sites_node == NULL)
      {
         /* add COMPND node */
         if((sites_node = xmlNewChild(root_node, NULL, 
                                      (xmlChar *)"entityCategory", 
                                      NULL))==NULL)
            XMLDIE(doc);
      }

      if((atom_node = xmlNewChild(sites_node, NULL, 
                                  (xmlChar *)"entity",
                                  NULL))==NULL){ XMLDIE(doc); }

      sprintf(buffer,"%d", i);
      xmlNewProp(atom_node, (xmlChar *)"id", (xmlChar *)buffer);

      if(strlen(compound.other))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"details",
                                (xmlChar *)compound.other))==NULL)
         { XMLDIE(doc); }
      }

      if(strlen(compound.molecule))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_description",
                                (xmlChar *)compound.molecule))==NULL)
         { XMLDIE(doc); }
      }

      if(strlen(compound.ec))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_ec",
                                (xmlChar *)compound.ec))==NULL)
         { XMLDIE(doc); }
      }

      if(strlen(compound.fragment))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_fragment",
                                (xmlChar *)compound.fragment))==NULL)
         { XMLDIE(doc); }
      }

      if(strlen(compound.mutation))
      {
         if((node = xmlNewChild(atom_node, NULL, 
                                (xmlChar *)"pdbx_mutation",
                                (xmlChar *)compound.mutation))==NULL)
         { XMLDIE(doc); }
      }

      /* set type to polymer */
      if((node = xmlNewChild(atom_node, NULL, 
                             (xmlChar *)"type",
                             (xmlChar *)"polymer"))==NULL)
      { XMLDIE(doc); }
   }

   /* add source nodes                                                  */
   sites_node = NULL;
   j = 0;
   for(i=1; i <= molid; i++)
   {
      if(blGetSpeciesWholePDBMolID(wpdb, i, &species))
      {
         if(sites_node == NULL)
         {
            /* add SOURCE node */
            if((sites_node = xmlNewChild(root_node, NULL, 
                                         (xmlChar *)"entity_src_genCategory", 
                                         NULL))==NULL)
               XMLDIE(doc);
         }

         j++; /* pdbx_src_id */
         
         if((atom_node = xmlNewChild(sites_node, NULL, 
                                  (xmlChar *)"entity_src_gen",
                                  NULL))==NULL){ XMLDIE(doc); }

         sprintf(buffer,"%d", i);
         xmlNewProp(atom_node,(xmlChar *)"entity_id",(xmlChar *)buffer);
         sprintf(buffer,"%d", j);
         xmlNewProp(atom_node,(xmlChar *)"pdbx_src_id",(xmlChar *)buffer);

         if(strlen(species.commonName))
         {
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"pdbx_gene_src_common_name",
                                   (xmlChar *)species.commonName))==NULL)
            { XMLDIE(doc); }
         }

         if(strlen(species.strain))
         {
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"pdbx_gene_src_strain",
                                   (xmlChar *)species.strain))==NULL)
            { XMLDIE(doc); }
         }

         if(species.taxid != 0)
         {
            sprintf(buffer,"%d", species.taxid);
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"pdbx_gene_src_ncbi_taxonomy_id",
                                   (xmlChar *)buffer))==NULL)
            { XMLDIE(doc); }
         }

         if(strlen(species.scientificName))
         {
            if((node = xmlNewChild(atom_node, NULL, 
                                   (xmlChar *)"pdbx_gene_src_scientific_name",
                                   (xmlChar *)species.scientificName))==NULL)
            { XMLDIE(doc); }
         }
      }
   }


   /* SEQRES nodes                                                      */
   /* get seqres chain ids from wpdb->header                            */
   seqres_chain = ReadSeqresChainLabelWholePDB(wpdb, &seqres_nchains);
   if(seqres_chain)
   {
      /* get seqres residues from wpdb->header                          */
      seqres_residues = ReadSeqresResidueListWholePDB(wpdb,
                                                      &seqres_nchains);
   
      if(seqres_residues)
      {
         /* add pdbx_poly_seq_schemeCategory node                       */
         if((sites_node = xmlNewChild(root_node, NULL, 
                           (xmlChar *)"pdbx_poly_seq_schemeCategory", 
                           NULL))==NULL){ XMLDIE(doc); }

         /* cycle through chains                                        */
         for(i=0;i<seqres_nchains;i++)
         {
            int res    = 1,                 /* reset residue count      */
                entity = 1;                 /* set entity_id to default */

            /* set entity based on chain id                             */
            if(chain_to_entity != NULL && 
               blHashKeyDefined(chain_to_entity, seqres_chain[i]))
            {
               entity = blGetHashValueInt(chain_to_entity, 
                                          seqres_chain[i]);
            }

             /* cycle through residues                                  */
            for(s=seqres_residues[i];s!=NULL;NEXT(s),res++)
            {
              /* add pdbx_poly_seq_scheme node                          */
              if((atom_node = xmlNewChild(sites_node, NULL,
                              (xmlChar *)"pdbx_poly_seq_scheme", NULL))
                                                                 == NULL)
              { XMLDIE(doc); }
              /* add attributes                                         */
              /* asym_id                                                */
              xmlNewProp(atom_node, (xmlChar *)"asym_id", 
                         (xmlChar *)seqres_chain[i]);
              /* entity_id                                              */
              sprintf(buffer, "%d", entity);
              xmlNewProp(atom_node, (xmlChar *)"entity_id", 
                         (xmlChar *)buffer);
              /* mon_id                                                 */
              xmlNewProp(atom_node, (xmlChar *)"mon_id", 
                         (xmlChar *)s->string);
              /* seq_id                                                 */
              sprintf(buffer, "%d", res);
              xmlNewProp(atom_node, (xmlChar *)"seq_id", 
                         (xmlChar *)buffer);

              /* add subnodes                                           */
              /* auth_mon_id                                            */
              sprintf(buffer,"%d", res);
              if((node = xmlNewChild(atom_node, NULL, 
                                     (xmlChar *)"auth_mon_id",
                                     (xmlChar *)s->string))==NULL)
                 XMLDIE(doc);
              /* ndb_seq_num                                            */
              sprintf(buffer,"%d", res);
              if((node = xmlNewChild(atom_node, NULL, 
                                     (xmlChar *)"ndb_seq_num",
                                     (xmlChar *)buffer))==NULL)
                 XMLDIE(doc);
              /* pdb_mon_id                                             */
              sprintf(buffer,"%d", res);
              if((node = xmlNewChild(atom_node, NULL, 
                                     (xmlChar *)"pdb_mon_id",
                                     (xmlChar *)s->string))==NULL)
                 XMLDIE(doc);
              /* pdb_strand_id                                          */
              sprintf(buffer,"%d", res);
              if((node = xmlNewChild(atom_node, NULL, 
                                     (xmlChar *)"pdb_strand_id",
                                     (xmlChar *)seqres_chain[i]))==NULL)
                 XMLDIE(doc);
            }

            FREELIST(seqres_residues[i],STRINGLIST);   /* free residues */
            free(seqres_chain[i]);                     /* free chain id */
         }
         free(seqres_residues);                  /* free residues array */
      }
      free(seqres_chain);                        /* free chain id array */
   }


   /* pdb entry */
   if(strlen(pdbcode))
   {
      if((sites_node = xmlNewChild(root_node, NULL, 
                                   (xmlChar *)"entryCategory", 
                                   NULL))==NULL)
         XMLDIE(doc);
            
      if((atom_node = xmlNewChild(sites_node, NULL, 
                                  (xmlChar *)"entry",
                                  (xmlChar *)""))==NULL)
         XMLDIE(doc);
      xmlNewProp(atom_node, (xmlChar *)"id", (xmlChar *)pdbcode);
   }

   /* title node */
   if(title != NULL && strlen(pdbcode))
   {
      if((sites_node = xmlNewChild(root_node, NULL, 
                                   (xmlChar *)"structCategory", 
                                   NULL))==NULL)
         XMLDIE(doc);
            
     if((atom_node = xmlNewChild(sites_node, NULL, 
                                 (xmlChar *)"struct",
                                 NULL))==NULL)
        XMLDIE(doc);
     xmlNewProp(atom_node, (xmlChar *)"entry_id", (xmlChar *)pdbcode);
     
     if((node = xmlNewChild(atom_node, NULL, 
                                 (xmlChar *)"title",
                                 (xmlChar *)title))==NULL)
        XMLDIE(doc);

   }
   if(title != NULL){ free(title); }

   /* header node */
   if(strlen(header) && strlen(pdbcode))
   {
      if((sites_node = xmlNewChild(root_node, NULL, 
                                   (xmlChar *)"struct_keywordsCategory", 
                                   NULL))==NULL)
         XMLDIE(doc);
            
     if((atom_node = xmlNewChild(sites_node, NULL, 
                                 (xmlChar *)"struct_keywords",
                                 NULL))==NULL)
        XMLDIE(doc);
     xmlNewProp(atom_node, (xmlChar *)"entry_id", (xmlChar *)pdbcode);
     
     if((node = xmlNewChild(atom_node, NULL, 
                                 (xmlChar *)"pdbx_keywords",
                                 (xmlChar *)header))==NULL)
        XMLDIE(doc);

   }


   /* Write to doc file pointer                                         */
   xmlDocFormatDump(fp,doc,1);

   /* Free Memory                                                       */
   xmlFreeDoc(doc);
   xmlCleanupParser();
   blFreeHash(chain_to_entity);

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
-  04.03.15 Added check on wpdb and wpdb->pdb being non-NULL
-  11.05.15 Updated to use blDoWritePDBAsPDBML().  By: CTP
*/
BOOL blWriteWholePDB(FILE *fp, WHOLEPDB *wpdb)
{
   int nter;

   if((wpdb==NULL) || (wpdb->pdb==NULL))
      return(FALSE);

   if((gPDBXMLForce == FORCEXML_XML) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == TRUE))
   {
#ifdef XML_SUPPORT
      /* Write PDBML file (including header and footer data)            */
      blDoWritePDBAsPDBML(fp, wpdb, TRUE);
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
-  06.08.15  Updated XML check. By: CTP
*/
void blWriteWholePDBHeader(FILE *fp, WHOLEPDB *wpdb)
{
   STRINGLIST *s;

   if((gPDBXMLForce == FORCEXML_PDB) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == FALSE))
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
-  06.08.15  Updated XML check. By: CTP
*/
void blWriteWholePDBTrailer(FILE *fp, WHOLEPDB *wpdb, int numTer)
{
   int nConect = 0;
   
   if((gPDBXMLForce == FORCEXML_PDB) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == FALSE))
   {
      /* Write the CONECT records                                       */
      PDB *p;
      for(p=wpdb->pdb; p!=NULL; NEXT(p))
      {
         if(p->nConect)
         {
            BOOL conectPrinted = FALSE;
            int  i, 
                 nPrinted,
                 width=0;
            char format[8];

            for(i=0, nPrinted=0; i<p->nConect; i++)
            {
               if(!(nPrinted%4))
               {
                  if(conectPrinted)
                     fprintf(fp, "\n");
                  fprintf(fp,"CONECT%5d", p->atnum);
                  nConect++;
                  width = 11;
                  conectPrinted = 1;
               }
               if(p->conect[i] != NULL)
               {
                  fprintf(fp,"%5d", p->conect[i]->atnum);
                  width += 5;
                  nPrinted++;
               }
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
-  02.04.15 Padded MASTER record to 80 columns.  By: CTP
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
   
   fprintf(fp,"MASTER    %5d    0%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d          \n",
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
      SEQADV
      MODRES
      HET
      SITE
      REMARK 465
      REMARK 470
      REMARK 475
      REMARK 480
      REMARK 525
      REMARK 610
      REMARK 615
      REMARK 620
      REMARK 630
      REMARK   3

-  02.03.15  Original   By: ACRM
-  09.03.15  Additionally skips SEQADV, MODRES, HET, SITE
-  10.03.15  Additionally skips REMARK 3,465,470,475,480,525,610,615,620,630
-  06.08.15  Updated XML check. By: CTP
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
                           "SEQADV",
                           "MODRES",
                           "HET   ",
                           "SITE  ",
                           "REMARK 465",
                           "REMARK 470",
                           "REMARK 475",
                           "REMARK 480",
                           "REMARK 525",
                           "REMARK 610",
                           "REMARK 615",
                           "REMARK 620",
                           "REMARK 630",
                           "REMARK   3",
                           NULL};
   
   if((gPDBXMLForce == FORCEXML_PDB) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == FALSE))
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

/************************************************************************/
/*>static BOOL blSetPDBMLDateField(char *pdbml_date, char *pdb_date)
   -----------------------------------------------------------------
*//**

   \param[out]    *pdbml_date    PDBML date string 'yyyy-mm-dd'
   \param[in]     *pdb_date      PDB date string   'dd-MTH-yy'
   \return                       Success?

   Convert pdb date format to pdbml date format.

-  12.05.15 Original. By: CTP

*/
static BOOL blSetPDBMLDateField(char *pdbml_date, char *pdb_date)
{
   char month_letter[12][4] = {"JAN","FEB","MAR","APR","MAY","JUN",
                               "JUL","AUG","SEP","OCT","NOV","DEC"};
   int day   = 0,
       month = 0,
       year  = 0,
       items = 0,
       i     = 0;

   char pdb_month[4] = "";
   
   /* parse pdb date */
   items = sscanf(pdb_date, "%2d-%3s-%2d", &day, pdb_month, &year);

   /* convert month_pdb */
   for(i=0; i<12; i++)
   {
      if(!strcmp(month_letter[i], pdb_month))
      {
         month = i + 1;
         break;
      }
   }

   /* convert 2-digit year to 4-digit year (74 == 2074, 75 == 1975)     */
   /* earliest searchable pdb structures date from 1976                 */
   if(year < 75)
   { year = year + 2000; }
   else
   { year = year + 1900; }
   
   /* error check                                                       */
   if(items != 3 || 
      year == 0 || month == 0 || day == 0 || 
      day   < 1 || day > 31   ||
      month < 1 || month > 12 ||
      year  < 1975)
   {
      /* conversion failed                                              */
      pdbml_date[0] = '\0';
      return FALSE;
   }
   
   /* set pdbml date                                                    */
   sprintf(pdbml_date, "%4d-%02d-%02d", year, month, day);

   return TRUE;
}

/************************************************************************/
/*>static HASHTABLE *blMapChainsToEntity(WHOLEPDB *wpdb)
   -----------------------------------------------------
*//**

   \param[in]     *wpdb    WHOLEPDB to parse.
   \return                 HASHTABLE mapping chain labels to MOL_ID

   Map chain labels to MOL_IDs from header COMPND records.
   
   Used to set entity_id for PDBML-format files if entity_id has not been
   set (eg if original input file was PDB-format).

-  18.06.14 Original. By: CTP

*/
static HASHTABLE *blMapChainsToEntity(WHOLEPDB *wpdb)
{
   HASHTABLE *hashtable = NULL;
   char      **chains   = NULL;
   int       nchains    =    0,
             i          =    0;
   COMPND    compnd;

   /* return if wpdb not set                                            */
   if(wpdb == NULL || wpdb->pdb == NULL || wpdb->header == NULL)
   { return NULL; }

   /* get chain array                                                   */
   chains = blGetPDBChainLabels(wpdb->pdb, &nchains);

   /* init hashtable                                                    */
   if((hashtable = blInitializeHash(0))==NULL)
   {
      fprintf(stderr,"No memory for hash table\n");
      return(NULL);
   }

   /* set hashtable                                                     */
   for(i=0;i<nchains;i++)
   {
      if(blGetCompoundWholePDBChain(wpdb, chains[i], &compnd))
      {
         if(compnd.molid)
            blSetHashValueInt(hashtable, chains[i], compnd.molid);
      }
   }

   /* free chains array                                                 */ 
   if(nchains)
   {
      for(i=0;i<nchains;i++){ free(chains[i]); }
      free(chains);
   }

   /* return hashtable                                                  */
   return hashtable;
}

/************************************************************************/
/*>static char **blReadSeqresChainLabelWholePDB(WHOLEPDB *wpdb,
                                                int *nchains)
   ------------------------------------------------------------
*//**

   \param[in]     wpdb      WHOLEPDB structure
   \param[out]    *nchains  Number of chains found
   \return                  Array of chain labels

   Reads the sequence from the SEQRES records from the PDB header
   stored in a WHOLEPDB structure. Creates an array of malloc()'d
   character arrays in which the chain label is stored.

-  09.07.15 Original based on blReadSeqresResidueWholePDB()   By: CTP
*/
static char **ReadSeqresChainLabelWholePDB(WHOLEPDB *wpdb, int *nchains)
{
   STRINGLIST *seqres = NULL, 
              *s;
   char       currchain[2] = " ",
              chain[2]     = " ",
              **chainid;
   int        chainnum = 0,
              i;

   *nchains = 0;
   
   /* First read the SEQRES records into a linked list                  */
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

   /* Return if no seqres records found                                 */
   if(seqres == NULL)
   {
      return(NULL);
   }

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
      chain labels
   */
   if((chainid=(char **)malloc((*nchains) * sizeof(char *)))==NULL)
   {
      FREELIST(seqres, STRINGLIST);
      return(NULL);
   }

   /* allocate memory for chain labels                                  */
   for(i=0;i<(*nchains);i++)
   {
      if((chainid[i] = (char *)malloc(8*sizeof(char))) == NULL)
      {
         FREELIST(seqres, STRINGLIST);
         free(chainid);
         return(NULL);
      }
   }

   /* SECOND PASS: Store the Chain labels                               */
   chainnum  = 0;
   strcpy(currchain,"");
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s",chain);
      if(!CHAINMATCH(chain,currchain))
      {
         /* store new chain id */
         strcpy(chainid[chainnum],chain);
         strcpy(currchain,chain);
         chainnum++;
      }
   }

   /* free memory and return                                            */
   FREELIST(seqres, STRINGLIST);
   return(chainid);
}



/************************************************************************/
/*>static STRINGLIST **ReadSeqresResidueListWholePDB(WHOLEPDB *wpdb,
                                                     int *nchains)
   -----------------------------------------------------------------
*//**

   \param[in]     wpdb      WHOLEPDB structure
   \param[out]    *nchains  Number of chains found
   \return                  Array of sequence strings

   Reads the sequence from the SEQRES records from the PDB header
   stored in a WHOLEPDB structure. Creates an array of malloc()'d
   STRINGLISTs in which the sequence is stored. Can therefore
   cope with any size of sequence information from the PDB file.

-  09.07.15 Original based on blReadSeqresResidueWholePDB()   By: CTP
*/
static STRINGLIST **ReadSeqresResidueListWholePDB(WHOLEPDB *wpdb, 
                                                  int *nchains)
{
   STRINGLIST *seqres = NULL, 
              *s,
              **residuelist = NULL;
   char       currchain[2] = " ",
              chain[2]     = " ",
              res[13][8];
   int        chainnum = 0,
              i;

   *nchains = 0;
   

   /* First read the SEQRES records into a linked list                  */
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

   /* Return if no seqres records found                                 */
   if(seqres == NULL)
   {
      return(NULL);
   }

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

   /* Allocate array for residue sequences                              */
   if((residuelist=
       (STRINGLIST **)malloc((*nchains) * sizeof(STRINGLIST *)))==NULL)
   {
      FREELIST(seqres, STRINGLIST);
      return(NULL);
   }

   /* set stringlist pointers to NULL                                   */
   for(i=0;i<(*nchains);i++)
   {
      residuelist[i] = NULL;
   }

   /* SECOND PASS: Store the sequence                                   */
   chainnum  = 0;
   strncpy(currchain,&(seqres->string[11]),1);
   for(s=seqres; s!=NULL; NEXT(s))
   {
      fsscanf(s->string,"%11x%1s%7x%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s",
              chain,res[0],res[1],res[2],res[3],res[4],res[5],res[6],
              res[7],res[8],res[9],res[10],res[11],res[12]);
      if(!CHAINMATCH(chain,currchain))
      {
         /* Start of new chain                                          */
         strcpy(currchain,chain);
         chainnum++;
      }
      
      /* Store these sequence data                                      */
      for(i=0; i<13; i++)
      {
         /* Break out if not all positions were filled in               */
         if(!strncmp(res[i],"    ",4))
            break;

         /* add to stringlist                                           */
         residuelist[chainnum] = blStoreString(residuelist[chainnum],
                                               res[i]);
      }
   }

   /* free memory and return                                            */
   FREELIST(seqres, STRINGLIST);
   return(residuelist);
}

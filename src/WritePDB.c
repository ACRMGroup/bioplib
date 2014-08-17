/************************************************************************/
/**

   \file       WritePDB.c
   
   \version    V1.15
   \date       16.08.14
   \brief      Write a PDB file from a linked list
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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
   structures of type pdb_entry. The strucure is set up by including the 
   file "pdb.h". For details of the structure, see this file.

**************************************************************************

   Usage:
   ======
   WritePDB(fp, pdb)
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

*************************************************************************/
/* Defines required for includes
*/
#define WRITEPDB_MAIN

/************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <libxml/tree.h>
#include <ctype.h>

#include "MathType.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/*>BOOL blWritePDB(FILE *fp, PDB *pdb)
   -----------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list...

-  21.06.14 Original By: CTP
*/
BOOL blWritePDB(FILE *fp,
                PDB  *pdb)
{
   if((gPDBXMLForce == FORCEXML_XML) ||
      (gPDBXMLForce == FORCEXML_NOFORCE && gPDBXML == TRUE))
   {
      /* Write PDBML file */
      blWriteAsPDBML(fp, pdb);
   }
   else
   {
      /* Check format */
      if(blFormatCheckWritePDB(pdb) == FALSE)
      {
         return FALSE;
      }

      /* Write whole PDB File */
      blWriteAsPDB(fp, pdb);
   }
   
   return TRUE;
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
/*>void blWriteAsPDB(FILE *fp, PDB *pdb)
   -------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list by calls to WritePDBRecord()

-  08.03.89 Original
-  01.06.92 ANSIed and autodoc'd
-  10.06.93 Uses NEXT macro; void type
-  08.07.93 Added insertion of TER cards
-  22.02.94 And a TER card at the end of the file
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  17.06.14 Renamed to blWriteAsPDB() By: CTP
-  07.07.14 Use blWritePDBRecord() By: CTP
*/
void blWriteAsPDB(FILE *fp,
                  PDB  *pdb)
{
   PDB   *p;
   char  PrevChain[8];
   
   strcpy(PrevChain,pdb->chain);

   for(p = pdb ; p ; NEXT(p))
   {
      if(!CHAINMATCH(PrevChain,p->chain))
      {
         /* Chain change, insert TER card                               */
         fprintf(fp,"TER   \n");
         strcpy(PrevChain,p->chain);
      }
      blWritePDBRecord(fp,p);
   }
   fprintf(fp,"TER   \n");
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

   fprintf(fp,"%-6s%5d %-4s%c%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\
          %2s%c%c\n",
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

   Write a PDB record

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
*/
void blWritePDBRecordAtnam(FILE *fp,
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
/*>void blWriteAsPDBML(FILE *fp, PDB *pdb)
   ---------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list in PDBML format.

-  02.06.14 Original. By: CTP
-  21.06.14 Renamed blWriteAsPDBML() and updated symbol handling. By: CTP
-  17.07.14 Use blSetElementSymbolFromAtomName() By: CTP
-  16.08.14 Use element and charge data. By: CTP

*/
void blWriteAsPDBML(FILE *fp, PDB  *pdb)
{
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
   
   /* Create doc */
   doc = xmlNewDoc((xmlChar *) "1.0");
   doc->encoding = xmlStrdup((xmlChar *) "UTF-8");
   
   /* Root node */
   root_node = xmlNewNode(NULL, (xmlChar *) "datablock");
   xmlDocSetRootElement(doc, root_node);
   pdbx = xmlNewNs(root_node, (xmlChar *) "null", (xmlChar *) "PDBx");
   xsi  = xmlNewNs(root_node, (xmlChar *) "null", (xmlChar *) "xsi");
   xmlSetNs(root_node,pdbx);
   
   
   /* Atom_sites node */
   sites_node = xmlNewChild(root_node, NULL,
                            (xmlChar *) "atom_siteCategory", NULL);
   
   /* Atom nodes */
   for(p = pdb ; p ; NEXT(p))
   {
      /* skip TER */
      if(!strncmp("TER",p->resnam,3))
      {
         continue;
      }

      /* Add atom node */
      atom_node = xmlNewChild(sites_node, NULL,
                              (xmlChar *) "atom_site", NULL);
      sprintf(buffer, "%d", p->atnum);
      xmlNewProp(atom_node, (xmlChar *) "id", (xmlChar *) buffer);
      
      /* Add atom data nodes */
      /* B value */
      sprintf(buffer,"%.2f", p->bval);
      node = xmlNewChild(atom_node, NULL, 
                         (xmlChar *) "B_iso_or_equiv",
                         (xmlChar *) buffer);

      /* coordinates */
      sprintf(buffer,"%.3f", p->x);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "Cartn_x",
                         (xmlChar *) buffer);

      sprintf(buffer,"%.3f", p->y);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "Cartn_y",
                         (xmlChar *) buffer);

      sprintf(buffer,"%.3f", p->z);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "Cartn_z",
                         (xmlChar *) buffer);

      /* author atom site labels */
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "auth_asym_id",
                         (xmlChar *) p->chain);

      strcpy(buffer,p->atnam);
      KILLTRAILSPACES(buffer);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "auth_atom_id",
                         (xmlChar *) buffer);

      strcpy(buffer,p->resnam);
      KILLTRAILSPACES(buffer);
      KILLLEADSPACES(buffer_ptr,buffer);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "auth_comp_id",
                         (xmlChar *) buffer_ptr);
      
      sprintf(buffer,"%d", p->resnum);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "auth_seq_id",
                         (xmlChar *) buffer);

      /* record type atom/hetatm */
      strcpy(buffer,p->record_type);
      KILLTRAILSPACES(buffer);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "group_PDB",
                         (xmlChar *) buffer);

      /* atom site labels */
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "label_alt_id",
                         NULL);
      if(p->altpos == ' ')
      {
         xmlNewNsProp(node, xsi, (xmlChar *) "nil", (xmlChar *) "true");
      }
      else
      {
         buffer[0] = p->altpos;
         buffer[1] = '\0';
         xmlNodeSetContent(node, (xmlChar *) buffer);
      }
      
      node = xmlNewChild(atom_node, NULL, 
                         (xmlChar *) "label_asym_id",
                         (xmlChar *) p->chain);

      strcpy(buffer,p->atnam);
      KILLTRAILSPACES(buffer);
      node = xmlNewChild(atom_node, NULL, 
                         (xmlChar *) "label_atom_id",
                         (xmlChar *) buffer);

      strcpy(buffer,p->resnam);
      KILLTRAILSPACES(buffer);
      KILLLEADSPACES(buffer_ptr,buffer);
      node = xmlNewChild(atom_node, NULL, 
                         (xmlChar *) "label_comp_id",
                         (xmlChar *) buffer_ptr);

      /* Note: Entity ID is not stored in PDB data structure. 
               Value set to 1 */
      node = xmlNewChild(atom_node, NULL,
                         (xmlChar *) "label_entity_id",
                         (xmlChar *) "1");
      
      sprintf(buffer,"%d", p->resnum);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "label_seq_id",
                         (xmlChar *) buffer);

      /* occupancy */
      sprintf(buffer,"%.2f", p->occ);
      node = xmlNewChild(atom_node, NULL, (xmlChar *) "occupancy",
                         (xmlChar *) buffer);
                         
      /* insertion code */
      /* Note: Insertion code node only included for residues with 
               insertion codes */
      if(strcmp(p->insert," "))
      {
         sprintf(buffer,"%s", p->insert);
         node = xmlNewChild(atom_node, NULL, 
                            (xmlChar *) "pdbx_PDB_ins_code",
                            (xmlChar *) buffer);
      }

      /* model number */
      /* Note: Model number is not stored in PDB data structure.
               Value set to 1 */
      node = xmlNewChild(atom_node, NULL,
                         (xmlChar *) "pdbx_PDB_model_num",
                         (xmlChar *) "1");

      /* formal charge */
      /* Note: Formal charge node not included for neutral atoms */
      if(p->formal_charge != 0)
      {
         sprintf(buffer,"%d", p->formal_charge);
         node = xmlNewChild(atom_node, NULL, 
                            (xmlChar *) "pdbx_formal_charge",
                            (xmlChar *) buffer);
      }

      /* atom symbol */
      /* Note: If the atomic symbol is not set in PDB data structure then
               the value set is based on columns 13-14 of pdb-formated
               text file.  */
      sprintf(buffer,"%s", p->element);
      KILLLEADSPACES(buffer_ptr,buffer);
      if(strlen(buffer_ptr))
      {
         node = xmlNewChild(atom_node, NULL, (xmlChar *) "type_symbol",
                            (xmlChar *) buffer_ptr);
      }
      else
      {
         blSetElementSymbolFromAtomName(buffer,p->atnam_raw);
         node = xmlNewChild(atom_node, NULL, (xmlChar *) "type_symbol",
                            (xmlChar *) buffer);
      }
   }

   /* Write to doc file pointer */
   xmlDocFormatDump(fp,doc,1);

   /* Free Memory */
    xmlFreeDoc(doc);
    xmlCleanupParser();

   return;
}

/************************************************************************/
/*>void blSetElementSymbolFromAtomName(char *element, char * atom_name)
   --------------------------------------------------------------------
*//**

   \param[out]    *element  Element symbol
   \param[in]     *atom     Atom name

   Set the element symbol (columns 77-78 of a pdb file) based on the 
   atom name (columns 13-16 of a pdb file). 

   The atom name is stored in the PDB data stucture as atnam_raw.

-  17.07.14 Original. By: CTP

*/
void blSetElementSymbolFromAtomName(char *element, char * atom_name)
{
   char  buffer[] = "  ",
         *buffer_ptr;

   /* copy first two letters of atom name to buffer */
   strncpy(buffer,atom_name,2);
   buffer[2] = '\0';
   KILLLEADSPACES(buffer_ptr,buffer);

   /* remove digits */
   if(strlen(buffer_ptr) == 2 && isdigit(buffer[1]))
   {
      buffer[1] = '\0';
   }
   else if(strlen(buffer_ptr) == 2 && !isalpha(buffer[0]))
   {
      buffer_ptr += 1;
   }

   /* fix hydrogens and carbons */
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

   /* copy atom symbol to element */
   strcpy(element, buffer_ptr);
}

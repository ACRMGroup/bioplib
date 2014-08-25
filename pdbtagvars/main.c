/************************************************************************/
/**

   \file       main.c
   
   \version    V0.2
   \date       25.08.14
   \brief      Demonstration of tag pinting code
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1988-2014
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
-  V0.1  06.08.14 Preliminary code
-  V0.2  25.08.14 Added pdbml write function. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

#include "pdbtagvars.h"


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int testInt(PDB *p);
char *testString_GetResID(PDB *p);


/************************************************************************/
/* Globals
*/

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *fp;
   
   blPDBAddXMLAccessTag();
   
   INIT_PDBTAGVAR(&testInt,             PDBTAGVAR_INT,    "test_int");
   INIT_PDBTAGVAR(&testString_GetResID, PDBTAGVAR_STRING, "pdbx_resid");

   if((fp=fopen("test.pdb", "r"))!=NULL)
   {
      PDB *pdb = NULL;
      int natoms;
      
      if((pdb = blReadPDB(fp, &natoms))!=NULL)
      {
         /* 
         blPrintAllTagVariables(pdb);
         */

         /* write pdbml file to stdout */
         testWriteAsPDBML(stdout,pdb);
         
         FREELIST(pdb, PDB);
      }
      else
      {
         fprintf(stderr,"Can't read coordinates from PDB file: test.pdb\n");
      }
      
      fclose(fp);
   }
   else
   {
      fprintf(stderr,"Can't open PDB file: test.pdb\n");
   }
   
   return(0);
}

/************************************************************************/
int testInt(PDB *p)
{
   return(9999);
}
/************************************************************************/
char *testString_GetResID(PDB *p)
{
   static char string[180];
   
   sprintf(string,"%s%d%-s", p->chain, p->resnum, p->insert);
   return(string);
}

/************************************************************************/
/*>void testWriteAsPDBML(FILE *fp, PDB *pdb)
   -----------------------------------------
*//**

   \param[in]     *fp   PDB file pointer to be written
   \param[in]     *pdb  PDB linked list to write

   Write a PDB linked list in PDBML format.
   
   This test function is based on the bioplib function blWriteAsPDBML(). 
   The function calls blAddTagVariablesNodes() which writes additional 
   user-defined tags for each atom.

   Tags are written if gPDBTagWrite is TRUE.

-  25.08.14 Original. By: CTP

*/
void testWriteAsPDBML(FILE *fp, PDB  *pdb)
{
   /* PDBML format supported */
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

      /* NEW CODE */
      /* user-defined tags */
      if(gPDBTagWrite && gNPDBTagFunctions)
      {
         blAddTagVariablesNodes(p,atom_node);
      }
   }

   /* Write to doc file pointer */
   xmlDocFormatDump(fp,doc,1);

   /* Free Memory */
    xmlFreeDoc(doc);
    xmlCleanupParser();

   return;
}



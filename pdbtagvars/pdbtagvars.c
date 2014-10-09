/************************************************************************/
/**

   \file       
   
   \version    V0.3
   \date       28.08.14
   \brief      Code for associating XML tags with additional 
               variables in the PDB structure
   
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
-  V0.2  25.08.14 Added blAddTagVariablesNodes() By: CTP
-  V0.3  28.08.14 Added blAddTagVariablesColumns() and updated 
                  blAddTagVariablesNodes() By: CTP

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

#define _PDBTAGVARS_CODE 1
#include "pdbtagvars.h"

/************************************************************************/
/* Prototypes
*/
void blPDBAddXMLAccessTag(void);

/************************************************************************/
/*>REAL blXMLGetPDBAccess(PDB *pdb)
   --------------------------------
*//**
   \param[in]   PDB  *  Pointer to PDB structure

   Extracts and returns the accessibility from a PDB structure

   This is used by blPDBAddXMLAccessTag()

   THIS ALSO SERVES AS A DEMONSTRATION FOR HOW TO EXTRACT AND RETURN A
   VALUE

-  06.08.14 Original   By: ACRM
*/
REAL blXMLGetPDBAccess(PDB *p)
{
   return(p->access);
}


/************************************************************************/
/*>void blPDBAddXMLAccessTag(void)
   -------------------------------
*//**
   Adds the blXMLGetPDBAccess() function with the <pdbx_accessibility>
   XML tag

   THIS ALSO SERVES AS A DEMONSTRATION FOR HOW TO ASSOCIATE A FUNCTION
   THAT EXTRACTS A VALUE FROM A PDB STRUCTURE WITH AN XML TAG NAME

-  06.08.14 Original   By: ACRM
*/
void blPDBAddXMLAccessTag(void)
{
   INIT_PDBTAGVAR(&blXMLGetPDBAccess, PDBTAGVAR_REAL, "pdbx_accessibility");
}


/************************************************************************/
/*>void blPrintAllTagVariables(PDB *pdb)
   -------------------------------------
*//**
   \param[in]   PDB  *  PDB linked list

   Runs through each item in a PDB linked list and calls the 
   blPrintTagVariables() routine to generate tags for that atom

   DEMO CODE THAT NEEDS REPLACING

-  06.08.14 Original   By: ACRM
*/
void blPrintAllTagVariables(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      blPrintTagVariables(p);
   }
}


/************************************************************************/
/*>void blPrintTagVariables(PDB *pdb)
   ----------------------------------
*//**
   \param[in]   PDB  *  Pointer to a PDB structure

   Dispatches out to each function that has been stored for printing data
   with an XML tag

   DEMO CODE THAT NEEDS REPLACING

-  06.08.14 Original   By: ACRM
*/
void blPrintTagVariables(PDB *p)
{
   REAL realVal = 0.0;
   int  intVal  = 0;
   char *stringVal = NULL;
   int  i;
   
   for(i=0; i<gNPDBTagFunctions; i++)
   {
      switch(gPDBTagFunctions[i].type)
      {
      case PDBTAGVAR_REAL:
         realVal = (*gPDBTagFunctions[i].realFunction)(p); 
         printf("<%s>%f</%s>\n", gPDBTagFunctions[i].tag, realVal, 
                gPDBTagFunctions[i].tag);
         break;
      case PDBTAGVAR_INT:
         intVal = (*gPDBTagFunctions[i].intFunction)(p);
         printf("<%s>%d</%s>\n", gPDBTagFunctions[i].tag, intVal, 
                gPDBTagFunctions[i].tag);
         break;
      case PDBTAGVAR_STRING:
         stringVal = (*gPDBTagFunctions[i].stringFunction)(p);
         printf("<%s>%s</%s>\n", gPDBTagFunctions[i].tag, stringVal, 
                gPDBTagFunctions[i].tag);
         break;
      }
   }
}

   
/************************************************************************/
/*>void blAddTagVariablesNodes(PDB *pdb, xmlNodePtr atom_node)
   -----------------------------------------------------------
*//**
   \param[in]   PDB        *  Pointer to a PDB structure
   \param[in]   atom_node  *  Pointer to an atom_site node

   Adds child nodes with user-defined data to pdbml atom_site node.

   Based on blPrintTagVariables().
   
-  25.08.14 Original   By: CTP
-  28.08.14 Added format for floating point and decimal data.
            Added max length for xmltag_data string. By: CTP
*/
void blAddTagVariablesNodes(PDB *pdb, xmlNodePtr atom_node)
{
   xmlNodePtr node = NULL;
   char       xmltag_name[MAXTAGNAME],
              xmltag_data[MAXTAGDATA];
   int        i;

   /* get functions for pdb list */
   for(i=0; i<gNPDBTagFunctions; i++)
   {
      /* get tag name */
      sprintf(xmltag_name,"%s",gPDBTagFunctions[i].tag);
   
      /* get tag data */
      switch(gPDBTagFunctions[i].type)
      {
      case PDBTAGVAR_REAL:
         sprintf(xmltag_data,"%8.3f",
                 (*gPDBTagFunctions[i].realFunction)(pdb));
         break;
      case PDBTAGVAR_INT:
         sprintf(xmltag_data,"%6d",
                 (*gPDBTagFunctions[i].intFunction)(pdb));
         break;
      case PDBTAGVAR_STRING:
         sprintf(xmltag_data,"%s",
                 (*gPDBTagFunctions[i].stringFunction)(pdb));
         break;
      }

      /* add child node */
      node = xmlNewChild(atom_node, NULL, (xmlChar *) xmltag_name, 
                         (xmlChar *) xmltag_data);
   }

   return;
}

/************************************************************************/
/*>char *blAddTagVariablesCols(PDB *pdb)
   -------------------------------------
*//**
   \param[in]   PDB  *  Pointer to a PDB structure

   Adds columns with user-defined data to pdb-format output.
   
   Floating point data is formatted as "%8.3" and decimal data is 
   formatted as "%d6".
   
   String data does not have a set column width giving the user the 
   freedom to set the data format from the tag function.

-  28.08.14 Original   By: CTP
*/
char *blAddTagVariablesCols(PDB *pdb)
{
   char *xmltag_string = NULL;
   char xmltag_data[MAXTAGDATA];
   int  i;

   /* get function values for pdb atom */
   for(i=0; i<gNPDBTagFunctions; i++)
   {   
      /* get tag data */
      switch(gPDBTagFunctions[i].type)
      {
      case PDBTAGVAR_REAL:
         sprintf(xmltag_data,"%8.3f",
                 (*gPDBTagFunctions[i].realFunction)(pdb));
         break;
      case PDBTAGVAR_INT:
         sprintf(xmltag_data,"%6d",
                 (*gPDBTagFunctions[i].intFunction)(pdb));
         break;
      case PDBTAGVAR_STRING:
         sprintf(xmltag_data,"%s",
                 (*gPDBTagFunctions[i].stringFunction)(pdb));
         break;
      }

      /* update xmltag_string with tag data */
      if(!xmltag_string)
      {
         /* first tag: allocate memory and copy data */
         if( (xmltag_string = (char *)malloc((strlen(xmltag_data) + 1)
                                             * sizeof(char))) == NULL )
         {
            return(NULL);
         }
         strcpy(xmltag_string, xmltag_data);
      }
      else
      {
         /* next tag: realloc memory and append data */
         if( (xmltag_string = realloc(xmltag_string,
                                      (strlen(xmltag_string)
                                       + strlen(xmltag_data) + 1)
                                       * sizeof(char))) == NULL )
         {
            return(NULL);
         }
         strcat(xmltag_string, xmltag_data);
      }
   }

   return xmltag_string;
}


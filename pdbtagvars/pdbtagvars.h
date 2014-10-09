/************************************************************************/
/**

   \file       pdbtagvars.h
   
   \version    V0.3
   \date       28.08.14
   \brief      Header file for associating XML tags with additional 
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
-  V0.3  28.08.14 Added blAddTagVariablesColumns() By: CTP

*************************************************************************/
/* Includes
*/
#include <libxml/tree.h>

/* Defines and macros
*/
#define MAXTAGNAME 360  /* Maximum length of an XML tag name            */
#define MAXTAGDATA 360  /* Maximum length of an XML tag data string     */

typedef struct 
{
   REAL (*realFunction)(PDB *);
   int (*intFunction)(PDB *);
   char *(*stringFunction)(PDB *);
   char tag[MAXTAGNAME];
   int  type;
}  PDBTAGVAR;

#define PDBTAGVAR_REAL   0
#define PDBTAGVAR_INT    1
#define PDBTAGVAR_STRING 2

/************************************************************************/
/* Macro to initialize the binding between a tag name and a function that
   can extract a variable from a PDB structure. The function takes a PDB
   pointer (PDB *) as input and returns REAL, int or (char *) as specified
   by the type (PDBTAGVAR_REAL, PDBTAGVAR_INT or PDBTAGVAR_STRING). The
   taglabel is a string specifying an XML tag that should be associated 
   with the value

   Called as:
   INIT_PDBTAGVAR(&myFunction, type, taglabel)

   06.08.14  Original   By: ACRM
*/
#define INIT_PDBTAGVAR(f, t, l)                                                  \
do {                                                                             \
   if(gPDBTagFunctions == NULL)                                                  \
   {                                                                             \
      gPDBTagFunctions=(PDBTAGVAR *)malloc(sizeof(PDBTAGVAR));                   \
   }                                                                             \
   else                                                                          \
   {                                                                             \
      gPDBTagFunctions=(PDBTAGVAR *)realloc(gPDBTagFunctions,                    \
                    (1+gNPDBTagFunctions)*sizeof(PDBTAGVAR));                    \
   }                                                                             \
   if(gPDBTagFunctions==NULL) break;                                             \
   strncpy(gPDBTagFunctions[gNPDBTagFunctions].tag, l, MAXTAGNAME);              \
   gPDBTagFunctions[gNPDBTagFunctions].realFunction = NULL;                      \
   gPDBTagFunctions[gNPDBTagFunctions].intFunction = NULL;                       \
   gPDBTagFunctions[gNPDBTagFunctions].stringFunction = NULL;                    \
   switch(t)                                                                     \
   {                                                                             \
   case PDBTAGVAR_REAL:                                                          \
      gPDBTagFunctions[gNPDBTagFunctions].realFunction = (REAL (*)(PDB *))f;     \
      break;                                                                     \
   case PDBTAGVAR_INT:                                                           \
      gPDBTagFunctions[gNPDBTagFunctions].intFunction = (int (*)(PDB *))f;       \
      break;                                                                     \
   case PDBTAGVAR_STRING:                                                        \
      gPDBTagFunctions[gNPDBTagFunctions].stringFunction = (char * (*)(PDB *))f; \
      break;                                                                     \
   }                                                                             \
   gPDBTagFunctions[gNPDBTagFunctions].type = t;                                 \
   gNPDBTagFunctions++;                                                          \
} while(0)


/************************************************************************/
/* Globals
*/
#ifdef _PDBTAGVARS_CODE
PDBTAGVAR *gPDBTagFunctions = NULL;
int       gNPDBTagFunctions = 0;
#else
extern PDBTAGVAR *gPDBTagFunctions;
extern int       gNPDBTagFunctions;
#endif


/************************************************************************/
/* Prototypes
*/
void blPDBAddXMLAccessTag(void);
void blPrintTagVariables(PDB *p);
void blPrintAllTagVariables(PDB *pdb);

void blAddTagVariablesNodes(PDB *pdb, xmlNodePtr atom_node);
char *blAddTagVariablesCols(PDB *pdb);

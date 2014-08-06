/************************************************************************/
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
/* Globals
*/

/************************************************************************/
void blPrintAllTagVariables(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      blPrintTagVariables(p);
   }
}
/************************************************************************/
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

   

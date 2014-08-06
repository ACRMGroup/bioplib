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

#include "pdbtagvars.h"


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL testReal_GetAccess(PDB *p);
int testInt(PDB *p);
char *testString_GetResID(PDB *p);


/************************************************************************/
/* Globals
*/

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *fp;
   
   INIT_PDBTAGVAR(&testReal_GetAccess,  PDBTAGVAR_REAL,   "pdbx_accessibility");
   INIT_PDBTAGVAR(&testInt,             PDBTAGVAR_INT,    "test_int");
   INIT_PDBTAGVAR(&testString_GetResID, PDBTAGVAR_STRING, "pdbx_resid");

   if((fp=fopen("test.pdb", "r"))!=NULL)
   {
      PDB *pdb = NULL;
      int natoms;
      
      if((pdb = blReadPDB(fp, &natoms))!=NULL)
      {
         blPrintAllTagVariables(pdb);
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
REAL testReal_GetAccess(PDB *p)
{
   return(p->access);
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



/************************************************************************/
/**

   \file       main.c
   
   \version    V0.1
   \date       06.08.14
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



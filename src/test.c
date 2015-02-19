#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   PDB *pdb, *p;
   FILE *in;
   int natoms;
   
   in=fopen("test.pdb", "r");

/*   
   pdb = ReadPDB(in, &natoms);
   WritePDB(stdout, pdb);
*/
   
   wpdb = ReadWholePDB(in);

/*  
   for(p=wpdb->pdb; p!=NULL; NEXT(p))
   {
      int i;
      if(p->nConect)
      {
         blWritePDBRecord(stdout, p);

         for(i=0; i<p->nConect; i++)
         {
            blWritePDBRecord(stdout, p->conect[i]);
         }

         printf("\n\n");
      }
   }
*/ 
   WriteWholePDB(stdout, wpdb);
   fprintf(stdout,"TEST  \n");

   blBuildConectData(wpdb->pdb);
   blWriteWholePDBTrailer(stdout, wpdb);



   return(0);
}

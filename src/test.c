#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   PDB *pdb, *p;
   FILE *in;
   int natoms;
   DISULPHIDE *dis = NULL, *d;
   BOOL       error;
   
   in=fopen("test.pdb", "r");

   wpdb = blReadWholePDB(in);

   dis = blReadDisulphidesWholePDB(wpdb, &error);
   for(d=dis; d!=NULL; NEXT(d))
   {
      printf("%s%d%s : %s%d%s\n", 
             d->chain1, d->res1, d->insert1,
             d->chain2, d->res2, d->insert2);
   }

   return(0);
}

#include <stdio.h>
#include "pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   FILE *in;
/*
   int nsec, nchains;
   DISULPHIDE *dis = NULL, *d;
   SECSTRUC   *sec = NULL, *s;
   BOOL       error;
   char       **sequences;
*/
   PDB        *p1, *p2;
   
   in=fopen("test.pdb", "r");

   wpdb = blReadWholePDB(in);

/*   p1=blFindResidueSpec(wpdb->pdb, "L23"); */
   p1=wpdb->pdb;
   p2=blFindResidueSpec(wpdb->pdb, "L35");
/*   wpdb->pdb = blDeleteAtomRangePDB(wpdb->pdb, p1, p2); */
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);
   NEXT(p2);
   blAddConect(p1,p2);

   p1->conect[3] = NULL;

   blWriteWholePDB(stdout, wpdb);

   return(0);
}

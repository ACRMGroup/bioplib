#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   FILE *in;
   int nsec, nchains;
   DISULPHIDE *dis = NULL, *d;
   SECSTRUC   *sec = NULL, *s;
   BOOL       error;
   char       **sequences;
   
   in=fopen("test.pdb", "r");

   wpdb = blReadWholePDB(in);

   if((dis = blReadDisulphidesWholePDB(wpdb, &error))!=NULL)
   {
      printf("Disulphides\n");
      for(d=dis; d!=NULL; NEXT(d))
      {
         printf("%s%d%s : %s%d%s\n", 
                d->chain1, d->res1, d->insert1,
                d->chain2, d->res2, d->insert2);
      }
   }

   if((sec = blReadSecWholePDB(wpdb, &nsec))!=NULL)
   {
      printf("Secondary Structure\n");
      for(s=sec; s!=NULL; NEXT(s))
      {
         printf("%s%d%s - %s%d%s : %c\n", 
                s->chain1, s->res1, s->insert1,
                s->chain2, s->res2, s->insert2,
                s->type);
      }
   }

   if((sequences = blReadSeqresWholePDB(wpdb, &nchains)) != NULL)
   {
      int i;
      printf("SEQRES\n");
      for(i=0; i<nchains; i++)
      {
         printf("%s\n", sequences[i]);
      }
   }

   return(0);
}

#include <stdio.h>
#include "bioplib/pdb.h"

int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   PDB *pdb;
   FILE *in;
   int natoms;
   
   in=fopen("/acrm/data/pdb/pdb1yqv.ent", "r");
   
   wpdb = ReadWholePDB(in);
   
   WriteWholePDB(stdout, wpdb);

   return(0);
}

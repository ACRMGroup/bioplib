/* Run Tests */

#include <stdlib.h>
#include <stdio.h>
#include <check.h>

/* Test Suites */
#include "parseresspec_suite.h"
#include "inpdbzone_suite.h"
#include "findzone_suite.h"
#include "getpdbchainlabels_suite.h"
#include "readpdbml_suite.h"
#include "writepdbml_suite.h"


int main(int argc, char **argv)
{
   BOOL verbose;
   int number_failed;
   SRunner *sr;

   printf("Bioplib Tests\n");
   
   /* Setup */
   verbose = (argc > 1 && !strncmp(argv[1],"-v",2)) ? TRUE : FALSE;
   sr = srunner_create(NULL);
   
   /* Add test Suites */
   srunner_add_suite(sr, parseresspec_suite());
   srunner_add_suite(sr, inpdbzone_suite());
   srunner_add_suite(sr, findzone_suite());
   srunner_add_suite(sr, getpdbchainlabels_suite());
   srunner_add_suite(sr, readpdbml_suite());
   srunner_add_suite(sr, writepdbml_suite());
                                                  /* add suites here... */


   /* Run tests */
   if(verbose)
   {
      srunner_run_all(sr, CK_VERBOSE);
   }
   else
   {
      srunner_run_all(sr, CK_NORMAL);
   }
   number_failed = srunner_ntests_failed (sr);
   srunner_free(sr);
   
   return((number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
}

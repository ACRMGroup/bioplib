/************************************************************************/
/**

   \file       main.c
   
   \version    V1.0
   \date       05.08.14
   \brief      Run test suites for BiopLib.

   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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

   Test suite runner for BiopLib. 

**************************************************************************

   Usage:
   ======

   Compile bioplib library from the bioplib/src directory with make:


    cd bioplib/src
    make
 
   Compile tests from the bioplib/src/TEST directory with make to generate
   executable file, run_tests.


    cd TEST
    make
    ./run_tests

   For verbose output use the -v option:


    ./run_tests -v
 

**************************************************************************

   Revision History:
   =================

-  V1.0  05.08.14 Original By: CTP

*************************************************************************/

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
#include "wholepdb_suite.h"


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
   srunner_add_suite(sr, wholepdb_suite());
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

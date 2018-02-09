/************************************************************************/
/**

   \file       findzone_suite.c
   
   \version    V1.0
   \date       05.08.14
   \brief      Test suite for blFindZonePDB().
   
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

   Test suite for blFindZonePDB().

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.08.14 Original By: CTP

*************************************************************************/

#include "findzone_suite.h"

/* Defines */
#define TEST_PDB_FILE "./data/test-deca-ala-01.pdb"

/* Globals */
PDB *pdb_in    = NULL,
    *pdb_start = NULL,
    *pdb_stop  = NULL;

int  start;
char startinsert[8];
int  stop;
char stopinsert[8];
char chain[8];
int  mode;

BOOL output, expected_output;
int expected_start_atom, expected_stop_atom;

/* Setup And Teardown */
void findzone_setup(void)
{
   FILE *fp;
   int natom = 0;
   
   fp = fopen(TEST_PDB_FILE,"r");

   if(fp == NULL)
   {
      fprintf(stderr, "Failed to open test pdb file!\n");
      return;
   }
   
   pdb_in = blReadPDB(fp,&natom);
   fclose(fp);
   
   if(pdb_in == NULL)
   {
      fprintf(stderr, "Failed to read test pdb file!\n");
   }
}

void findzone_teardown(void)
{
   /* Free PDB */
   FREELIST(pdb_in,PDB);
}

/* PDB Data Read Test */
START_TEST(test_read_01)
{
   ck_assert_msg(pdb_in != NULL, "No data read from test file.");
}
END_TEST


/* Core tests */
START_TEST(test_01)
{
   /* Input */
   start               =  2  ;
   strcpy(startinsert,   " ");
   stop                =  5  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "A");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    6;
   expected_stop_atom  =   26;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
}
END_TEST

START_TEST(test_02)
{
   /* Input */
   start               =  2  ;
   strcpy(startinsert,   " ");
   stop                =  3  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "B");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =   36;
   expected_stop_atom  =   46;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
}
END_TEST

START_TEST(test_03)
{
   /* Input */
   start                =  2  ;
   strcpy(startinsert,    " ");
   stop                 =  5  ;
   strcpy(stopinsert,     " ");
   strcpy(chain,          " ");
   mode = ZONE_MODE_SEQUENTIAL;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    6;
   expected_stop_atom  =   26;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
}
END_TEST

START_TEST(test_04)
{
   /* Input */
   start                =  2  ;
   strcpy(startinsert,    " ");
   stop                 =  9  ;
   strcpy(stopinsert,     " ");
   strcpy(chain,          " ");
   mode = ZONE_MODE_SEQUENTIAL;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    6;
   expected_stop_atom  =   46;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
}
END_TEST

/* limits tests */
START_TEST(test_limit_01)
{
   /* Input */
   start               = -999  ;
   strcpy(startinsert,   " ");
   stop                = -999  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         " ");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    1;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         == NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );

}
END_TEST

START_TEST(test_limit_02)
{
   /* Input */
   start               = -999  ;
   strcpy(startinsert,   " ");
   stop                = -999  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "A");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    1;
   expected_stop_atom  =   31;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
}
END_TEST

START_TEST(test_limit_03)
{
   /* Input */
   start             = -999  ;
   strcpy(startinsert,   " ");
   stop              = -999  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "B");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =   31;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         == NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );

}
END_TEST

START_TEST(test_limit_04)
{
   /* Input */
   start             = -999  ;
   strcpy(startinsert,   " ");
   stop              = -999  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,   "nochain");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = FALSE;
   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        == NULL                );
   ck_assert(pdb_stop         == NULL                );
}
END_TEST

START_TEST(test_limit_05)
{
   /* Input */
   start             = -999  ;
   strcpy(startinsert,   " ");
   stop              =    5  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "A");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    1;
   expected_stop_atom  =   26;
   

   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
   
}
END_TEST

START_TEST(test_limit_06)
{
   /* Input */
   start             =    1  ;
   strcpy(startinsert,   " ");
   stop              =  -999 ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "A");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    1;
   expected_stop_atom  =   31;
   

   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
   
}
END_TEST

START_TEST(test_limit_07)
{
   /* Input */
   start             = -999  ;
   strcpy(startinsert,   " ");
   stop              =    3  ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "B");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =    1;
   expected_stop_atom  =   46;
   

   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         != NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   ck_assert(pdb_stop->atnum  == expected_stop_atom  );
   
}
END_TEST

START_TEST(test_limit_08)
{
   /* Input */
   start             =    1  ;
   strcpy(startinsert,   " ");
   stop              =  -999 ;
   strcpy(stopinsert,    " ");
   strcpy(chain,         "B");
   mode =    ZONE_MODE_RESNUM;
   
   /* Expected Output */
   expected_output     = TRUE;
   expected_start_atom =   31;
   

   
   /* Run Test */
   ck_assert_msg(pdb_in != NULL,"No pdb data found.");
   output =  blFindZonePDB(pdb_in, start, startinsert, stop, stopinsert, 
                           chain, mode, &pdb_start, &pdb_stop);
                           
   ck_assert(output           == expected_output     );
   ck_assert(pdb_start        != NULL                );
   ck_assert(pdb_stop         == NULL                );
   ck_assert(pdb_start->atnum == expected_start_atom );
   
}
END_TEST



/* Create Suite */
Suite *findzone_suite(void)
{
   Suite *s        = suite_create("FindZonePDB");
   TCase *tc_read  = tcase_create("Read");
   TCase *tc_core  = tcase_create("Core");   
   TCase *tc_limit = tcase_create("Limits");


   /* Check read from test file */
   tcase_add_checked_fixture(tc_read, findzone_setup, findzone_teardown);
   tcase_add_test(tc_read, test_read_01);
   suite_add_tcase(s, tc_read);   
   
   /* Core test case */
   tcase_add_checked_fixture(tc_core, findzone_setup, findzone_teardown);
   tcase_add_test(tc_core, test_01);
   tcase_add_test(tc_core, test_02);
   tcase_add_test(tc_core, test_03);
   tcase_add_test(tc_core, test_04);
   suite_add_tcase(s, tc_core);

   /* Limits test case */
   tcase_add_checked_fixture(tc_limit, findzone_setup, findzone_teardown);
   tcase_add_test(tc_limit, test_limit_01);
   tcase_add_test(tc_limit, test_limit_02);
   tcase_add_test(tc_limit, test_limit_03);
   tcase_add_test(tc_limit, test_limit_04);
   tcase_add_test(tc_limit, test_limit_05);
   tcase_add_test(tc_limit, test_limit_06);
   tcase_add_test(tc_limit, test_limit_07);
   tcase_add_test(tc_limit, test_limit_08);
   suite_add_tcase(s, tc_limit);


   return(s);
}

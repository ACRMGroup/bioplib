#include "inpdbzone_suite.h"

/* Globals */
PDB  *pdb;
char chain[8],
     insert1[8], insert2[8],
     resspec1[8], resspec2[8],
     single_letter_chain,
     single_letter_insert1,
     single_letter_insert2;

int  resnum1, resnum2;
BOOL output, expected_output;


/* Setup And Teardown */
void inpdbzone_setup(void)
{
   /* Input PDB */
   pdb = (PDB *) malloc(sizeof(PDB));
   strcpy(pdb->record_type, "ATOM");
   pdb->atnum =                  2 ;
   strcpy(pdb->atnam,       " CA ");
   strcpy(pdb->atnam_raw,   " CA ");
   strcpy(pdb->resnam,       "ALA");
   strcpy(pdb->chain,          "A");
   pdb->resnum =                 1 ;
   strcpy(pdb->insert,         " ");
   pdb->x =                  0.000 ;
   pdb->y =                  0.000 ;
   pdb->z =                  0.000 ;
   pdb->occ =                 1.00 ;
   pdb->altpos =               ' ' ;
   pdb->next =                NULL ;
    
   /* Range chain, resnum, insert */
   strcpy(chain,  "A");
   resnum1 =        1 ;
   strcpy(insert1," ");
   resnum2 =        1 ;
   strcpy(insert2," ");
}

void inpdbzone_teardown(void)
{
   /* Free PDB */
   free(pdb);
}

/* Core tests */
START_TEST(test_01)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_02)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =        5 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_03)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       25 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST


START_TEST(test_04)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "B");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_05)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_06)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_07)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "A");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_08)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "A");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_09)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_10)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_11)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output = TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_12)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           20 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

/* range == single resnum */
START_TEST(test_single_resnum_01)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           10 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_02)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           10 ;
   strcpy(insert2,    " ");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_03)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           10 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_04)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    " ");
   resnum2 =           10 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_05)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "A");
   resnum2 =           10 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_06)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "A");
   resnum2 =           10 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_07)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "A");
   resnum2 =           10 ;
   strcpy(insert2,    "A");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"B");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_single_resnum_08)
{
   /* range */
   strcpy(chain,      "A");
   resnum1 =           10 ;
   strcpy(insert1,    "B");
   resnum2 =           10 ;
   strcpy(insert2,    "B");

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output = FALSE;
   
   output = blInPDBZone(pdb, chain, resnum1, insert1, resnum2, insert2);
   ck_assert( output == expected_output );
}
END_TEST

/* Wrapper functions */

/* InPDBZoneSpec() */
START_TEST(test_wrapper_a_01)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = ' ';
   resnum2               =  20;
   single_letter_insert2 = ' ';

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_a_02)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = 'A';
   resnum2               =  20;
   single_letter_insert2 = ' ';

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_a_03)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = 'A';
   resnum2               =  20;
   single_letter_insert2 = ' ';

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       10 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_a_04)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = ' ';
   resnum2               =  20;
   single_letter_insert2 = ' ';

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output = FALSE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_a_05)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = ' ';
   resnum2               =  20;
   single_letter_insert2 = 'A';

   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       20 ;
   strcpy(pdb->insert,"A");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_a_06)
{
   /* range */
   single_letter_chain   = 'A';
   resnum1               =  10;
   single_letter_insert1 = ' ';
   resnum2               =  20;
   single_letter_insert2 = ' ';

   /* pdb */
   strcpy(pdb->chain, "B");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = InPDBZone(pdb, single_letter_chain, 
                      resnum1, single_letter_insert1, 
                      resnum2, single_letter_insert2);
   ck_assert( output == expected_output );
}
END_TEST


/* InPDBZoneSpec() */
START_TEST(test_wrapper_b_01)
{
   /* range */
   strcpy(resspec1, "A10");
   strcpy(resspec2, "A20");
   
   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output =  TRUE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_b_02)
{
   /* range */
   strcpy(resspec1, "A10");
   strcpy(resspec2, "B20");
   
   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = FALSE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_b_03)
{
   /* range */
   strcpy(resspec1,   "*");
   strcpy(resspec2,    "");
   
   /* pdb */
   strcpy(pdb->chain, " ");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = TRUE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_b_04)
{
   /* range */
   strcpy(resspec1,  "A*");
   strcpy(resspec2,    "");
   
   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = TRUE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_b_05)
{
   /* range */
   strcpy(resspec1,  ".*");
   strcpy(resspec2,    "");
   
   /* pdb */
   strcpy(pdb->chain, " ");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = TRUE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST

START_TEST(test_wrapper_b_06)
{
   /* range */
   strcpy(resspec1, "A.*");
   strcpy(resspec2,    "");
   
   /* pdb */
   strcpy(pdb->chain, "A");
   pdb->resnum =       15 ;
   strcpy(pdb->insert," ");
   
   /* expected output */
   expected_output = TRUE;
   
   output = InPDBZoneSpec(pdb, resspec1, resspec2);
   ck_assert( output == expected_output );
}
END_TEST


/* Create Suite */
Suite *inpdbzone_suite(void)
{
   Suite *s = suite_create("InPDBZone");
   
   /* blInPDBZone() */
   TCase *tc_core = tcase_create("Core");
   tcase_add_checked_fixture(tc_core, inpdbzone_setup, 
                             inpdbzone_teardown);
   tcase_add_test(tc_core, test_01);
   tcase_add_test(tc_core, test_02);
   tcase_add_test(tc_core, test_03);
   tcase_add_test(tc_core, test_04);
   tcase_add_test(tc_core, test_05);
   tcase_add_test(tc_core, test_06);
   tcase_add_test(tc_core, test_07);
   tcase_add_test(tc_core, test_08);
   tcase_add_test(tc_core, test_09);
   tcase_add_test(tc_core, test_10);
   tcase_add_test(tc_core, test_11);
   tcase_add_test(tc_core, test_12);
   suite_add_tcase(s, tc_core);

   TCase *tc_single_resnum = tcase_create("Single_Resnum");
   tcase_add_checked_fixture(tc_single_resnum, inpdbzone_setup, 
                             inpdbzone_teardown);
   tcase_add_test(tc_single_resnum, test_single_resnum_01);
   tcase_add_test(tc_single_resnum, test_single_resnum_02);
   tcase_add_test(tc_single_resnum, test_single_resnum_03);
   tcase_add_test(tc_single_resnum, test_single_resnum_04);
   tcase_add_test(tc_single_resnum, test_single_resnum_05);
   tcase_add_test(tc_single_resnum, test_single_resnum_06);
   tcase_add_test(tc_single_resnum, test_single_resnum_07);
   tcase_add_test(tc_single_resnum, test_single_resnum_08);
   suite_add_tcase(s, tc_single_resnum);

   /* InPDBZone() */
   TCase *tc_wrap_a = tcase_create("Wrapper_A");
   tcase_add_checked_fixture(tc_wrap_a, inpdbzone_setup, 
                             inpdbzone_teardown);
   tcase_add_test(tc_wrap_a, test_wrapper_a_01);
   tcase_add_test(tc_wrap_a, test_wrapper_a_02);
   tcase_add_test(tc_wrap_a, test_wrapper_a_03);
   tcase_add_test(tc_wrap_a, test_wrapper_a_04);
   tcase_add_test(tc_wrap_a, test_wrapper_a_05);
   tcase_add_test(tc_wrap_a, test_wrapper_a_06);   
   suite_add_tcase(s, tc_wrap_a);
   
   /* InPDBZoneSpec() */
   TCase *tc_wrap_b = tcase_create("Wrapper_B");
   tcase_add_checked_fixture(tc_wrap_b, inpdbzone_setup, 
                             inpdbzone_teardown);
   tcase_add_test(tc_wrap_b, test_wrapper_b_01);
   tcase_add_test(tc_wrap_b, test_wrapper_b_02);
   tcase_add_test(tc_wrap_b, test_wrapper_b_03);
   tcase_add_test(tc_wrap_b, test_wrapper_b_04);
   tcase_add_test(tc_wrap_b, test_wrapper_b_05);
   tcase_add_test(tc_wrap_b, test_wrapper_b_06);   
   suite_add_tcase(s, tc_wrap_b);
   return(s);
}

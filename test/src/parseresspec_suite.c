/************************************************************************/
/**

   \file       parseresspec.c
   
   \version    V1.0
   \date       05.08.14
   \brief      Test suite for blParseResSpec().
   
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

   Test suite for blParseResSpec().

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.08.14 Original By: CTP

*************************************************************************/

#include "parseresspec_suite.h"

/* Globals */
char chain[8], insert[8], spec[16];
int  resnum;

/* Setup And Teardown */
void parseresspec_setup(void)
{
   /* No setup required */
}

void parseresspec_teardown(void)
{
   /* No teardown required */
}

/* Core tests */
START_TEST(test_std_run)
{
   BOOL output = FALSE;
   strcpy(spec,"A1 ");
   
   output = blParseResSpec(spec, chain, &resnum, insert);
   ck_assert_msg(output,"blParseResSpec() returned FALSE for 'A1 '.");
}
END_TEST

START_TEST(test_std_format)
{
   strcpy(spec,"A1 ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "A");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_dot_run)
{
   BOOL output = FALSE;
   strcpy(spec,"A.1 ");
   
   output = blParseResSpec(spec, chain, &resnum, insert);
   ck_assert_msg(output,"blParseResSpec() returned FALSE for 'A.1 '.");
}
END_TEST

START_TEST(test_dot_format)
{
   strcpy(spec,"A.1 ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "A");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_multi_run)
{
   BOOL output = FALSE;
   strcpy(spec,"Ab1 ");
   
   output = blParseResSpec(spec, chain, &resnum, insert);
   ck_assert_msg(output,"blParseResSpec() returned FALSE for 'Ab1 '.");
}
END_TEST

START_TEST(test_multi_format)
{
   strcpy(spec,"Ab1 ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "Ab");
   ck_assert_int_eq(resnum,  1 );
   ck_assert_str_eq(insert, " ");
}
END_TEST

/* Standard format "A123Y" */
START_TEST(test_std_01)
{
   strcpy(spec,"a1 ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "a");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_std_02)
{
   strcpy(spec," 1 ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, " ");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_std_03)
{
   strcpy(spec," 1");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, " ");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_std_04)
{
   strcpy(spec,"1");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, " ");
   ck_assert_int_eq(resnum, 1 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_std_05)
{
   strcpy(spec,"A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "A");
   ck_assert_int_eq(resnum, 0 );
   ck_assert_str_eq(insert," ");
}
END_TEST

START_TEST(test_std_06)
{
   strcpy(spec,"A123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,   "A");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  "A");
}
END_TEST

START_TEST(test_std_07)
{
   strcpy(spec,"A-123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,    "A");
   ck_assert_int_eq(resnum, -123 );
   ck_assert_str_eq(insert,   "A");
}
END_TEST

START_TEST(test_std_08)
{
   strcpy(spec,"-123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,    " ");
   ck_assert_int_eq(resnum, -123 );
   ck_assert_str_eq(insert,   "A");
}
END_TEST


/* Dot format "A.123A" */
START_TEST(test_dot_01)
{
   strcpy(spec,"A.123");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,   "A");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST

START_TEST(test_dot_02)
{
   strcpy(spec,"A.123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,   "A");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  "A");
}
END_TEST

START_TEST(test_dot_03)
{
   strcpy(spec,".123");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,   " ");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST

START_TEST(test_dot_04)
{
   strcpy(spec,"1.123");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,   "1");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST

START_TEST(test_dot_05)
{
   strcpy(spec,"A.-123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,    "A");
   ck_assert_int_eq(resnum, -123 );
   ck_assert_str_eq(insert,   "A");
}
END_TEST

START_TEST(test_dot_06)
{
   strcpy(spec,".-123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,    " ");
   ck_assert_int_eq(resnum, -123 );
   ck_assert_str_eq(insert,   "A");
}
END_TEST

/* Multi-letter chain "Abc.123Y"*/
START_TEST(test_multi_01)
{
   strcpy(spec,"Abc.123");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "Abc");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST

START_TEST(test_multi_02)
{
   strcpy(spec,"Abc.123A");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "Abc");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  "A");
}
END_TEST

START_TEST(test_multi_03)
{
   strcpy(spec,"Ab1.123");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "Ab1");
   ck_assert_int_eq(resnum, 123 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST

START_TEST(test_multi_04)
{
   strcpy(spec,"Abc");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain, "Abc");
   ck_assert_int_eq(resnum,   0 );
   ck_assert_str_eq(insert,  " ");
}
END_TEST


/* Limits */
START_TEST(test_limits_01)
{
   strcpy(spec,"");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,  " ");
   ck_assert_int_eq(resnum,  0 );
   ck_assert_str_eq(insert, " ");
}
END_TEST

START_TEST(test_limits_02)
{
   strcpy(spec," ");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,  " ");
   ck_assert_int_eq(resnum,  0 );
   ck_assert_str_eq(insert, " ");
}
END_TEST

START_TEST(test_limits_03)
{
   strcpy(spec,".");
   
   ck_assert(blParseResSpec(spec, chain, &resnum, insert));
   
   ck_assert_str_eq(chain,  " ");
   ck_assert_int_eq(resnum,  0 );
   ck_assert_str_eq(insert, " ");
}
END_TEST


/* Create Suite */
Suite *parseresspec_suite(void)
{
   Suite *s         = suite_create("ParseResSpec");
   TCase *tc_core   = tcase_create("Core"),
         *tc_std    = tcase_create("Std"),
         *tc_dot    = tcase_create("Dot"),
         *tc_multi  = tcase_create("Multi"),
         *tc_limits = tcase_create("Limits");

   /* Core test case */
   tcase_add_checked_fixture(tc_core, parseresspec_setup, 
                             parseresspec_teardown);
   tcase_add_test(tc_core, test_std_run);
   tcase_add_test(tc_core, test_std_format);
   tcase_add_test(tc_core, test_dot_run);
   tcase_add_test(tc_core, test_dot_format);
   tcase_add_test(tc_core, test_multi_run);
   tcase_add_test(tc_core, test_multi_format);
   suite_add_tcase(s, tc_core);

   /* Standard format test case */
   tcase_add_test(tc_std, test_std_01);
   tcase_add_test(tc_std, test_std_02);
   tcase_add_test(tc_std, test_std_03);
   tcase_add_test(tc_std, test_std_04);
   tcase_add_test(tc_std, test_std_05);
   tcase_add_test(tc_std, test_std_06);
   tcase_add_test(tc_std, test_std_07);
   tcase_add_test(tc_std, test_std_08);
   suite_add_tcase(s, tc_std);
   
   /* Dot format test case */
   tcase_add_test(tc_dot, test_dot_01);
   tcase_add_test(tc_dot, test_dot_02);
   tcase_add_test(tc_dot, test_dot_03);
   tcase_add_test(tc_dot, test_dot_04);
   tcase_add_test(tc_dot, test_dot_05);
   tcase_add_test(tc_dot, test_dot_06);
   suite_add_tcase(s, tc_dot);

   /* Multi-letter chain test case */
   tcase_add_test(tc_multi, test_multi_01);
   tcase_add_test(tc_multi, test_multi_02);
   tcase_add_test(tc_multi, test_multi_03);
   tcase_add_test(tc_multi, test_multi_04);
   suite_add_tcase(s, tc_multi);

   /* Limits test case */
   tcase_add_test(tc_limits, test_limits_01);
   tcase_add_test(tc_limits, test_limits_02);
   tcase_add_test(tc_limits, test_limits_03);
   suite_add_tcase(s, tc_limits);

   return s;
}

#include "readpdbml_suite.h"

/* Globals */
PDB  *pdb;
FILE  *fp;
int   natoms = 0;

/* Setup And Teardown */
void readpdbml_setup(void)
{
   /* Input PDB */
   pdb = NULL;
}

void readpdbml_teardown(void)
{
   /* Free PDB */
   FREELIST(pdb, PDB);
}

/* Core tests */
START_TEST(test_read_pdb)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon.pdb";
   
   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* tests */
   ck_assert_msg(pdb != NULL, "Failed to read PDB file.");
   ck_assert_int_eq(natoms, 1);
}
END_TEST

START_TEST(test_read_pdbml)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* tests */
   ck_assert_msg(pdb != NULL, "Failed to read PDBML file.");
   ck_assert_int_eq(natoms, 1);
}
END_TEST

/* PDB Tests */
START_TEST(test_read_pdb_data_01)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdb_data_02)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_zinc.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "ZN  ");
   ck_assert_str_eq( pdb->atnam_raw, "ZN  ");
   ck_assert_str_eq( pdb->resnam,    " ZN ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdb_data_03)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_insert.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       "A");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdb_data_04)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_01.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       0.75);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST

START_TEST(test_read_pdb_data_05)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_02.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        4.000);
   ck_assert(        pdb->y ==        5.000);
   ck_assert(        pdb->z ==        6.000);
   ck_assert(        pdb->occ ==       0.75);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST

START_TEST(test_read_pdb_data_06)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_03.pdb";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        4.000);
   ck_assert(        pdb->y ==        5.000);
   ck_assert(        pdb->z ==        6.000);
   ck_assert(        pdb->occ ==       0.50);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST


/* PDBML Tests */
START_TEST(test_read_pdbml_data_01)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdbml_data_02)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_zinc.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "ZN  ");
   ck_assert_str_eq( pdb->atnam_raw, "ZN  ");
   ck_assert_str_eq( pdb->resnam,    " ZN ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdbml_data_03)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_insert.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       "A");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       1.00);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
}
END_TEST

START_TEST(test_read_pdbml_data_04)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_01.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        1.000);
   ck_assert(        pdb->y ==        2.000);
   ck_assert(        pdb->z ==        3.000);
   ck_assert(        pdb->occ ==       0.75);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST

START_TEST(test_read_pdbml_data_05)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_02.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        4.000);
   ck_assert(        pdb->y ==        5.000);
   ck_assert(        pdb->z ==        6.000);
   ck_assert(        pdb->occ ==       0.75);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST

START_TEST(test_read_pdbml_data_06)
{
   /* set filename */
   char filename[] = "data/readpdbml_suite/test_alpha_carbon_alt_03.xml";

   /* read test file */
   fp = fopen(filename,"r");
   pdb = ReadPDB(fp, &natoms);
   fclose(fp);

   /* check data */
   ck_assert(        pdb->atnum ==        1);
   ck_assert_str_eq( pdb->atnam,     "CA  ");
   ck_assert_str_eq( pdb->atnam_raw, " CA ");
   ck_assert_str_eq( pdb->resnam,    "ALA ");
   ck_assert_str_eq( pdb->chain,        "A");
   ck_assert(        pdb->resnum ==       1);
   ck_assert_str_eq( pdb->insert,       " ");
   ck_assert(        pdb->x ==        4.000);
   ck_assert(        pdb->y ==        5.000);
   ck_assert(        pdb->z ==        6.000);
   ck_assert(        pdb->occ ==       0.50);
   ck_assert(        pdb->altpos ==     ' ');
   ck_assert(        pdb->bval ==     20.00);
   ck_assert(        natoms ==            1);
}
END_TEST


/* Create Suite */
Suite *readpdbml_suite(void)
{
   Suite *s = suite_create("ReadPDBML");

   /* Core test case */
   TCase *tc_core = tcase_create("Core");
   tcase_add_checked_fixture(tc_core, readpdbml_setup, 
                             readpdbml_teardown);
   tcase_add_test(tc_core, test_read_pdb);
   tcase_add_test(tc_core, test_read_pdbml);
   suite_add_tcase(s, tc_core);
   
   TCase *tc_pdb = tcase_create("PDB");
   tcase_add_checked_fixture(tc_pdb, readpdbml_setup, 
                             readpdbml_teardown);
   tcase_add_test(tc_pdb, test_read_pdb_data_01);
   tcase_add_test(tc_pdb, test_read_pdb_data_02);
   tcase_add_test(tc_pdb, test_read_pdb_data_03);
   tcase_add_test(tc_pdb, test_read_pdb_data_04);
   tcase_add_test(tc_pdb, test_read_pdb_data_05);
   tcase_add_test(tc_pdb, test_read_pdb_data_06);
   suite_add_tcase(s, tc_pdb);
   
   TCase *tc_pdbml = tcase_create("PDBML");
   tcase_add_checked_fixture(tc_pdbml, readpdbml_setup, 
                             readpdbml_teardown);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_01);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_02);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_03);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_04);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_05);
   tcase_add_test(tc_pdbml, test_read_pdbml_data_06);
   suite_add_tcase(s, tc_pdbml);

   return s;
}
#include "writepdbml_suite.h"

/* Globals */
static char test_output_filename[]     = "/tmp/test-XXXXX",
            test_example_filename[160] = "data/writepdbml_suite/";

static FILE *fp             =  NULL;
static PDB  *pdb_out        =  NULL,
            *pdb_in         =  NULL;
static BOOL files_identical = FALSE;


/* Compare file function */
static BOOL compare_files(char *filename_a, char *filename_b)
{
   char command[120];
   
   /* set command */
   sprintf(command,"cmp %s %s > /dev/null", filename_a, filename_b);
   
   /* return TRUE if files match */
   return system(command) == 0 ? TRUE:FALSE;
}


/* Setup And Teardown */
void writepdbml_setup(void)
{
   /* Output PDB */
   pdb_out = (PDB *) malloc(sizeof(PDB));
   strcpy(pdb_out->record_type, "ATOM");
   pdb_out->atnum =                  1 ;
   strcpy(pdb_out->atnam,       "CA  ");
   strcpy(pdb_out->atnam_raw,   " CA ");
   strcpy(pdb_out->resnam,      "ALA ");
   strcpy(pdb_out->chain,          "A");
   pdb_out->resnum =                 1 ;
   strcpy(pdb_out->insert,         " ");
   pdb_out->x =                  1.000 ;
   pdb_out->y =                  2.000 ;
   pdb_out->z =                  3.000 ;
   pdb_out->occ =                 1.00 ;
   pdb_out->altpos =               ' ' ;
   pdb_out->bval =               20.00 ;
   pdb_out->next =                NULL ;

   /* Set temp file name */
   mkstemp(test_output_filename);
}

void writepdbml_teardown(void)
{
   /* Free PDB */
   FREELIST(pdb_out,PDB);
   FREELIST(pdb_in ,PDB);
}

/* Core tests */
START_TEST(test_write_pdb)
{
   int natoms = 0;

  /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDB(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,         "Failed to open file.");

   /* read test file */
   pdb_in = ReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in != NULL,     "Failed to read file.");
   ck_assert_msg(natoms ==  1,   "No atoms read from file.");
}
END_TEST

START_TEST(test_write_pdbml)
{
   int natoms = 0;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,         "Failed to open file.");

   /* read test file */
   pdb_in = ReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in != NULL,     "Failed to read file.");
   ck_assert_msg(natoms ==  1,   "No atoms read from file.");
}
END_TEST


/* PDB */
START_TEST(test_write_pdb_data_01)
{
   /* set message and example file */
   char test_error_msg[] = "Write alpha carbon.";
   char filename[]       = "test_alpha_carbon.pdb";

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdb_data_02)
{
   /* set message and example file */
   char test_error_msg[] = "Set insert code.";
   char filename[]       = "test_alpha_carbon_insert.pdb";

   /* Update PDB */
   strcpy(pdb_out->insert, "A");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdb_data_03)
{
   /* set message and example file */
   char test_error_msg[] = "Set alt position.";
   char filename[]       = "test_alpha_carbon_alt.pdb";

   /* Update PDB */
   pdb_out->altpos = 'A' ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdb_data_04)
{
   /* set message and example file */
   char test_error_msg[] = "Write zinc atom.";
   char filename[]       = "test_zinc.pdb";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "ZN  ");
   strcpy(pdb_out->atnam_raw,     "ZN  ");
   strcpy(pdb_out->resnam,        " ZN ");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


/* PDBML */
START_TEST(test_write_pdbml_data_01)
{
   /* set message and example file */
   char test_error_msg[] = "Write alpha carbon.";
   char filename[]       = "test_alpha_carbon.xml";

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_02)
{
   /* set message and example file */
   char test_error_msg[] = "Set insert code.";
   char filename[]       = "test_alpha_carbon_insert.xml";

   /* Update PDB */
   strcpy(pdb_out->insert, "A");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_03)
{
   /* set message and example file */
   char test_error_msg[] = "Set alt position.";
   char filename[]       = "test_alpha_carbon_alt.xml";

   /* Update PDB */
   pdb_out->altpos = 'A' ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_04)
{
   /* set message and example file */
   char test_error_msg[] = "Write zinc atom.";
   char filename[]       = "test_zinc.xml";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "ZN  ");
   strcpy(pdb_out->atnam_raw,     "ZN  ");
   strcpy(pdb_out->resnam,        " ZN ");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST



START_TEST(test_write_pdbml_data_05)
{
   /* set message and example file */
   char test_error_msg[] = "Multi-letter chain name.";
   char filename[]       = "test_alpha_carbon_mlchain.xml";

   /* Update PDB */
   strcpy(pdb_out->chain, "Light_A");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


START_TEST(test_write_pdbml_data_06)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 01";
   char filename[]       = "test_H_leu_beta.xml";

   /* Update PDB */
   strcpy(pdb_out->atnam,         "HB2 ");
   strcpy(pdb_out->atnam_raw,     " HB2");
   strcpy(pdb_out->resnam,        "LEU ");
   pdb_out->bval =                 0.00  ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_07)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 02";
   char filename[]       = "test_H_leu_gamma.xml";

   /* Update PDB */
   strcpy(pdb_out->atnam,         "HG  ");
   strcpy(pdb_out->atnam_raw,     " HG ");
   strcpy(pdb_out->resnam,        "LEU ");
   pdb_out->bval =                 0.00  ;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


/* NEED FIX CODE Gives -HD- */
START_TEST(test_write_pdbml_data_08)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 03";
   char filename[]       = "test_H_leu_delta.xml";

   /* Update PDB */
   strcpy(pdb_out->atnam,         "HD11");
   strcpy(pdb_out->atnam_raw,     "HD11");
   strcpy(pdb_out->resnam,        "LEU ");
   pdb_out->bval =                 0.00  ;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_09)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 04";
   char filename[]       = "test_H_dna.xml";

   /* Update PDB */
   strcpy(pdb_out->atnam,         "HO3'");
   strcpy(pdb_out->atnam_raw,     "HO3'");
   strcpy(pdb_out->resnam,        " DG ");
   pdb_out->bval =                 0.00  ;


   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_10)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 05";
   char filename[]       = "test_H_Hg.xml";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "HG  ");
   strcpy(pdb_out->atnam_raw,     "HG  ");
   strcpy(pdb_out->resnam,        " HG ");


   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_11)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 06";
   char filename[]       = "test_H_1ZK.xml";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "CD21");
   strcpy(pdb_out->atnam_raw,     "CD21");
   strcpy(pdb_out->resnam,        "1ZK ");


   /* write test file */
   fp = fopen(test_output_filename,"w");
   WritePDBML(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = compare_files(test_example_filename, 
                                   test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

/* Create Suite */
Suite *writepdbml_suite(void)
{
   Suite *s = suite_create("WritePDBML");

   /* Core test case */
   TCase *tc_core = tcase_create("Core");
   tcase_add_checked_fixture(tc_core, writepdbml_setup, 
                             writepdbml_teardown);
   tcase_add_test(tc_core, test_write_pdb);
   tcase_add_test(tc_core, test_write_pdbml);
   suite_add_tcase(s, tc_core);

   TCase *tc_pdb = tcase_create("PDB");
   tcase_add_checked_fixture(tc_pdb, writepdbml_setup, 
                             writepdbml_teardown);
   tcase_add_test(tc_pdb, test_write_pdb_data_01);
   tcase_add_test(tc_pdb, test_write_pdb_data_02);
   tcase_add_test(tc_pdb, test_write_pdb_data_03);
   tcase_add_test(tc_pdb, test_write_pdb_data_04);
   suite_add_tcase(s, tc_pdb);
   
   TCase *tc_pdbml = tcase_create("PDBML");
   tcase_add_checked_fixture(tc_pdbml, writepdbml_setup, 
                             writepdbml_teardown);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_01);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_02);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_03);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_04);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_05);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_06);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_07);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_08);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_09);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_10);
   tcase_add_test(tc_pdbml, test_write_pdbml_data_11);
   suite_add_tcase(s, tc_pdbml);

   return s;
}

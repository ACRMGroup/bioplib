#include "writepdbml_suite.h"

/* Globals */
//static char test_output_filename[]     = "/tmp/test-XXXXX",
static char test_output_filename[]     = "tmp_out/test-XXXXX",
            test_input_filename[160]   = "data/wholepdb_suite/",
            test_example_filename[160] = "data/wholepdb_suite/";

static FILE      *fp              =             NULL;
static WHOLEPDB  *wpdb            =             NULL;
static BOOL      files_identical  =            FALSE,
                 read_pdbml_flag  =            FALSE;
static int       force_pdbml_flag = FORCEXML_NOFORCE;

/* Compare file function */
static BOOL wholepdb_compare_files(char *filename_a, char *filename_b)
{
   char command[120];
   
   /* set command */
   sprintf(command,"cmp %s %s > /dev/null", filename_a, filename_b);
   
   /* return TRUE if files match */
   return system(command) == 0 ? TRUE:FALSE;
}

/* Setup And Teardown */
void wholepdb_setup(void)
{
   /* Set PDB/PDBML flags to default */
   force_pdbml_flag = gPDBXMLForce;
   read_pdbml_flag  = gPDBXML;
   gPDBXMLForce     = FORCEXML_NOFORCE;
}

void wholepdb_teardown(void)
{
   /* Reset PDB/PDBML flags */
   gPDBXML      = read_pdbml_flag;
   gPDBXMLForce = force_pdbml_flag;

   /* Free WPDB */
   if(wpdb) FreeWholePDB(wpdb);
}

/* Core tests */
START_TEST(test_read_pdb)
{
   /* set filename */
   char filename_in[] = "test_alanine_in.pdb";
   
   /* set gPDBXML flag */
   gPDBXML = TRUE;
   
   /* read test file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* tests */
   ck_assert_msg(wpdb != NULL,          "Failed to read PDB file.");
   ck_assert_msg(wpdb->pdb != NULL, "No atom data read from file.");
   ck_assert_msg(wpdb->natoms > 0,           "Atom count not set.");
   ck_assert_msg(gPDBXML == FALSE,        "PDB read flag not set.");
}
END_TEST

START_TEST(test_write_pdb_01)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.pdb",
        filename_example[] = "test_alanine_out_01.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST

START_TEST(test_write_pdbml_01)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.pdb",
        filename_example[] = "test_alanine_out.xml",
        test_message[]     = "Output PDBML does not match example file.";
        
   /* force write PDB */
   FORCEXML;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST

START_TEST(test_read_pdbml)
{
   /* set filename */
   char filename_in[] = "test_alanine_in.xml";
   
   /* set gPDBXML flag */
   gPDBXML = FALSE;
   
   /* read test file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* tests */
   ck_assert_msg(wpdb != NULL,          "Failed to read PDB file.");
   ck_assert_msg(wpdb->pdb != NULL, "No atom data read from file.");
   ck_assert_msg(wpdb->natoms > 0,           "Atom count not set.");
   ck_assert_msg(gPDBXML == TRUE,         "PDB read flag not set.");
}
END_TEST

START_TEST(test_write_pdb_02)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.xml",
        filename_example[] = "test_alanine_out_02.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST

START_TEST(test_write_pdbml_02)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.xml",
        filename_example[] = "test_alanine_out.xml",
        test_message[]     = "Output PDBML does not match example file.";
        
   /* force write PDB */
   FORCEXML;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST

START_TEST(test_read_write_pdb)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.pdb",
        filename_example[] = "test_alanine_out_01.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* Set Default */
   gPDBXMLForce = FORCEXML_NOFORCE;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST

START_TEST(test_read_write_pdbml)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in.xml",
        filename_example[] = "test_alanine_out.xml",
        test_message[]     = "Output PDBML does not match example file.";
        
   /* Set Default */
   gPDBXMLForce = FORCEXML_NOFORCE;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = ReadWholePDB(fp);
   fclose(fp);

   /* write output file */
   mkstemp(test_output_filename);
   fp = fopen(test_output_filename,"w");
   blWriteWholePDB(fp, wpdb);
   fclose(fp);

   /* compare output file to example file */
   strcat(test_example_filename, filename_example);
   files_identical = wholepdb_compare_files(test_example_filename, 
                                            test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_message);
}
END_TEST



/* Create Suite */
Suite *wholepdb_suite(void)
{
   Suite *s = suite_create("WholePDB");

   /* Core test case */
   TCase *tc_core = tcase_create("Core");
   tcase_add_checked_fixture(tc_core, 
                             wholepdb_setup, 
                             wholepdb_teardown);
   tcase_add_test(tc_core, test_read_pdb);
   tcase_add_test(tc_core, test_write_pdb_01);
   tcase_add_test(tc_core, test_write_pdbml_01);
   tcase_add_test(tc_core, test_read_pdbml);
   tcase_add_test(tc_core, test_write_pdb_02);
   tcase_add_test(tc_core, test_write_pdbml_02);
   tcase_add_test(tc_core, test_read_write_pdb);
   tcase_add_test(tc_core, test_read_write_pdbml);   
   suite_add_tcase(s, tc_core);

   return s;
}

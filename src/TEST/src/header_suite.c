/************************************************************************/
/**

   \file       header_suite.c
   
   \version    V1.0
   \date       05.05.15
   \brief      Test suite for header data for pdbml.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2015
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

   Test suite for reading and writing header data for pdbml files.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.05.15 Original By: CTP

*************************************************************************/

#include "header_suite.h"

/* Globals */
static char test_output_filename[]     = "tmp/test-XXXXX",
            test_example_basename[]    = "data/header_suite/",
            test_input_filename[160]   = "",
            test_example_filename[160] = "";

static FILE      *fp              =             NULL;
static WHOLEPDB  *wpdb            =             NULL;
static BOOL      files_identical  =            FALSE,
                 read_pdbml_flag  =            FALSE;
static int       force_pdbml_flag = FORCEXML_NOFORCE;

/* Compare file function */
static BOOL wholepdb_compare_files(char *filename_a, char *filename_b)
{
   char command[120];

#ifndef MS_WINDOWS

   /* compare files command */
   sprintf(command,"cmp %s %s > /dev/null", filename_a, filename_b);

#else   

   /* convert output file format from dos to unix */
   sprintf(command,"dos2unix -q %s", filename_b);
   system(command);
   strcpy(command,""); /* reset command */

   /* compare files command */
   sprintf(command,"cmp %s %s", filename_a, filename_b);

#endif

   /* return TRUE if files match */
   return system(command) == 0 ? TRUE:FALSE;
}

/* Setup And Teardown */
static void wholepdb_setup(void)
{
   /* Set PDB/PDBML flags to default */
   force_pdbml_flag = gPDBXMLForce;
   read_pdbml_flag  = gPDBXML;
   gPDBXMLForce     = FORCEXML_NOFORCE;

   /* Copy base name to input and example file name */
   strcpy(test_example_filename, test_example_basename);
   strcpy(test_input_filename,   test_example_basename);
}

static void wholepdb_teardown(void)
{
   /* Reset PDB/PDBML flags */
   gPDBXML      = read_pdbml_flag;
   gPDBXMLForce = force_pdbml_flag;

   /* Free WPDB */
   if(wpdb) blFreeWholePDB(wpdb);
}

/* Core tests */
/* Parse PDB and Write PDBML*/
START_TEST(test_read_pdb)
{
   /* set filename */
   char filename_in[] = "test_alanine_in_01.pdb";
   
   /* set gPDBXML flag */
   gPDBXML = TRUE;
   
   /* read test file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);

   /* tests */
   ck_assert_msg(wpdb != NULL,          "Failed to read PDB file.");
   ck_assert_msg(wpdb->pdb != NULL, "No atom data read from file.");
   ck_assert_msg(wpdb->natoms > 0,           "Atom count not set.");
   ck_assert_msg(gPDBXML == FALSE,        "PDB read flag not set.");
}
END_TEST

START_TEST(test_write_pdbml_01)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_01.pdb",
        filename_example[] = "test_alanine_out_01.xml",
        test_message[]     = "Output PDBML does not match example file.";
        
   /* force write PDBML */
   FORCEXML;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* Fix: Set flag for XML */
   gPDBXML = TRUE;
   
   /* write output file */
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


/* CONVERT MULTI-CHAIN PDB COMPND TO PDBML */
/* NEED TEST FILES */
START_TEST(test_write_pdbml_02)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_01.pdb",
//   char filename_in[]      = "1ZEH.pdb",
        filename_example[] = "test_alanine_out_01.xml",
        test_message[]     = "Output PDBML does not match example file.";

   /* force write PDBML */
   FORCEXML;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* Fix: Set flag for XML */
   gPDBXML = TRUE;
   
   /* write output file */
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


/* Parse PDBML - Write PDB */
START_TEST(test_read_pdbml)
{
   /* set filename */
   char filename_in[] = "test_alanine_in_01.xml";
   
   /* set gPDBXML flag */
   gPDBXML = FALSE;
   
   /* read test file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);

   /* tests */
   ck_assert_msg(wpdb != NULL,          "Failed to read PDB file.");
   ck_assert_msg(wpdb->pdb != NULL, "No atom data read from file.");
   ck_assert_msg(wpdb->natoms > 0,           "Atom count not set.");
   ck_assert_msg(gPDBXML == TRUE,         "PDB read flag not set.");
}
END_TEST


START_TEST(test_write_pdb_01)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_01.xml",
        filename_example[] = "test_alanine_out_01.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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

/* TEST FOR MULTI-CHAIN COMPND */
/* REQUIRES TEST FILES */
START_TEST(test_write_pdb_02)
{
   /* get pdb data */
//   char filename_in[]      = "1ZEH.xml",
   char filename_in[]      = "test_alanine_in_01.xml",
        filename_example[] = "test_alanine_out_01.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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


/* Title tests */
/* =========== */

/* single-line title pdb -> pdb */
START_TEST(test_title_01)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_02.pdb",
        filename_example[] = "test_alanine_in_02.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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

/* multi-line title pdb -> pdb */
START_TEST(test_title_02)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_03.pdb",
        filename_example[] = "test_alanine_in_03.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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

/* multi-line title xml -> pdb */
START_TEST(test_title_03)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_02.xml",
        filename_example[] = "test_alanine_out_02.pdb",
        test_message[]     = "Output PDB does not match example file.";
        
   /* force write PDB */
   FORCEPDB;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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

/* multi-line title pdb -> xml */
START_TEST(test_title_04)
{
   /* get pdb data */
   char filename_in[]      = "test_alanine_in_03.pdb",
        filename_example[] = "test_alanine_out_02.xml",
        test_message[]     = "Output PDBML does not match example file.";
        
   /* force write PDB */
   FORCEXML;
   
   /* read input file */
   strcat(test_input_filename,filename_in);
   fp = fopen(test_input_filename,"r");
   wpdb = blReadWholePDB(fp);
   fclose(fp);
   ck_assert_msg(wpdb != NULL, "Failed to read PDB file.");

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif

   /* write output file */
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
Suite *header_suite(void)
{
   Suite *s = suite_create("Header");
   TCase *tc_core  = tcase_create("Core");
   TCase *tc_title = tcase_create("Title");

   /* Core test case */
   tcase_add_checked_fixture(tc_core, 
                             wholepdb_setup, 
                             wholepdb_teardown);

   /* read pdb and write pdbml */
   tcase_add_test(tc_core, test_read_pdb);
   tcase_add_test(tc_core, test_write_pdbml_01);
   tcase_add_test(tc_core, test_write_pdbml_02);

   /* read pdbml and write pdb */
   tcase_add_test(tc_core, test_read_pdbml);
   tcase_add_test(tc_core, test_write_pdb_01);
   tcase_add_test(tc_core, test_write_pdb_02);

   suite_add_tcase(s, tc_core);

   /* Title test case */
   tcase_add_checked_fixture(tc_title, 
                             wholepdb_setup, 
                             wholepdb_teardown);

   /* Check title */
   tcase_add_test(tc_title, test_title_01);
   tcase_add_test(tc_title, test_title_02);
   tcase_add_test(tc_title, test_title_03);
   tcase_add_test(tc_title, test_title_04);

   suite_add_tcase(s, tc_title);


   return s;
}

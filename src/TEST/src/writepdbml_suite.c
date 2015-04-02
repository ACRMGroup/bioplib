/************************************************************************/
/**

   \file       writepdbml_suite.c
   
   \version    V1.5
   \date       02.04.15
   \brief      Test suite for writing pdb and pdbml data to file.
   
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

   Test suite for writing pdb and pdbml data to file.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.08.14 Original By: CTP
-  V1.1  16.08.14 Test write charge and element to file. By: CTP
-  V1.2  26.08.14 Set ATOM record type to "ATOM  ". By: CTP
-  V1.3  29.08.14 Check file format before reading as format check no 
                  longer rewinds. By: CTP
-  V1.4  12.09.14 Update tests for MS Windows. By: CTP
-  V1.5  02.04.15 Update tests for segment ID. By: CTP

*************************************************************************/

#include "writepdbml_suite.h"

/* Globals */
static char test_output_filename[]     = "tmp/test-XXXXX",
            test_example_basename[]    = "data/writepdbml_suite/",
            test_example_filename[160] = "";

static FILE *fp               =             NULL;
static PDB  *pdb_out          =             NULL,
            *pdb_in           =             NULL;
static BOOL files_identical   =            FALSE,
            read_pdbml_flag   =            FALSE;
static int  force_pdbml_flag  = FORCEXML_NOFORCE;


/* Compare file function */
static BOOL writepdbml_compare_files(char *filename_a, char *filename_b)
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

/* Set test data function */
static void writepdbml_set_test_data(void)
{
   /* Output PDB */
   pdb_out = (PDB *) malloc(sizeof(PDB));
   strcpy(pdb_out->record_type, "ATOM  ");
   pdb_out->atnum =                    1 ;
   strcpy(pdb_out->atnam,         "CA  ");
   strcpy(pdb_out->atnam_raw,     " CA ");
   strcpy(pdb_out->resnam,        "ALA ");
   strcpy(pdb_out->chain,            "A");
   pdb_out->resnum =                   1 ;
   strcpy(pdb_out->insert,           " ");
   pdb_out->x =                    1.000 ;
   pdb_out->y =                    2.000 ;
   pdb_out->z =                    3.000 ;
   pdb_out->occ =                   1.00 ;
   pdb_out->altpos =                 ' ' ;
   pdb_out->bval =                 20.00 ;
   strcpy(pdb_out->segid,         "    ");
   strcpy(pdb_out->element,          "C");
   pdb_out->formal_charge =            0 ;
   pdb_out->partial_charge =         0.0 ;
   pdb_out->next =                  NULL ;

   /* Copy base name to example file name */
   strcpy(test_example_filename, test_example_basename);

#ifndef MS_WINDOWS   
   /* Set temp file name */
   mkstemp(test_output_filename);
#endif
}

/* Setup And Teardown */
static void writepdbml_setup_default(void)
{
   /* Set PDB/PDBML flags to default */
   force_pdbml_flag = gPDBXMLForce;
   read_pdbml_flag  = gPDBXML;
   gPDBXMLForce     = FORCEXML_NOFORCE;
   
   /* Set test data */
   writepdbml_set_test_data();
}

static void writepdbml_setup_pdb(void)
{
   /* Set PDB/PDBML flags to force PDB */
   force_pdbml_flag = gPDBXMLForce;
   read_pdbml_flag  = gPDBXML;
   gPDBXMLForce     = FORCEXML_PDB;
   
   /* Set test data */
   writepdbml_set_test_data();
}

static void writepdbml_setup_pdbml(void)
{
   /* Set PDB/PDBML flags to force PDBML */
   force_pdbml_flag = gPDBXMLForce;
   read_pdbml_flag  = gPDBXML;
   gPDBXMLForce     = FORCEXML_XML;
   
   /* Set test data */
   writepdbml_set_test_data();
}

static void writepdbml_teardown(void)
{
   /* Reset PDB/PDBML flags */
   gPDBXML      = read_pdbml_flag;
   gPDBXMLForce = force_pdbml_flag;

   /* Free PDB */
   FREELIST(pdb_out,PDB);
   FREELIST(pdb_in ,PDB);
}

/* Core tests */
START_TEST(test_write_pdb)
{
   int natoms        =     0;
   BOOL pdbml_format = FALSE;

   /* set format */
   FORCEPDB;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,                "Failed to open file.");

   /* check test file format */
   pdbml_format = blCheckFileFormatPDBML(fp);

   /* read test file */
   pdb_in = blReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in       !=  NULL,     "Failed to read file.");
   ck_assert_msg(natoms       ==     1, "No atoms read from file.");
   ck_assert_msg(pdbml_format == FALSE,   "File is not pdb format");
}
END_TEST

START_TEST(test_write_pdbml)
{
   int natoms        =     0;
   BOOL pdbml_format = FALSE;

   /* set format */
   FORCEXML;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,                "Failed to open file.");
   
   /* check test file format */
   pdbml_format = blCheckFileFormatPDBML(fp);

   /* read test file */
   pdb_in = blReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in       !=  NULL,     "Failed to read file.");
   ck_assert_msg(natoms       ==     1, "No atoms read from file.");
   ck_assert_msg(pdbml_format ==  TRUE, "File is not pdbml format");
}
END_TEST


START_TEST(test_write_default_pdb_in)
{
   int natoms        =     0;
   BOOL pdbml_format = FALSE;

   /* set format */
   gPDBXML      = FALSE;
   gPDBXMLForce = FORCEXML_NOFORCE;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,                "Failed to open file.");

   /* check test file format */
   pdbml_format = blCheckFileFormatPDBML(fp);

   /* read test file */
   pdb_in = blReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in       !=  NULL,     "Failed to read file.");
   ck_assert_msg(natoms       ==     1, "No atoms read from file.");
   ck_assert_msg(pdbml_format == FALSE,   "File is not pdb format");
}
END_TEST

START_TEST(test_write_default_pdbml_in)
{
   int natoms        =     0;
   BOOL pdbml_format = FALSE;

   /* set format */
   gPDBXML      = TRUE; /* pdb input file */
   gPDBXMLForce = FORCEXML_NOFORCE;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);
   
   /* open test file */
   fp = fopen(test_output_filename,"r");
   ck_assert_msg(fp != NULL,                "Failed to open file.");

   /* check test file format */
   pdbml_format = blCheckFileFormatPDBML(fp);

   /* read test file */
   pdb_in = blReadPDB(fp, &natoms);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* tests */
   ck_assert_msg(pdb_in       !=  NULL,     "Failed to read file.");
   ck_assert_msg(natoms       ==     1, "No atoms read from file.");
   ck_assert_msg(pdbml_format ==  TRUE, "File is not pdbml format");
}
END_TEST


/* Format Check */
START_TEST(test_pdb_format_check_valid)
{
   /* Set format */
   FORCEPDB;

   /* Update PDB */
   strcpy(pdb_out->chain, "X"); /* Single-letter chain ID */

   /* remove unused output file */
   remove(test_output_filename);

   /* Test */
   ck_assert_msg(blFormatCheckWritePDB(pdb_out) == TRUE,
                 "Format check for single-letter chain id.");
}
END_TEST

START_TEST(test_pdb_format_check_invalid)
{
   /* Set format */
   FORCEPDB;

   /* Update PDB */
   strcpy(pdb_out->chain, "XXX"); /* Multi-letter chain ID */

   /* remove unused output file */
   remove(test_output_filename);

   /* Test */
   ck_assert_msg(blFormatCheckWritePDB(pdb_out) == FALSE,
                 "Format check for multi-letter chain id.");
}
END_TEST

START_TEST(test_pdb_format_error_valid)
{
   BOOL return_value;

   /* Set format */
   FORCEPDB;

   /* Update PDB */
   strcpy(pdb_out->chain, "X"); /* Multi-letter chain ID */

   /* Attempt to write test file */
   fp = fopen(test_output_filename,"w");
   return_value = blWritePDB(fp, pdb_out);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* Test */
   ck_assert(return_value == TRUE);
}
END_TEST

START_TEST(test_pdb_format_error_invalid)
{
   BOOL return_value;

   /* Set format */
   FORCEPDB;

   /* Update PDB */
   strcpy(pdb_out->chain, "XXX"); /* Multi-letter chain ID */

   /* Attempt to write test file */
   fp = fopen(test_output_filename,"w");
   return_value = blWritePDB(fp, pdb_out);
   fclose(fp);

   /* remove output file */
   remove(test_output_filename);

   /* Test */
   ck_assert(return_value == FALSE);
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "ZN");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdb_data_05)
{
   /* set message and example file */
   char test_error_msg[] = "Write element (pdb)";
   char filename[]       = "test_alpha_carbon_element.pdb";

   /* Update PDB */
   strcpy(pdb_out->element,         "XX");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


/* PDB Ion */
START_TEST(test_write_pdb_ion_data_01)
{
   /* set message and example file */
   char test_error_msg[] = "Write chloride (pdb)";
   char filename[]       = "test_chloride.pdb";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "CL  ");
   strcpy(pdb_out->atnam_raw,     "CL  ");
   strcpy(pdb_out->resnam,        " CL ");
   strcpy(pdb_out->element,         "CL");
   pdb_out->formal_charge =           -1 ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdb_ion_data_02)
{
   /* set message and example file */
   char test_error_msg[] = "Write heme iron (pdb)";
   char filename[]       = "test_heme_iron.pdb";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "FE  ");
   strcpy(pdb_out->atnam_raw,     "FE  ");
   strcpy(pdb_out->resnam,        "HEM ");
   strcpy(pdb_out->element,         "FE");
   pdb_out->formal_charge =            3 ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "ZN");

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "  ");
   pdb_out->bval =                 0.00  ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "  ");
   pdb_out->bval =                 0.00  ;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


START_TEST(test_write_pdbml_data_08)
{
   /* set message and example file */
   char test_error_msg[] = "Write hydrogen 03";
   char filename[]       = "test_H_leu_delta.xml";

   /* Update PDB */
   strcpy(pdb_out->atnam,         "HD11");
   strcpy(pdb_out->atnam_raw,     "HD11");
   strcpy(pdb_out->resnam,        "LEU ");
   strcpy(pdb_out->element,         "  ");
   pdb_out->bval =                 0.00  ;
   
   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "  ");
   pdb_out->bval =                 0.00  ;


   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "  ");


   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   strcpy(pdb_out->element,         "  ");


   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_data_12)
{
   /* set message and example file */
   char test_error_msg[] = "Write element (pdbml)";
   char filename[]       = "test_alpha_carbon_element.xml";

   /* Update PDB */
   strcpy(pdb_out->element,         "XX");


   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST


/* PDBML Ion */
START_TEST(test_write_pdbml_ion_data_01)
{
   /* set message and example file */
   char test_error_msg[] = "Write chloride (pdbml)";
   char filename[]       = "test_chloride.xml";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "CL  ");
   strcpy(pdb_out->atnam_raw,     "CL  ");
   strcpy(pdb_out->resnam,        " CL ");
   strcpy(pdb_out->element,         "CL");
   pdb_out->formal_charge =           -1 ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
                                              test_output_filename);

   /* remove output file */
   remove(test_output_filename);
  
   /* return test result */
   ck_assert_msg(files_identical, test_error_msg);
}
END_TEST

START_TEST(test_write_pdbml_ion_data_02)
{
   /* set message and example file */
   char test_error_msg[] = "Write heme iron (pdbml)";
   char filename[]       = "test_heme_iron.xml";

   /* Update PDB */
   strcpy(pdb_out->record_type, "HETATM");
   strcpy(pdb_out->atnam,         "FE  ");
   strcpy(pdb_out->atnam_raw,     "FE  ");
   strcpy(pdb_out->resnam,        "HEM ");
   strcpy(pdb_out->element,         "FE");
   pdb_out->formal_charge =            3 ;

   /* write test file */
   fp = fopen(test_output_filename,"w");
   blWritePDB(fp, pdb_out);
   fclose(fp);

   /* compare to example file */
   strcat(test_example_filename, filename);
   files_identical = writepdbml_compare_files(test_example_filename, 
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
   Suite *s            = suite_create("WritePDBML");
   TCase *tc_core      = tcase_create("Core"),
         *tc_format    = tcase_create("Format"),
         *tc_pdb       = tcase_create("PDB"),
         *tc_pdb_ion   = tcase_create("PDB_ion"),
         *tc_pdbml     = tcase_create("PDBML"),
         *tc_pdbml_ion = tcase_create("PDBML_ion");


   /* Core test case */
   tcase_add_checked_fixture(tc_core, 
                             writepdbml_setup_default, 
                             writepdbml_teardown);
   tcase_add_test(tc_core, test_write_pdb);
   tcase_add_test(tc_core, test_write_pdbml);
   tcase_add_test(tc_core, test_write_default_pdb_in);
   tcase_add_test(tc_core, test_write_default_pdbml_in);
   suite_add_tcase(s, tc_core);

   /* Format test case */
   tcase_add_checked_fixture(tc_format, 
                             writepdbml_setup_default, 
                             writepdbml_teardown);
   tcase_add_test(tc_format, test_pdb_format_check_valid);
   tcase_add_test(tc_format, test_pdb_format_check_invalid);
   tcase_add_test(tc_format, test_pdb_format_error_valid);
   tcase_add_test(tc_format, test_pdb_format_error_invalid);
   suite_add_tcase(s, tc_format);

   /* PDB test case */
   tcase_add_checked_fixture(tc_pdb, writepdbml_setup_pdb, 
                             writepdbml_teardown);
   tcase_add_test(tc_pdb, test_write_pdb_data_01);
   tcase_add_test(tc_pdb, test_write_pdb_data_02);
   tcase_add_test(tc_pdb, test_write_pdb_data_03);
   tcase_add_test(tc_pdb, test_write_pdb_data_04);
   tcase_add_test(tc_pdb, test_write_pdb_data_05);
   suite_add_tcase(s, tc_pdb);
   
   /* PDB_ion test case */
   tcase_add_checked_fixture(tc_pdb_ion, writepdbml_setup_pdb, 
                             writepdbml_teardown);
   tcase_add_test(tc_pdb_ion, test_write_pdb_ion_data_01);
   tcase_add_test(tc_pdb_ion, test_write_pdb_ion_data_02);
   suite_add_tcase(s, tc_pdb_ion);

   /* PDBML test case */
   tcase_add_checked_fixture(tc_pdbml, writepdbml_setup_pdbml, 
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
   tcase_add_test(tc_pdbml, test_write_pdbml_data_12);
   suite_add_tcase(s, tc_pdbml);

   /* PDBML_ion test case */
   tcase_add_checked_fixture(tc_pdbml_ion, writepdbml_setup_pdbml, 
                             writepdbml_teardown);
   tcase_add_test(tc_pdbml_ion, test_write_pdbml_ion_data_01);
   tcase_add_test(tc_pdbml_ion, test_write_pdbml_ion_data_02);
   suite_add_tcase(s, tc_pdbml_ion);

   return s;
}

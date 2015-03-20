/************************************************************************/
/**

   \file       ReadPDB.c
   
   \version    V3.3
   \date       20.03.15
   \brief      Read coordinates from a PDB file 
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1988-2015
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

\code
   pdb = ReadPDB(fp,natom)
\endcode

   This subroutine will read a .PDB file
   of any size and form a linked list of the protein structure.
   This list is contained in a linked set of structures of type
   pdb_entry. The structure is set up by including the file
   "pdb.h". For details of the structure, see this file.

   To free the space created by this routine, call FREELIST(pdb,PDB).

   The parameters passed to the subroutine are:

   -   fp    - A pointer to type FILE in which the .PDB file is stored.
   -   pdb   - A pointer to type PDB.
   -   natom - A pointer to type integer in which the number of atoms
               found is stored.

   As of V2.3, the routine makes provision for partial occupancies. If 
   the occupancies are 1.0 or 0.0, the atoms are read verbatim. If not,
   only the highest occupancy atoms are read and the atom names are 
   corrected to remove alternative labels. This behaviour can be 
   overridden by calling one of the ...OccRank() routines to read lower 
   occupancy atoms. If any partial occupancy atoms are read the global
   flag gPDBPartialOcc is set to TRUE.

   The various PDB reading routines set the following global flags:
   gPDBPartialOcc    - the PDB file contained multiple occupancies
   gPDBMultiNMR      - the PDB file contained multiple models
   gPDBXML           - the file was in PDBML (XML) format
   gPDBModelNotFound - the requested model was not found
   

NOTE:  Although some of the fields are represented by a single character,
       they are still stored in character arrays.

BUGS:  The subroutine cannot read files with VAX Fortran carriage control!
       It just sits there and page faults like crazy.

BUGS:  The multiple occupancy code assumes that all positions for a given
       atom in consecutive records of the file

BUGS:  25.01.05 Note the multiple occupancy code won't work properly for
       3pga where atoms have occupancies of zero and one

**************************************************************************

   Usage:
   ======

\code
   pdb = blReadPDB(fp,natom)
\endcode

   \param[in]      *fp      A pointer to type FILE in which the
                            .PDB file is stored.

   \param[out]     *natom   Number of atoms read.
   
   \return         *pdb     A pointer to the first allocated item of
                            the PDB linked list

**************************************************************************

   Revision History:
   =================
-  V1.0  04.11.88 Original
-  V1.1  07.02.89 Now ignores any records from the .PDB file which 
                  don't start with ATOM or HETATM.
-  V1.2  28.03.90 Some fields altered to match the exact specifications 
                  of the PDB. The only differences from the standard 
                  are:
                  1. The residue name is 4 characters rather than 3 
                     (allowing LYSH, HISA, etc.).
                  2. The atom name starts one column later than the 
                     standard and is four columns wide encompasing the 
                     standard's `alternate' field. These two 
                     differences from the standard reflect the common
                     usage.
-  V1.2a 28.06.90 Buffer size increased to 85 chars.
-  V1.2b 15.02.91 Simply changed comment header to match new standard.
-  V1.3  07.01.92 Corrected small bug in while() loop. Now ignores 
                  blank lines properly
-  V1.4  11.05.92 Added check on EOF in while() loop and memset() of 
                  buffer. ANSIfied.
-  V1.5  01.06.92 Documented for autodoc
-  V1.7  01.10.92 Changed to use fgets()
-  V1.6  19.06.92 Corrected use of stdlib
-  V1.8  08.12.92 SAS/C V6 now defines atof() in stdlib
-  V1.9  10.06.93 Returns TRUE or FALSE rather than exiting on failure
-  V2.0  17.06.93 Rewritten to use fsscanf()
-  V2.1  08.07.93 Modified to give ReadPDB() and ReadPDBAtoms()
-  V2.2  09.07.93 Modified to return the PDB pointer rather than a BOOL.
                  There is now no need to initialise the structure first.
                  Rewrote allocation scheme.
-  V2.3  17.03.94 Handles partial occupancies. If occupancies are not
                  1.0 or 0.0, the normal routine now reads only the 
                  highest occupancy atoms and corrects the atoms names 
                  to remove alternative labels. This behaviour can be 
                  overridden by calling one of the ...OccRank()
                  routines to read lower occupancy atoms. 
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Added ReadPDBOccRank() & ReadPDBAtomsOccRank()
                  Sets gPDBPartialOcc flag.
-  V2.4  06.04.94 With atom names which start in column 13, now checks
                  if the first character is a digit. If so, moves it
                  to the end of the atom name. Thus, 1HH1 becomes HH11
                  and 2HH1 becomes HH12.
-  V2.5  04.10.94 Fixed partial occ when resnum changes as well as atom
                  name. Fixed bug when MAXPARTIAL exceeded.
-  V2.6  03.11.94 Simply Corrected description. No code changes
-  V2.7  06.03.95 Now reads just the first NMR model by default
                  doReadPDB() no longer static
                  Sets gPDBMultiNMR if ENDMDL records found.
-  V2.8  13.01.97 Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
-  V2.9  25.02.98 Added transparent reading of gzipped PDB files if
                  GUNZIP_SUPPORT is defined
-  V2.10 18.08.98 Added cast to popen() for SunOS
-  V2.11 08.10.99 Initialised some variables
-  V2.12 15.02.01 Added atnam_raw into PDB structure
-  V2.13 30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V2.14 27.04.05 Fixed bug in atnam_raw for multiple occupancies
-  V2.15 03.06.05 Added altpos field to PDB structure. The massaged atom
                  name no longer contains the alternate indicator and
                  atnam_raw has only the atom name with altpos having the
                  alternate indicator (as it should!)
-  V2.16 14.10.05 Fixed a problem in StoreOccRankAtom() when a lower
                  occupancy atom has (erroneously) been set to occupancy
                  of zero and you want to pull out that atom
-  V2.17 25.01.06 Added calls to RemoveAlternates()
-  V2.18 03.02.06 Added prototypes for popen() and pclose()
-  V2.19 05.06.07 Added support for Unix compress'd files
-  V2.20 29.06.07 popen() and pclose() prototypes now skipped for MAC OSX
                  which defines them differently
-  V2.21 17.03.09 popen() prototype skipped for Windows. By: CTP
-  V2.22 21.12.11 doReadPDB() modified for cases where atoms are single
                  occupancy but occupancy is < 1.0
-  V2.23 04.02.14 Use CHAINMATCH macro. By: CTP
-  V2.24 22.04.14 Added PDBML parsing with doReadPDBML() and 
                  CheckFileFormatPDBML(). By CTP
-  V2.25 02.06.14 Updated doReadPDBML(). By: CTP
-  V2.26 09.06.14 Set gPDBXML flag. By: CTP
-  V2.27 07.07.14 Renaming of functions with "bl" prefix. By: CTP
-  V2.28 04.08.14 blReadPDB() and blReadPDBML() get element and charge.
                  Set access and radius to 0.0. Set atomType to NULL.
                  Added blProcessElementField() and blProcessChargeField()
                  By: CTP
-  V2.29 15.08.14 Updated blDoReadPDB() and blDoReadPDBML() to use 
                  CLEAR_PDB().  By: CTP
-  V2.30 16.08.14 Replaced charge with formal_charge and partial_charge 
                  for PDB structure. By: CTP
-  V2.31 18.08.14 Added XML_SUPPORT option allowing compilation without 
                  support for PDBML format. By: CTP
-  V2.32 26.08.14 blDoReadPDBML() pads record type to six chars. By: CTP
-  V2.33 29.08.14 Rewrote blCheckFileFormatPDBML() to take sample from 
                  input steam then push sample back on stream. By: CTP
-  V2.34 31.08.14 Fixed bug in blCheckFileFormatPDBML(). By: CTP
-  V2.35 09.09.14 Updated blCheckFileFormatPDBML() for non-unix systems.
                  Decreased size of XML_SAMPLE.
                  Reading of gzipped files with gunzip not supported for 
                  MS Windows. By: CTP
-  V2.36 29.09.14 Allow single character check for filetype where ungetc()
                  fails after pushback of single character. Updates to 
                  blCheckFileFormatPDBML() and blDoReadPDB(). By: CTP
-  V2.37 17.02.15 Added segid support   By: ACRM
-  V3.0  20.02.15 Merged functionality of ReadWholePDB() into this
-  V3.1  25.02.15 Enhanced error checking and safety of blDoReadPDBML()
-  V3.2  05.03.15 Fixed core dump in StoreConectRecords() when alternates
                  have been deleted so CONECT specifies an atom that
                  doesn't exist any more
-  V3.3  20.03.15 Fixed reading of whole PDB when not the first model
                  Now sets gPDBModelNotFound flag if the requested model
                  isn't found

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP File IO

   #KEYFUNCTION blReadPDB() 
   Main way of reading a PDB file into a linked list, reading just the
   highest occupancy atoms

   #FUNCTION blReadPDBAll() 
   Reads a PDB file into a linked list, reading all multiple
   occupancy atoms

   #FUNCTION blReadPDBAtoms() 
   Reads only ATOM records from a PDB file into a linked list, reading 
   just the highest occupancy atoms

   #FUNCTION blReadPDBOccRank() 
   Reads the specified ranking of occupancy (e.g. the second most 
   populated coordinates) from a PDB file into a linked list

   #FUNCTION blReadPDBAtomsOccRank() 
   Reads only the ATOM records for the specified ranking of occupancy 
   (e.g. the second most populated coordinates) from a PDB file into a 
   linked list

   #FUNCTION blDoReadPDB() 
   A lower level routine giving full control over reading all or only
   ATOM records, occupancy rankings and model numbers.

   #FUNCTION blDoReadPDBML() 
   A lower level routine giving full control over reading all or only
   ATOM records, occupancy rankings and model numbers from a PDBML XML
   file.

   #FUNCTION blCheckFileFormatPDBML() 
   A simple test to detect whether a file is a PDBML-formatted PDB file.

   #KEYFUNCTION  blReadWholePDB()
   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   #FUNCTION  blFreeWholePDB()
   Frees the header, trailer and atom content from a WHOLEPDB structure

   #FUNCTION  blReadWholePDBAtoms()
   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Only reads the ATOM record for 
   coordinates


   #SUBGROUP Atom names and elements
   #FUNCTION blFixAtomName()
   Fixes an atom name by removing leading spaces, or moving a leading
   digit to the end of the string.


   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION blRemoveAlternates()
   Removes alternate occupancy atoms. This may be useful after
   blReadPDBAll()
*/

/************************************************************************/
/* Defines required for includes
*/
#define READPDB_MAIN

/************************************************************************/
/* Includes
*/
#include "port.h"    /* Required before stdio.h                         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#ifdef XML_SUPPORT /* Required to read PDBML files                      */
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif

#include "SysDefs.h"
#include "MathType.h"
#include "pdb.h"
#include "macros.h"
#include "fsscanf.h"
#include "general.h"

#define MAXPARTIAL 8
#define SMALL      0.000001
#define XML_BUFFER 1024
#define XML_SAMPLE 256
#define MAXBUFF    160

#define LOCATION_HEADER      0
#define LOCATION_COORDINATES 1
#define LOCATION_TRAILER     2

/************************************************************************/
/* Prototypes
*/
static BOOL blStoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                               int NPartial, PDB **ppdb, PDB **pp, 
                               int *natom);
static void blProcessElementField(char *element, char *element_field);
static void blProcessChargeField(int *charge, char *charge_field);
static void StoreConectRecords(WHOLEPDB *wpdb, char *buffer);
#ifdef XML_SUPPORT
static BOOL blSetPDBDateField(char *pdb_date, char *pdbml_date);
static STRINGLIST *blParseHeaderPDBML(FILE *fpin);
#endif



#if !defined(__APPLE__) && !defined(MS_WINDOWS)
FILE *popen(char *, char *);
#endif
#ifndef __APPLE__
int  pclose(FILE *);
#endif

/************************************************************************/
/*>PDB *blReadPDB(FILE *fp, int *natom)
   ------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           .PDB file is stored.
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDB file into a PDB linked list

-  08.07.93 Written as entry for doReadPDB()
-  09.07.93 Modified to return pointer to PDB
-  17.03.94 Modified to handle OccRank
-  06.03.95 Added value for NMR model to read (1 = first)
-  25.01.06 Added call to RemoveAlternates() - this deals with odd
            cases where alternate atom positions don't appear where
            they should!
-  25.01.06 Added call to RemoveAlternates(). This deals with odd uses
            of multiple occupancies like 3pga and the instance where
            the alternates are all grouped at the end of the file.
-  07.07.14 Renamed to blReadPDB() By: CTP
*/
PDB *blReadPDB(FILE *fp,
               int  *natom)
{
   PDB *pdb;
   WHOLEPDB *wpdb;
   *natom=(-1);

   if((wpdb = blDoReadPDB(fp, TRUE, 1, 1, FALSE))!=NULL)
   {
      blFreeStringList(wpdb->header);
      blFreeStringList(wpdb->trailer);
      *natom = wpdb->natoms;
      pdb = wpdb->pdb;
      free(wpdb);

      pdb = blRemoveAlternates(pdb);
   }
   
   return(pdb);
}

/************************************************************************/
/*>PDB *blReadPDBAll(FILE *fp, int *natom)
   ---------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           .PDB file is stored.
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDB file into a PDB linked list. Reads all partial occupancy
   atoms. Reads both ATOM and HETATM records.

-  04.10.94 Original    By: ACRM
-  06.03.95 Added value for NMR model to read (0 = all)   
-  07.07.14 Renamed to blReadPDBAll() By: CTP
*/
PDB *blReadPDBAll(FILE *fp,
             int  *natom)
{
   PDB *pdb;
   WHOLEPDB *wpdb;
   *natom=(-1);

   if((wpdb = blDoReadPDB(fp, TRUE, 0, 0, FALSE))!=NULL)
   {
      blFreeStringList(wpdb->header);
      blFreeStringList(wpdb->trailer);
      *natom = wpdb->natoms;
      pdb = wpdb->pdb;
      free(wpdb);
   }
   
   return(pdb);
}

/************************************************************************/
/*>PDB *blReadPDBAtoms(FILE *fp, int *natom)
   -----------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           .PDB file is stored.
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDB file into a PDB linked list. Atoms only (no HETATM cards).

-  08.07.93 Written as entry for doReadPDB()
-  09.07.93 Modified to return pointer to PDB
-  17.03.94 Modified to handle OccRank
-  06.03.95 Added value for NMR model to read (1 = first)
-  25.01.06 Added call to RemoveAlternates(). This deals with odd uses
            of multiple occupancies like 3pga and the instance where
            the alternates are all grouped at the end of the file.
-  07.07.14 Renamed to blReadPDBAtoms() By: CTP
*/
PDB *blReadPDBAtoms(FILE *fp,
                    int  *natom)
{
   PDB *pdb;
   WHOLEPDB *wpdb;
   *natom=(-1);

   if((wpdb = blDoReadPDB(fp, FALSE, 1, 1, FALSE))!=NULL)
   {
      blFreeStringList(wpdb->header);
      blFreeStringList(wpdb->trailer);
      *natom = wpdb->natoms;
      pdb = wpdb->pdb;
      free(wpdb);

      pdb = blRemoveAlternates(pdb);
   }
   
   return(pdb);
}

/************************************************************************/
/*>PDB *blReadPDBOccRank(FILE *fp, int *natom, int OccRank)
   --------------------------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           .PDB file is stored.
   \param[in]     OccRank  Occupancy ranking (>=1)
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDB file into a PDB linked list selecting the OccRank'th
   highest occupancy atoms

-  17.03.94 Original    By: ACRM
-  06.03.95 Added value for NMR model to read (1 = first)
-  07.07.14 Renamed to blDoReadPDB() By: CTP
*/
PDB *blReadPDBOccRank(FILE *fp, int *natom, int OccRank)
{
   PDB *pdb;
   WHOLEPDB *wpdb;
   *natom=(-1);

   if((wpdb = blDoReadPDB(fp, TRUE, OccRank, 1, FALSE))!=NULL)
   {
      blFreeStringList(wpdb->header);
      blFreeStringList(wpdb->trailer);
      *natom = wpdb->natoms;
      pdb = wpdb->pdb;
      free(wpdb);

      pdb = blRemoveAlternates(pdb);
   }
   
   return(pdb);
}

/************************************************************************/
/*>PDB *blReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank)
   -------------------------------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           .PDB file is stored.
   \param[in]     OccRank  Occupancy ranking (>=1)
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDB file into a PDB linked list ignoring HETATM records
   and selecting the OccRank'th highest occupancy atoms

-  17.03.94 Original    By: ACRM
-  06.03.95 Added value for NMR model to read (1 = first)
-  07.07.14 Renamed to blReadPDBAtomsOccRank() By: CTP
*/
PDB *blReadPDBAtomsOccRank(FILE *fp, int *natom, int OccRank)
{
   PDB *pdb;
   WHOLEPDB *wpdb;
   *natom=(-1);

   if((wpdb = blDoReadPDB(fp, FALSE, OccRank, 1, FALSE))!=NULL)
   {
      blFreeStringList(wpdb->header);
      blFreeStringList(wpdb->trailer);
      *natom = wpdb->natoms;
      pdb = wpdb->pdb;
      free(wpdb);

      pdb = blRemoveAlternates(pdb);
   }
   
   return(pdb);
}

/************************************************************************/
/*>WHOLEPDB *blDoReadPDB(FILE *fpin, BOOL AllAtoms, int OccRank,
                         int ModelNum, BOOL DoWhole)
   --------------------------------------------------------------------
*//**

   \param[in]     *fpin    A pointer to type FILE in which the
                           .PDB file is stored.
   \param[in]     AllAtoms TRUE:  ATOM & HETATM records
                           FALSE: ATOM records only
   \param[in]     OccRank  Occupancy ranking
   \param[in]     ModelNum NMR Model number (0 = all)
   \param[in]     DoWhole  Read the whole PDB file rather than just 
                           the ATOM/HETATM records.
   \return                 A pointer to a malloc'd WHOLEPDB structure

   Reads a PDB file into a PDB linked list. The OccRank value indicates
   occupancy ranking to read for partial occupancy atoms.
   If any partial occupancy atoms are read the global flag 
   gPDBPartialOcc is set to TRUE.

-  04.11.88 V1.0  Original
-  07.02.89 V1.1  Ignore records which aren't ATOM or HETATM
-  28.03.90 V1.2  Altered field widths to match PDB standard better
                  See notes above for deviations
-  28.06.90 V1.2a Buffer size increased to 85 chars.
-  15.02.91 V1.2b Changed comment header to match new standard.
-  07.01.92 V1.3  Ignores blank lines properly
-  11.05.92 V1.4  Check on EOF in while() loop, memset() buffer. 
                  ANSIed.
-  01.06.92 V1.5  Documented for autodoc
-  19.06.92 V1.6  Corrected use of stdlib
-  01.10.92 V1.7  Changed to use fgets()
-  10.06.93 V1.9  Returns 0 on failure rather than exiting
                  Replaced SIZE with sizeof(PDB) directly
-  17.06.93 V2.0  Rewritten to use fsscanf()
-  08.07.93 V2.1  Split from ReadPDB()
-  09.07.93 V2.2  Modified to return pointer to PDB. Rewrote allocation
                  scheme.
-  17.03.94 V2.3  Handles partial occupancies
                  Sets natom to -1 if there was an error to distinguish 
                  from no atoms.
                  Handles atom names which start in column 13 rather
                  than column 14. This is allowed in the standard, but
                  very rare.
                  Sets flag for partials.
-  06.04.94 V2.4  Atom names starting in column 13 have their first
                  character moved to the end if it is a digit.
-  03.10.94 V2.5  Check residue number as well as atom name when running
                  through alternative atoms for partial occupancy
                  Moved increment of NPartial, so only done if there
                  is space in the array. If OccRank is 0, all atoms are
                  read regardless of occupancy.
-  06.03.95 V2.7  Added value for NMR model to read (0 = all)
                  No longer static. Sets gPDBMultiNMR if ENDMDL records
                  found.
-  13.01.97 V2.8  Added check on return from fsscanf. Blank lines used
                  to result in duplication of the previous line since
                  fsscanf() does not reset the variables on receiving
                  a blank line. Also fixed in fsscanf().
-  25.02.98 V2.9  Added code to read gzipped PDB files transparently
                  when GUNZIP_SUPPORT is defined
-  17.08.98 V2.10 Added case to popen() for SunOS
-  08.10.99 V2.11 Initialise CurIns and CurRes
-  15.02.01 V2.12 Added atnam_raw
-  27.04.05 V2.14 Added another atnam_raw for multiple occupancies
-  03.06.05 V2.15 Added altpos
-  14.10.05 V2.16 Modified detection of partial occupancy. handles
                  residues like 1zeh/B16 where a lower partial is
                  erroneously set to zero
-  05.06.07 V2.19 Added support for Unix compress'd files
-  21.12.11 V2.22 Modified for cases of single occupancy < 1.0
-  22.04.14 V2.24 Call doReadPDBML() for PDBML-formatted PDB file. By: CTP
-  02.06.14 V2.25 Updated doReadPDBML(). By: CTP
-  09.06.14 V2.26 Set gPDBXML flag. By: CTP
-  07.07.14 V2.27 Renamed to blDoReadPDB() By: CTP
-  15.08.14 V2.29 Use CLEAR_PDB() to set default values. By: CTP
-  16.08.14 V2.30 Replaced charge with formal_charge and partial_charge 
                  for PDB structure. By: CTP
-  18.08.14 V2.31 Added XML_SUPPORT option allowing BiopLib to be compiled
                  without support for PDBML format. By: CTP
-  09.09.14 V2.35 Reading of gzipped files with gunzip not supported for 
                  MS Windows. By: CTP
-  29.09.14 V2.36 Allow single character filetype check for gzipped files.
                  By: CTP
-  17.02.15 V2.37 Added segid support   By: ACRM
-  20.02.15 V3.0  NOT COMPATIBLE WITH PREVIOUS VERSIONS. The functionality
                  of the old ReadWholePDB() is now integrated into this
                  function.
-  20.03.15 V3.3  Fixed behaviour with reading other than the first model
                  and now sets a global error flag if the requested model
                  was not found. Uses MODEL records rather than ENDMDL
                  records in counting models

We need to deal with freeing wpdb if we are returning null.
Also need to deal with some sort of error code
*/
WHOLEPDB *blDoReadPDB(FILE *fpin,
                      BOOL AllAtoms,
                      int  OccRank,
                      int  ModelNum,
                      BOOL DoWhole)
{
   char     record_type[8],
            atnambuff[8],
            *atnam,
            atnam_raw[8],
            resnam[8],
            chain[4],
            insert[4],
            segid[8],
            buffer[160],
            CurAtom[8],
            cmd[80],
            CurIns = ' ',
            altpos,
            element_buff[4] = "",
            charge_buff[4]  = "",
            element[4]      = "";
   int      atnum,
            resnum,
            CurRes = 0,
            NPartial,
            ModelCount = 0,
            charge = 0,
            inLocation = LOCATION_HEADER;
   FILE     *fp = fpin;
   double   x,y,z,
            occ,
            bval;
   PDB      *p,
            multi[MAXPARTIAL];   /* Temporary storage for partial occ   */
   WHOLEPDB *wpdb = NULL;
   BOOL     pdbml_format;
   

#if defined(GUNZIP_SUPPORT) && !defined(MS_WINDOWS)
   int      signature[3],
            ch;
   BOOL     gzipped_file = FALSE;
#  ifndef SINGLE_CHAR_FILECHECK
   int      i;
#  endif
#endif

   if((wpdb=(WHOLEPDB *)malloc(sizeof(WHOLEPDB)))==NULL)
      return(NULL);

   wpdb->pdb         = NULL;
   wpdb->header      = NULL;
   wpdb->trailer     = NULL;
   
   wpdb->natoms      = 0;
   CurAtom[0]        = '\0';
   NPartial          = 0;
   gPDBPartialOcc    = FALSE;
   gPDBMultiNMR      = 0;
   cmd[0]            = '\0';
   gPDBXML           = FALSE;
   gPDBModelNotFound = TRUE;  /* Assume we haven't found the model      */

#if defined(GUNZIP_SUPPORT) && !defined(MS_WINDOWS)
   /* See whether this is a gzipped file                                */
#  ifndef SINGLE_CHAR_FILECHECK
   /* Default three character filetype check                            */
   for(i=0; i<3; i++)
      signature[i] = fgetc(fpin);
   for(i=2; i>=0; i--)
      ungetc(signature[i], fpin);
   if(((signature[0] == (int)0x1F) &&    /* gzip                        */
       (signature[1] == (int)0x8B) &&
       (signature[2] == (int)0x08)) ||
      ((signature[0] == (int)0x1F) &&    /* 05.06.07 compress           */
       (signature[1] == (int)0x9D) &&
       (signature[2] == (int)0x90)))
   {
      gzipped_file = TRUE;
   }
#  else
   /* Single character filetype check                                   */
   signature[0] = fgetc(fpin);
   ungetc(signature[0], fpin);
   if(signature[0] == (int)0x1F) gzipped_file = TRUE;
#  endif

   if(gzipped_file)
   {
      /* It is gzipped so we'll open gunzip as a pipe and send the data
         through that into a temporary file
      */
      sprintf(cmd,"gunzip >/tmp/readpdb_%d",(int)getpid());
      if((fp = (FILE *)popen(cmd,"w"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
      while((ch=fgetc(fpin))!=EOF)
         fputc(ch, fp);
      pclose(fp);

      /* We now reopen the temporary file as our PDB input file         */
      sprintf(cmd,"/tmp/readpdb_%d",(int)getpid());
      if((fp = fopen(cmd,"r"))==NULL)
      {
         wpdb->natoms = (-1);
         return(NULL);
      }
   }
#endif   


   /* Check file format                                                 */
   pdbml_format = blCheckFileFormatPDBML(fp);
   
   /* If it's PDBML then call the appropriate parser                    */
   if(pdbml_format)
   {
#ifdef XML_SUPPORT
      /* Parse PDBML-formatted PDB file                                 */
      if(DoWhole)
         wpdb->header = blParseHeaderPDBML(fp);
      
      wpdb->pdb = blDoReadPDBML(fp,&wpdb->natoms,AllAtoms,
                                OccRank,ModelNum);
      if(DoWhole)
         wpdb->trailer = blStoreString(wpdb->trailer, "END   \n");
      xmlCleanupParser();     /* free globals set by parser             */
      if(cmd[0]) unlink(cmd); /* delete tmp file                        */
      return(wpdb);           /* return PDB list                        */
#else
      /* PDBML format not supported.                                    */
      if(cmd[0]) unlink(cmd); /* delete tmp file                        */
      wpdb->natoms = (-1);    /* Indicate error                         */
      return(NULL);           /* return NULL list                       */
#endif
   }

   inLocation = LOCATION_HEADER;
   
   while(fgets(buffer,159,fp))
   {
      /*** Deal with counting model numbers                           ***/
      if(ModelNum != 0)          /* We are interested in model numbers  */
      {
         if(!strncmp(buffer,"MODEL ",6))
         {
            ModelCount++;
            gPDBMultiNMR++;
         }

         /* See if we are in the right model                            */
         if(inLocation == LOCATION_COORDINATES)
         {
            if((ModelCount != ModelNum) && (ModelCount != 0))
               continue;
            else
               gPDBModelNotFound = FALSE;
         }
      }
      else
      {
         gPDBModelNotFound = FALSE;
      }
      

      if(!strncmp(buffer, "ATOM  ", 6) ||
         !strncmp(buffer, "HETATM", 6) ||
         !strncmp(buffer, "MODEL ", 6))
      {
         inLocation = LOCATION_COORDINATES;
      }
      else if(!strncmp(buffer, "CONECT", 6) ||
              !strncmp(buffer, "MASTER", 6) ||
              !strncmp(buffer, "END   ", 6))
      {
         inLocation = LOCATION_TRAILER;
      }
      
      /* If we are in the header, just store it                         */
      if(inLocation == LOCATION_HEADER)
      {
         if(DoWhole)
         {
            if((wpdb->header = blStoreString(wpdb->header, buffer))==NULL)
               return(NULL);
         }
         continue;
      }
      if(inLocation == LOCATION_TRAILER)
      {
         if(DoWhole)
         {
            wpdb->trailer = blStoreString(wpdb->trailer, buffer);
            if(!strncmp(buffer, "CONECT", 6))
               StoreConectRecords(wpdb, buffer);
         }
         
         continue;
      }

      /* Read a record                                                  */
      if(fsscanf(buffer,
                 "%6s%5d%1x%5s%4s%1s%4d%1s%3x%8lf%8lf%8lf%6lf%6lf%6x%4s%2s%2s",
                 record_type,&atnum,atnambuff,resnam,chain,&resnum,insert,
                 &x,&y,&z,&occ,&bval,segid,element_buff,charge_buff) 
         != EOF)
      {
         if((!strncmp(record_type,"ATOM  ",6)) || 
            (!strncmp(record_type,"HETATM",6) && AllAtoms))
         {
            /* Copy the raw atom name                                   */
            /* 03.06.05 Note: this reads the alternate atom position as 
               well as the atom name - changes in FixAtomName() now strip
               that
               We now copy only the first 4 characters into atnam_raw and
               put the 5th character into altpos
            */
            strncpy(atnam_raw, atnambuff, 4);
            atnam_raw[4] = '\0';
            altpos = atnambuff[4];

            /* Fix the atom name accounting for start in column 13 or 14*/
            atnam = blFixAtomName(atnambuff, occ);
            
            /* Set element and charge                                   */
            blProcessElementField(element, element_buff);
            blProcessChargeField(&charge, charge_buff);
            
            /* Set element from atom name if not in input file          */
            if(strlen(element) == 0)
            {
               blSetElementSymbolFromAtomName(element, atnam_raw);
            }

            /* Check for full occupancy. If occupancy is 0.0 assume that 
               it is actually fully occupied; the column just hasn't been
               filled in correctly
               
               04.10.94 Read all atoms if OccRank is 0

               14.10.05 Now takes an atom as full occupancy:
                           if occ==1.0
                           if occ==0.0 and altpos==' '
                           if OccRank==0
                        This fixes problems where a lower (partial)
                        occupancy has erroneously been set to zero
               21.12.11 Now only worries about partial occupancy if altpos
                        is a space. The first line of the if() statement
                        here would assume single occupancy if altpos was
                        a space and occupancy was zero:
                        if(((altpos == ' ') && (occ < (double)SMALL)) ||
                        - it now assumes single occupancy if altpos is a
                        space regardless of the actual occupancy. This
                        deals with cases like 1ap2 ZN A112 and 1ces ZN
                        A238 where these HETATMs are single occupancy
                        but with occupancy < 1.0
            */
            if((altpos == ' ') ||
               (occ > (double)0.999) || 
               (OccRank == 0))
            {
               /* Trim the atom name to 4 characters                    */
               atnam[4] = '\0';
               
               if(NPartial != 0)
               {
                  if(!blStoreOccRankAtom(OccRank,multi,NPartial,
                                         &wpdb->pdb,&p,&(wpdb->natoms)))
                  {
                     if(wpdb->pdb != NULL) FREELIST(wpdb->pdb, PDB);
                     wpdb->natoms = (-1);
                     if(cmd[0]) unlink(cmd);
                     return(NULL);
                  }
                  
                  /* Set partial occupancy counter to 0                 */
                  NPartial = 0;
               }
               
               /* Allocate space in the linked list                     */
               if(wpdb->pdb == NULL)
               {
                  INIT(wpdb->pdb, PDB);
                  p = wpdb->pdb;
               }
               else
               {
                  ALLOCNEXT(p, PDB);
               }
               
               /* Failed to allocate space; free up list so far & return*/
               if(p==NULL)
               {
                  if(wpdb->pdb != NULL) FREELIST(wpdb->pdb, PDB);
                  wpdb->natoms = (-1);
                  if(cmd[0]) unlink(cmd);
                  return(NULL);
               }
               
               /* Increment the number of atoms                         */
               (wpdb->natoms)++;
               
               /* Store the information read                            */
               CLEAR_PDB(p);
               p->atnum  = atnum;
               p->resnum = resnum;
               p->x      = (REAL)x;
               p->y      = (REAL)y;
               p->z      = (REAL)z;
               p->occ    = (REAL)occ;
               p->bval   = (REAL)bval;
               p->altpos = altpos;    /* 03.06.05 Added this one        */
               p->formal_charge  = charge;
               p->partial_charge = (REAL)charge;
               p->access = 0.0;
               p->radius = 0.0;
               p->atomType = NULL;
               p->next   = NULL;
               strcpy(p->record_type, record_type);
               strcpy(p->atnam,       atnam);
               strcpy(p->atnam_raw,   atnam_raw);
               strcpy(p->resnam,      resnam);
               strcpy(p->chain,       chain);
               strcpy(p->insert,      insert);
               strcpy(p->element,     element);
               strcpy(p->segid,       segid);
            }
            else   /* Partial occupancy                                 */
            {
               /* Set flag to say we've got a partial occupancy atom    */
               gPDBPartialOcc = TRUE;
               
               /* First in a group, store atom name                     */
               if(NPartial == 0)
               {
                  CurIns = insert[0];
                  CurRes = resnum;
                  strncpy(CurAtom,atnam,8);
               }
               
               if(strncmp(CurAtom,atnam,strlen(CurAtom)-1) || 
                  resnum != CurRes || 
                  CurIns != insert[0])
               {
                  /* Atom name has changed 
                     Select and store the OccRank highest occupancy atom
                  */
                  if(!blStoreOccRankAtom(OccRank,multi,NPartial,
                                         &wpdb->pdb,&p,&wpdb->natoms))
                  {
                     if(wpdb->pdb != NULL) FREELIST(wpdb->pdb, PDB);
                     wpdb->natoms = (-1);
                     if(cmd[0]) unlink(cmd);
                     return(NULL);
                  }
                  
                  /* Reset the partial atom counter                     */
                  NPartial = 0;
                  strncpy(CurAtom,atnam,8);
                  CurRes = resnum;
                  CurIns = insert[0];
               }
               
               if(NPartial < MAXPARTIAL)
               {
                  /* Store the partial atom data                        */
                  CLEAR_PDB((multi + NPartial));
                  multi[NPartial].atnum  = atnum;
                  multi[NPartial].resnum = resnum;
                  multi[NPartial].x      = (REAL)x;
                  multi[NPartial].y      = (REAL)y;
                  multi[NPartial].z      = (REAL)z;
                  multi[NPartial].occ    = (REAL)occ;
                  multi[NPartial].bval   = (REAL)bval;
                  multi[NPartial].formal_charge = charge;
                  multi[NPartial].partial_charge = (REAL)charge;
                  multi[NPartial].access = 0.0;
                  multi[NPartial].radius = 0.0;
                  multi[NPartial].atomType = NULL;
                  multi[NPartial].next   = NULL;
                  strcpy(multi[NPartial].record_type, record_type);
                  strcpy(multi[NPartial].atnam,       atnam);
                  /* 27.04.05 - added this line                         */
                  strcpy(multi[NPartial].atnam_raw,   atnam_raw);
                  strcpy(multi[NPartial].resnam,      resnam);
                  strcpy(multi[NPartial].chain,       chain);
                  strcpy(multi[NPartial].insert,      insert);
                  strcpy(multi[NPartial].element,     element);
                  /* 03.06.05 - added this line                         */
                  multi[NPartial].altpos = altpos;
                  /* 17.02.15 - added this line                         */
                  strcpy(multi[NPartial].segid,       segid);

                  NPartial++;
               }
            }
         }
         charge_buff[0] = '\0';
         charge = 0;
      }
   }

   if(NPartial != 0)
   {
      if(!blStoreOccRankAtom(OccRank,multi,NPartial,&wpdb->pdb,&p,
                             &wpdb->natoms))
      {
         if(wpdb->pdb != NULL) FREELIST(wpdb->pdb, PDB);
         wpdb->natoms = (-1);
         if(cmd[0]) unlink(cmd);
         return(NULL);
      }
   }

   if(cmd[0]) unlink(cmd);

   /* Return pointer to start of linked list                            */
   return(wpdb);
}

/************************************************************************/
/*>static BOOL blStoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                                  int NPartial, PDB **ppdb, PDB **pp, 
                                  int *natom)
   ------------------------------------------------------------------
*//**

   \param[in]     OccRank     Occupancy ranking required (>=1)
   \param[in]     multi[]     Array of PDB records for alternative atom
                             positions
   \param[in]     NPartial    Number of items in multi array
   \param[in,out] **ppdb      Start of PDB linked list (or NULL)
   \param[in,out] **pp        Current position in PDB linked list (or NULL)
   \param[in,out] *natom      Number of atoms read
   \return                     Memory allocation success

   Takes an array of PDB records which represent alternative atom 
   positions for an atom. Select the OccRank'th highest occupancy and
   add this one into the PDB linked list.

   To be called by doReadPDB().

-  17.03.94 Original    By: ACRM
-  08.10.99 Initialise IMaxOcc and MaxOcc
-  27.04.05 Added atnam_raw
-  03.06.05 Added altpos
-  14.10.05 Modified the flag value from 0.0 to -1.0 so that erroneous
            lower occupancies of 0.0 are read properly and written back
            with their occupancy (0.0) rather than the next higher
            occupancy. Handles residues like 1zeh/B16
-  04.08.14 Read charge and element. By: CTP
-  16.08.14 Read formal charge and set partial charge. By: CTP
-  17.02.15 Added segid support   By: ACRM
*/
static BOOL blStoreOccRankAtom(int OccRank, PDB multi[MAXPARTIAL], 
                               int NPartial, PDB **ppdb, PDB **pp, 
                               int *natom)
{
   int  i,
        j,
        IMaxOcc = 0;
   REAL MaxOcc  = (REAL)0.0,
        LastOcc = (REAL)0.0;
   
   if(OccRank < 1) OccRank = 1;
   
   for(i=0; i<OccRank; i++)
   {
      MaxOcc  = (REAL)0.0;
      IMaxOcc = 0;
      
      for(j=0; j<NPartial; j++)
      {
         if(multi[j].occ >= MaxOcc)
         {
            MaxOcc  = multi[j].occ;
            IMaxOcc = j;
         }
      }
      /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
         of zero are treated properly
      */
      multi[IMaxOcc].occ = (REAL)-1.0;

      /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
         of zero are treated properly
      */
      if(MaxOcc < (REAL)0.0) break;
      LastOcc = MaxOcc;
   }

   /* If we ran out of rankings, take the last one to be found          */
   /* 14.10.05 Changed flag value to -1 so that erroneous occupancies
      of zero are treated properly
   */
   if(MaxOcc < (REAL)0.0)
      MaxOcc = LastOcc;

   /* Store this atom
      Allocate space in the linked list
   */
   if(*ppdb == NULL)
   {
      INIT((*ppdb), PDB);
      *pp = *ppdb;
   }
   else
   {
      ALLOCNEXT(*pp, PDB);
   }
            
   /* Failed to allocate space; error return.                           */
   if(*pp==NULL)
      return(FALSE);
               
   /* Increment the number of atoms                                     */
   (*natom)++;
               
   /* Store the information read                                        */
   (*pp)->atnum  = multi[IMaxOcc].atnum;
   (*pp)->resnum = multi[IMaxOcc].resnum;
   (*pp)->x      = multi[IMaxOcc].x;
   (*pp)->y      = multi[IMaxOcc].y;
   (*pp)->z      = multi[IMaxOcc].z;
   (*pp)->occ    = MaxOcc;
   (*pp)->bval   = multi[IMaxOcc].bval;
   (*pp)->formal_charge = multi[IMaxOcc].formal_charge;
   (*pp)->partial_charge = multi[IMaxOcc].partial_charge;
   (*pp)->access = multi[IMaxOcc].access;
   (*pp)->radius = multi[IMaxOcc].radius;
   (*pp)->atomType = NULL;
   (*pp)->next   = NULL;
   /* 03.06.05 Added this line                                          */
   (*pp)->altpos = multi[IMaxOcc].altpos;
   strcpy((*pp)->record_type, multi[IMaxOcc].record_type);
   strcpy((*pp)->atnam,       multi[IMaxOcc].atnam);
   /* 27.04.05 Added this line                                          */
   strcpy((*pp)->atnam_raw,   multi[IMaxOcc].atnam_raw);
   strcpy((*pp)->resnam,      multi[IMaxOcc].resnam);
   strcpy((*pp)->chain,       multi[IMaxOcc].chain);
   strcpy((*pp)->insert,      multi[IMaxOcc].insert);
   strcpy((*pp)->element,     multi[IMaxOcc].element);
   /* 17.02.15 Added this line                                          */
   strcpy((*pp)->segid,       multi[IMaxOcc].segid);

   /* Patch the atom name to remove the alternate letter                */
   if(strlen((*pp)->atnam) > 4)
      ((*pp)->atnam)[4] = '\0';
   else
      ((*pp)->atnam)[3] = ' ';

   return(TRUE);
}

/************************************************************************/
/*>char *blFixAtomName(char *name, REAL occup)
   -------------------------------------------
*//**

   \param[in]     *name     Atom name read from file
   \param[in]     occup     Occupancy to allow fixing of partial occupancy
                            atom names
   \return                  Fixed atom name (pointer into name)

   Fixes an atom name by removing leading spaces, or moving a leading
   digit to the end of the string. Used by doReadPDB()

-  06.04.94 Original    By: ACRM
-  01.03.01 No longer static
-  03.06.05 The name passed in has always contained the column which is
            officially the alternate atom position indicator, but is 
            used by some programs as part of the atom name. Thus the
            properly constructed variable coming into the routine should
            be something like '1HG1 ' or '1HG1A' for an alternate atom
            position. However some programs use ' HG11'. Therefore we
            now check for a character in the last position and replace
            it with a space if there is a space in the preceeding
            position (e.g. ' CA A' -> ' CA  ') or if there is a 
            character in the first position (e.g. '1HG1A' -> '1HG1 ')
            or if the occupancy is not zero/one
            NOTE!!! To support this, the routine now has a second 
            parameter: REAL occup
-  07.07.14 Renamed to blFixAtomName() By: CTP
*/
char *blFixAtomName(char *name, REAL occup)
{
   char *newname;
   int  len;

   /* Default behaviour, just return the input string                   */
   newname = name;

   if(name[0] == ' ')       /* Name starts in column 14                 */
   {
      /* remove leading spaces                                          */
      KILLLEADSPACES(newname,name);
      /* 03.06.05 If the last-but-one position is a space, force the last
         position (the alternate atom indicator) to be a space
      */
      if(newname[2] == ' ')
      {
         newname[3] = ' ';
      }
   }
   else                     /* Name starts in column 13                 */
   {
      /* 03.06.05 The last character is the alternate atom indicator, 
         so force it to be a space
      */
      name[4] = ' ';
      
      /* If the first character is a digit, move it to the end          */
      if(isdigit(name[0]))
      {
         if((len = blChindex(name,' ')) == (-1))
         {
            /* We didn't find a space in the name, so add the character
               onto the end of the string and re-terminate
            */
            len         = strlen(name);
            newname     = name+1;
            name[len]   = name[0];
            name[len+1] = '\0';
         }
         else
         {
            /* We did find a space in the name, so put the first
               character there
            */
            newname     = name+1;
            name[len]   = name[0];
         }
      }
   }
   return(newname);
}

/************************************************************************/
/*>PDB *blRemoveAlternates(PDB *pdb)
   ---------------------------------
*//**

   \param[in,out] *pdb       PDB 
   \return                   Ammended linked list (in case start has
                             changed)

   Remove alternate atoms - we keep only the highest occupancy or the 
   first if there are more than one the same.

-  25.01.05 Original based on code written for Inpharmatica   By: ACRM
-  04.02.14 Use CHAINMATCH macro. By: CTP
-  07.07.14 Renamed to blRemoveAlternates() Use blWritePDBRecord()
            Use bl prefix for functions By: CTP

*/
PDB *blRemoveAlternates(PDB *pdb)
{
   PDB   *p, 
         *q, 
         *r, 
         *s, 
         *s_prev,
         *r_prev,
         *a_prev,
         *next,
         *alts[MAXPARTIAL];
   int   i, 
         altCount, 
         highest;
   
   
   /* Step through residues                                             */
   r_prev=NULL;
   for(p=pdb; p!=NULL; p=q)
   {
      q=blFindNextResidue(p);
      
      /* Step through atoms                                             */
      for(r=p; r!=q; NEXT(r))
      {
         if(r->altpos != ' ')
         {
#ifdef DEBUG
            fprintf(stderr,"\n\nAlt pos found for record:\n");
            blWritePDBRecord(stderr, r);
#endif
            /* We have an alternate, store it and search for the other 
               ones  
            */
            altCount=0;
            alts[altCount++] = r;
            /* Search through this residue for the alternates.
               This will work for 99.9% of files where the alternates are
               with the main atoms
            */
            for(s=r->next; s!=q; NEXT(s))
            {
               if(!strcmp(s->atnam_raw, alts[0]->atnam_raw))
               {
                  if(altCount < MAXPARTIAL)
                  {
                     alts[altCount++] = s;
#ifdef DEBUG
                     fprintf(stderr,"Partner atom found in res:\n");
                     blWritePDBRecord(stderr, s);
#endif
                  }
                  else
                  {
                     fprintf(stderr,"Warning==> More than %d alternative \
conformations in\n", MAXPARTIAL);
                     fprintf(stderr,"           residue %c%d%c atom %s. \
Increase MAXPARTIAL in ReadPDB.c\n", s->chain[0],
                                     s->resnum,
                                     s->insert[0],
                                     s->atnam);
                  }
               }
            }
            /* If we didn't find the alternates within the residue, then
               we search the rest of the records.
               This covers the known entry where the alternates are shoved
               on the end instead! 
            */
            if(altCount<2)
            {
#ifdef DEBUG
               fprintf(stderr,"No partner found in residue\n");
#endif

               s_prev = NULL;
               for(s=q; s!=NULL; NEXT(s))
               {
                  if((s->resnum    == alts[0]->resnum) &&
                     (s->insert[0] == alts[0]->insert[0]) &&
                     CHAINMATCH(s->chain,alts[0]->chain) &&
                     !strcmp(s->atnam_raw, alts[0]->atnam_raw))
                  {
                     if(altCount < MAXPARTIAL)
                     {
                        alts[altCount++] = s;
#ifdef DEBUG
                        fprintf(stderr,"Partner found outside \
residue:\n");
                        blWritePDBRecord(stderr, s);
#endif
                     }
                     else
                     {
                        fprintf(stderr,"Warning==> More than %d \
alternative conformations in\n", MAXPARTIAL);
                        fprintf(stderr,"           residue %c%d%c atom \
%s. Increase MAXPARTIAL in ReadPDB.c\n", s->chain[0],
                                         s->resnum,
                                         s->insert[0],
                                         s->atnam);

                        /* Move this record to the correct position in the
                           linked list

                           First unlink s from its old position
                        */
                        if(s_prev != NULL)
                           s_prev->next = s->next;
                     
                        /* Now link it back in where it should be       */
                        next = r->next;
                        r->next = s;
                        s->next = next;
                     }
                  }
                  s_prev = s;
               }
            }

            if(altCount < 2)
            {
#ifdef DEBUG
               fprintf(stderr,"No alternates found. Resetting ALT \
flag\n\n");
#endif
               alts[0]->altpos = ' ';
               
            }
            else
            {
               /* Find the highest occupancy, defaulting to the first   */
               highest = 0;
               for(i=0; i<altCount; i++)
               {
                  if(alts[i]->occ > alts[highest]->occ)
                     highest = i;
               }
               
               /* Delete the unwanted alternates                        */
               for(i=0; i<altCount; i++)
               {
                  if(i==highest) /* For the highest remove the ALT flag */
                  {
#ifdef DEBUG
                     fprintf(stderr,"Highest occupancy selected:\n");
                     blWritePDBRecord(stderr, alts[i]);
#endif
                     alts[i]->altpos = ' ';
                  }
                  else
                  {
                     /* If we are deleting the current record pointer, 
                        then we need to update it
                     */
                     if(alts[i] == r)
                     {
#ifdef DEBUG
                        fprintf(stderr,"Deleting current record \
pointer\n");
#endif
                        
                        if(r_prev == NULL)
                        {
                           r_prev = r;
                           NEXT(r);
                           /* We are deleting the head of the list so we 
                              must update the main list pointer
                           */
                           pdb = r;
                        }
                        else
                        {
                           r = r_prev;
                           FINDPREV(r_prev, pdb, r);
                        }
                     }
                     
                     /* Delete the alternate we don't need              */
#ifdef DEBUG
                     fprintf(stderr,"Deleting Alt pos record:\n");
                     blWritePDBRecord(stderr, alts[i]);
#endif
                     
                     FINDPREV(a_prev, pdb, alts[i]);
                     if(a_prev != NULL)
                        a_prev->next = alts[i]->next;
                     free(alts[i]);
                     
                  }  /* Not the highest, so we delete it                */
               }  /* Stepping through the alternates                    */
            }
            
         }  /* We have an alternate                                     */
         r_prev = r;
      }  /* Stepping through the atoms of this residue                  */
   }  /* Stepping through the residues                                  */
   return(pdb);
}



/************************************************************************/
/*>PDB *blDoReadPDBML(FILE *fp, int *natom, BOOL AllAtoms, int OccRank,
                      int ModelNum)
   --------------------------------------------------------------------
*//**

   \param[in]     *fpin    A pointer to type FILE in which the
                           .PDB file is stored.
   \param[in]     AllAtoms TRUE:  ATOM & HETATM records
                           FALSE: ATOM records only
   \param[in]     OccRank  Occupancy ranking
   \param[in]     ModelNum NMR Model number (0 = all)
   \param[out]    *natom   Number of atoms read. -1 if error.
   \return                 A pointer to the first allocated item of
                           the PDB linked list

   Reads a PDBML-formatted PDB file into a PDB linked list.
   
   The OccRank value indicates occupancy ranking  to read for partial 
   occupancy atoms. If any partial occupancy atoms are read the global 
   flag gPDBPartialOcc is set to TRUE.
   
   The global multiple-models flag is set to true if more than one model 
   is found.

-  22.04.14 Original By: CTP
-  02.06.14 Updated setting atnam_raw and parsing data from PDB atom site 
            labels (label_seq_id, etc.) if author-defined labels are 
            omitted. By: CTP
-  09.06.14 Set gPDBXML flag. By: CTP
-  07.07.14 Renamed to blDoReadPDBML() By: CTP
-  04.08.14 Read element and formal charge. By: CTP
-  15.08.14 Use CLEAR_PDB() to set default values. By: CTP
-  16.08.14 Read formal and partial charges. Use blCopyPDB() to copy data
            for partial occupancy atoms. By: CTP
-  18.08.14 Added XML_SUPPORT option. Return error if XML_SUPPORT not 
            defined By: CTP
-  26.08.14 Pad record_type to six characters. By: CTP
-  17.02.15 Added segid support   By: ACRM
-  25.02.15 Added some checks on potential NULL pointers
            Changed all strcpy()s to strncpy()s
            Initialized numeric content variable before each sscanf()
            in case it fails
            Calls blRenumAtomsPDB() at the end since we don't use the
            atom site IDs for atom numbers
-  13.03.15 Cosmetic changes
*/
PDB *blDoReadPDBML(FILE *fpin,
                   int  *natom,
                   BOOL AllAtoms,
                   int  OccRank,
                   int  ModelNum)
{
#ifndef XML_SUPPORT

   /* PDBML format not supported.                                       */
   *natom = (-1);
   return( NULL );

#else

   /* Parse PDBML-formatted file.                                       */
   xmlParserCtxtPtr ctxt;
   xmlDoc  *document;
   xmlNode *root_node  = NULL, 
           *sites_node = NULL, 
           *atom_node  = NULL, 
           *n          = NULL;
   int     size_t;
   char    xml_buffer[XML_BUFFER];
   xmlChar *content;
   double  content_lf;

   PDB     *pdb      = NULL,
           *curr_pdb = NULL,
           *end_pdb  = NULL,
           multi[MAXPARTIAL];

   int     NPartial       =  0,
           model_number   =  0;
   char    store_atnam[8] = "",
           pad_resnam[8]  = "";
       

   /* Zero natoms and reset flags                                       */
   gPDBXML        = TRUE;  /* global PDBML-fornmat flag                 */
   gPDBPartialOcc = FALSE; /* global partial occupancy flag             */
   gPDBMultiNMR   = FALSE; /* global multiple models flag               */
   *natom = 0;             /* atoms stored                              */

   /* Generate Document From Filehandle                                 */
   size_t = fread(xml_buffer, 1, XML_BUFFER, fpin);
   ctxt = xmlCreatePushParserCtxt(NULL, NULL, xml_buffer, size_t, "file");
   while ((size_t = fread(xml_buffer, 1, XML_BUFFER, fpin)) > 0) 
   {
      xmlParseChunk(ctxt, xml_buffer, size_t, 0);
   }
   xmlParseChunk(ctxt, xml_buffer, 0, 1);
   document = ctxt->myDoc;
   xmlFreeParserCtxt(ctxt);
   
   if(document == NULL)
   {
      /* Error: Failed to parse file                                    */
      *natom = -1;
      return(NULL);
   }
   

   /* Parse Document Tree                                               */
   root_node = xmlDocGetRootElement(document);   
   if(root_node == NULL)                   /* 25.02.15                  */
   {
      *natom = -1;
      return(NULL);
   }
      
   for(n=root_node->children; n!=NULL; NEXT(n))
   {
      /* Find Atom Sites Node                                           */
      if(!strcmp("atom_siteCategory",(char *) n->name))
      {
         /* Found Atom Sites                                            */
         sites_node = n;
         break;
      }
   }
   
   if(sites_node == NULL)
   {
      /* Error: Failed to find atom sites                               */
      xmlFreeDoc(document);
      *natom = -1;
      return(NULL);
   }


   /* Scan through atom nodes and populate PDB list.                    */
   for(atom_node = sites_node->children; atom_node; 
       atom_node = atom_node->next)
   {
      if(!strcmp("atom_site",(char *) atom_node->name))
      {
         /* Current PDB                                                 */
         INIT(curr_pdb,PDB);
         
         if(curr_pdb == NULL)
         {
            /* Error: Failed to store atom in pdb list                  */
            xmlFreeDoc(document);
            if(pdb != NULL) FREELIST(pdb,PDB);
            *natom = -1;
            return(NULL);
         }

         /* Set default values                                          */
         CLEAR_PDB(curr_pdb);
         strcpy(curr_pdb->chain,   "");
         strcpy(curr_pdb->atnam,   "");
         strcpy(curr_pdb->resnam,  "");
         strcpy(curr_pdb->insert, " ");
         strcpy(curr_pdb->element, "");
         strcpy(curr_pdb->segid,   "");

         /* Scan atom node children                                     */
         for(n=atom_node->children; n!=NULL; NEXT(n))
         {
            if(n->type != XML_ELEMENT_NODE){ continue; }
            content = xmlNodeGetContent(n);
            if(content == NULL)            /* 25.02.15                  */
            {
               FREELIST(pdb, PDB);
               *natom = -1;
               return(NULL);
            }
            
            /* Set PDB values                                           */
            if(!strcmp((char *)n->name, "B_iso_or_equiv"))
            {
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->bval = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "Cartn_x"))
            {
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->x = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "Cartn_y"))
            {
               sscanf((char *) content,"%lf",&content_lf);
               curr_pdb->y = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "Cartn_z"))
            {
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->z = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "auth_asym_id"))
            {
               strcpy(curr_pdb->chain, (char *)content);
            }
            else if(!strcmp((char *)n->name, "auth_atom_id"))
            {
               strcpy(curr_pdb->atnam, (char *)content);
            }
            else if(!strcmp((char *)n->name, "auth_comp_id"))
            {
               strcpy(curr_pdb->resnam, (char *) content);
            }
            else if(!strcmp((char *)n->name, "auth_seq_id"))
            {
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->resnum = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "pdbx_PDB_ins_code"))
            {
               /* set insertion code
                  25.02.15 Changed to strncpy()  By: ACRM
               */
               strncpy(curr_pdb->insert, (char *)content, 8);
            }
            else if(!strcmp((char *)n->name, "group_PDB"))
            {
               /* 25.02.15 Changed to strncpy()  By: ACRM               */
               strncpy(curr_pdb->record_type, (char *)content, 8);
               PADMINTERM(curr_pdb->record_type, 6);
            }
            else if(!strcmp((char *)n->name, "occupancy"))
            {
               content_lf = (REAL)0.0;     /* 25.02.15                  */
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->occ = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "label_alt_id"))
            {
               /* Use strlen as test for alt position                   */
               curr_pdb->altpos = strlen((char *)content) ? content[0]:' ';
            }
            else if(!strcmp((char *)n->name, "pdbx_PDB_model_num"))
            {
               content_lf = (REAL)0.0;     /* 25.02.15                  */
               sscanf((char *)content, "%lf", &content_lf);
               model_number = (int)content_lf;
            }
            else if(!strcmp((char *)n->name, "type_symbol"))
            {
               /* 25.02.15 Changed to strncpy()  By: ACRM               */
               strncpy(curr_pdb->element, (char *)content, 8);
            }
            else if(!strcmp((char *)n->name, "label_asym_id"))
            {
               if(strlen(curr_pdb->chain) == 0)
               {
                  /* 25.02.15 Changed to strncpy()  By: ACRM            */
                  strncpy(curr_pdb->chain, (char *)content, 8);
               }
            }
            else if(!strcmp((char *)n->name, "label_atom_id"))
            {
               if(strlen(curr_pdb->atnam) == 0)
               {
                  /* 25.02.15 Changed to strncpy()  By: ACRM            */
                  strncpy(curr_pdb->atnam, (char *)content, 8);
               }
            }
            else if(!strcmp((char *) n->name, "label_comp_id"))
            {
               if(strlen(curr_pdb->resnam) == 0)
               {
                  /* 25.02.15 Changed to strncpy()  By: ACRM            */
                  strncpy(curr_pdb->resnam, (char *)content, 8);
               }
            }
            else if(!strcmp((char *)n->name, "label_seq_id"))
            {
               if((curr_pdb->resnum == 0) && 
                  (strlen((char *)content) > 0))
               {
                  content_lf = (REAL)0.0;  /* 25.02.15                  */
                  sscanf((char *)content, "%lf", &content_lf);
                  curr_pdb->resnum = (REAL)content_lf;
               }
            }
            else if(!strcmp((char *)n->name, "pdbx_formal_charge"))
            {
               content_lf = (REAL)0.0;     /* 25.02.15                  */
               sscanf((char *)content, "%lf", &content_lf);
               curr_pdb->formal_charge = (int)content_lf;
               curr_pdb->partial_charge = (REAL)content_lf;
            }
            else if(!strcmp((char *)n->name, "seg_id"))  /* 17.02.15    */
            {
               if(strlen(curr_pdb->segid) == 0)
               {
                  /* 25.02.15 Changed to strncpy()  By: ACRM            */
                  strncpy(curr_pdb->segid, (char *)content, 8);
               }
            }

            xmlFree(content);           
         }
         
         /* Set raw atom name
            Note: The text pdb format uses columns 13-16 to store the atom
                  name. By convention, columns 13-14 contain the 
                  right-justified element symbol for the atom.
                  
                  The raw atom name is equivalent to colums 13-16 of a
                  pdb-formatted text file .                             
         */

         if(strlen(curr_pdb->atnam) == 1)
         {
            /* copy 1-letter name atnam_raw                             */
            strcpy((curr_pdb->atnam_raw), " ");
            /* 25.02.15 Changed to strncpy()  By: ACRM                  */
            strncpy((curr_pdb->atnam_raw)+1, curr_pdb->atnam, 7);
         }
         if(strlen(curr_pdb->atnam) == 4)
         {
            /* copy 4-letter name atnam_raw                             */
            /* 25.02.15 Changed to strncpy()  By: ACRM                  */
            strncpy(curr_pdb->atnam_raw, curr_pdb->atnam, 8);
         }
         else if(strlen(curr_pdb->element) == 1)
         {
            strcpy((curr_pdb->atnam_raw),               " ");
            /* 25.02.15 Changed to strncpy()  By: ACRM                  */
            strncpy((curr_pdb->atnam_raw)+1, curr_pdb->atnam, 7);
         }
         else
         {
            /* 25.02.15 Changed to strncpy()  By: ACRM                  */
            strncpy(curr_pdb->atnam_raw, curr_pdb->atnam, 4);
         }
         
         /* Pad atom names to 4 characters                              */
         PADMINTERM(curr_pdb->atnam,     4);
         PADMINTERM(curr_pdb->atnam_raw, 4);
         
         /* Pad Residue Name
            Note: The text pdb format uses columns 18-20 to store the 
                  residue name (right-justified).
                  
                  curr_pdb->resnam is is equivalent to colums 18-21 of a
                  pdb-formatted text file.                              
         */
         sprintf(pad_resnam, "%3s", curr_pdb->resnam);
         PADMINTERM(pad_resnam, 4);
         /* 25.02.15 Changed to strncpy()  By: ACRM                     */
         strncpy(curr_pdb->resnam, pad_resnam, 8);
         
         /* Set chain to " " if not already set                         */
         if(strlen(curr_pdb->chain) == 0)
         {
            strcpy(curr_pdb->chain, " ");
         }

         /* Pad the segment id                                          */
         PADMINTERM(curr_pdb->segid, 4);

         /* Set multi-model flag                                        */
         if(model_number > 1)
         {
            gPDBMultiNMR = TRUE;
         }

         /* Filter: Model Number                                        */
         if(model_number != ModelNum)
         {
            /* Free curr_pdb                                            */
            FREELIST(curr_pdb,PDB);
            curr_pdb = NULL;
            
            if(model_number > ModelNum)
            {
               break;    /* skip rest of tree                           */
            }
            else
            {
               continue; /* filter                                      */
            }
         }


         /* Filter: All Atoms                                           */
         if(!AllAtoms && strncmp(curr_pdb->record_type, "ATOM  ", 6))
         {
            /* Free curr_pdb and skip atom                              */
            FREELIST(curr_pdb,PDB);
            curr_pdb = NULL;
            continue;    /* filter                                      */
         }


         /* Add partial occ atom from temp storage to output PDB list   */
         if((NPartial != 0) && strcmp(curr_pdb->atnam,store_atnam))
         {
            /* Store atom                                               */
            if(blStoreOccRankAtom(OccRank,multi,NPartial,&pdb,&end_pdb,
                                  natom))
            {
               LAST(end_pdb);
               NPartial = 0;
            }
            else
            {
               /* Error: Failed to store partial occ atom               */
               xmlFreeDoc(document);
               if(curr_pdb != NULL) FREELIST(curr_pdb,PDB);
               if(pdb      != NULL) FREELIST(pdb,     PDB);
               *natom = -1;
               return(NULL);
            }
         }


         /* Set atom number
            Note: Cannot use atom site id for atom number so base atnum on
                  number of atoms stored 

            25.02.15 We will renumber afterwards
         */
         curr_pdb->atnum = *natom + 1;
         
         
         /* Add partial occupancy atom to temp storage                  */
         if((curr_pdb->altpos != ' ') && (NPartial < MAXPARTIAL))
         {
            /* Copy the partial atom data to storage                    */
            blCopyPDB(&multi[NPartial], curr_pdb);

            /* Set global partial occupancy flag                        */
            gPDBPartialOcc = TRUE;
            
            /* Store current atom name                                  */
            /* 25.02.15 Changed to strncpy()  By: ACRM                  */
            strncpy(store_atnam, curr_pdb->atnam, 8);
            NPartial++;

            /* Free curr_pdb and continue                               */
            FREELIST(curr_pdb,PDB);
            curr_pdb = NULL;
            continue;
         }

         
         /* Store Atom                                                  */
         if(pdb == NULL)
         {
            pdb      = curr_pdb;
            end_pdb  = curr_pdb;
            curr_pdb = NULL;
            *natom   = 1;
         }
         else
         {
            end_pdb->next = curr_pdb;
            end_pdb       = curr_pdb;
            curr_pdb      = NULL;
            (*natom)++;
         }
      }
   }
      
   /* Free document                                                     */
   xmlFreeDoc(document);

   /* Store final atom (if partial occupancy)                           */
   if(NPartial != 0)
   {
      if(!blStoreOccRankAtom(OccRank,multi,NPartial,&pdb,&end_pdb,natom))
      {
         /* Error: Failed to store atom in pdb list                     */
         if(pdb != NULL) FREELIST(pdb,PDB);
         *natom = -1;
         return(NULL);
      }
   }
   
   /* Check atoms have been stored                                      */
   if(pdb == NULL || *natom == 0)
   {
      /* Error: pdb list empty or no atoms stored                       */
      if(pdb != NULL) FREELIST(pdb,PDB);
      *natom = -1;
   }

   /* 25.02.15 Renumber atoms since we don't use the atom site IDs      */
   blRenumAtomsPDB(pdb, 1);
    
   /* Return PDB linked list                                            */
   return(pdb);

#endif
}


/************************************************************************/
/*>BOOL blCheckFileFormatPDBML(FILE *fp)
   -------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE.
   \return                 File is in PDBML format?

   Simple test to detect PDBML-formatted pdb file.
   
   Todo: Consider replacement with general function to detect file format
         for uncompressed file returning file type (eg pdb/pdbml/unknown).
   

-  22.04.14 Original By: CTP
-  07.07.14 Renamed to blCheckFileFormatPDBML() By: CTP
-  29.08.14 Function re-written to take sample from the input stream then
            reset the stream with ungetc. By: CTP
-  31.08.14 Bugfix: Check for 'PDBx:datablock' tag skipped if blank line 
            before xml tag. By: CTP
-  09.09.14 Use rewind() for DOS instead of pushing sample back on stream 
            with ungetc(). By: CTP
-  29.09.14 Use single character check for pdbml files for Windows or 
            systems where ungetc() fails after pushback of singe char. 
            By: CTP

*/
BOOL blCheckFileFormatPDBML(FILE *fp)
{
#if !defined(SINGLE_CHAR_FILECHECK) && !defined(MS_WINDOWS)

   /* Default Filetype Check                                            */
   char buffer[XML_SAMPLE];
   int  i, c;
   BOOL found_xml  = FALSE,
        found_pdbx = FALSE;

   /* store sample from stream                                          */
   for(i = 0; i < (XML_SAMPLE - 1); i++)
   {
      c = fgetc(fp);
      if(c == EOF || feof(fp)) break;
      buffer[i] = (char)c;      
   }
   buffer[i] = '\0'; /* terminate string                                */

   /* push sample back on input stream                                  */
   for(i = strlen(buffer) - 1; i >= 0; i--)
   {
      ungetc(buffer[i], fp);
   }

   /* check first line                                                  */
   if(!strncmp(buffer,"<?xml ",6)) found_xml  = TRUE;
   
   /* check remaining lines                                             */
   for(i = 0; i < strlen(buffer); i++)
   {
      if(buffer[i] != '\n') continue;

      /*i++;*/
      if(!strncmp(&buffer[i+1],"<?xml ",6))            found_xml  = TRUE;
      if(!strncmp(&buffer[i+1],"<PDBx:datablock ",16)) found_pdbx = TRUE;
   }

   return ((found_xml && found_pdbx) ? TRUE : FALSE);

#else

   /* Single Character Filetype Check                                   */
   int c;

   /* get single char from input stream                                 */
   c = fgetc(fp);
   if(c == EOF || feof(fp)) return FALSE;

   /* pushback character                                                */
   ungetc(c, fp);

   /* detect filetype                                                   */
   return (((char)c == '<') ? TRUE:FALSE);

#endif 
}

/************************************************************************/
/*>static void blProcessElementField(char *element_field, char *element)
   ---------------------------------------------------------------------
*//**

   \param[in]  element_field   Columns 77 to 78 of the ATOM/HETATM record 
                               of pdb file.
   \param[out] element         Element symbol.
   
   Get element symbol for ATOM/HETATM record.
   
-  04.08.14 Original By: CTP
   
*/
static void blProcessElementField(char *element, char *element_field)
{
   char element_sym[4] = "";
   char *element_ptr   = NULL;

   /* Get element                                                       */
   if(strlen(element_field) >= 2)
   {
      strncpy(element_sym, element_field, 2);
      element_sym[2] = '\0';
      KILLLEADSPACES(element_ptr, element_sym);
      strcpy(element,element_ptr);
   }

   return;
}

/************************************************************************/
/*>static void blProcessChargeField(char *element_charge, int *charge)
   -------------------------------------------------------------------
*//**

   \param[in]  element_charge  Columns 79 to 80 of the ATOM/HETATM record 
                               of pdb file.
   \param[out] charge          Charge on the atom.
   
   Get formal charge for ATOM/HETATM record.

-  04.08.14 Original By: CTP
   
*/
static void blProcessChargeField(int *charge, char *charge_field)
{
   /* Get charge magnitude                                              */
   if(strlen(charge_field) >= 2 && isdigit(charge_field[0]))
   {
      *charge = charge_field[0] - '0';
   }

   /* Get charge sign                                                   */
   if(strlen(charge_field) >= 2 && charge_field[1] == '-')
   {
      *charge = *charge * -1;
   }

   return;
}

/************************************************************************/
/*>void blFreeWholePDB(WHOLEPDB *wpdb)
   ---------------------------------
*//**

   \param[in]     *wpdb    WHOLEPDB structure to be freed

   Frees the header, trailer and atom content from a WHOLEPDB structure

-  30.05.02  Original   By: ACRM
-  07.07.14  Renamed to blFreeWholePDB() By: CTP
*/
void blFreeWholePDB(WHOLEPDB *wpdb)
{
   blFreeStringList(wpdb->header);
   blFreeStringList(wpdb->trailer);
   FREELIST(wpdb->pdb, PDB);
   free(wpdb);
}


/************************************************************************/
/*>WHOLEPDB *blReadWholePDB(FILE *fpin)
   ------------------------------------
*//**

   \param[in]     *fpin     File pointer
   \return                  Whole PDB structure containing linked
                            list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

-  07.03.07 Made into a wrapper to doReadWholePDB()
-  07.07.14 Use blDoReadWholePDB() Renamed to blReadWholePDB() By: CTP
*/
WHOLEPDB *blReadWholePDB(FILE *fpin)
{
   WHOLEPDB *wpdb;
   wpdb = blDoReadPDB(fpin, TRUE, 1, 1, TRUE);
   wpdb->pdb = blRemoveAlternates(wpdb->pdb);
   return(wpdb);
}

/************************************************************************/
/*>WHOLEPDB *blReadWholePDBAtoms(FILE *fpin)
   -----------------------------------------
*//**

   \param[in]     *fpin     File pointer
   \return                  Whole PDB structure containing linked
                            list to PDB coordinate data

   Reads a PDB file, storing the header and trailer information as
   well as the coordinate data. Can read gzipped files as well as
   uncompressed files.

   Coordinate data is accessed as linked list of type PDB as follows:
   
   WHOLEPDB *wpdb;
   PDB      *p;
   wpdb = ReadWholePDB(fp);
   for(p=wpdb->pdb; p!=NULL; p=p->next)
   {
      ... Do something with p ...
   }

-  07.03.07 Made into a wrapper to doReadWholePDB()
-  07.07.14 Use blDoReadWholePDB() Renamed to blReadWholePDBAtoms() 
            By: CTP
*/
WHOLEPDB *blReadWholePDBAtoms(FILE *fpin)
{
   WHOLEPDB *wpdb;
   wpdb = blDoReadPDB(fpin, FALSE, 1, 1, TRUE);
   wpdb->pdb = blRemoveAlternates(wpdb->pdb);
   return(wpdb);
}


/************************************************************************/
/*>static void StoreConectRecords(WHOLEPDB *wpdb, char *buffer)
   ------------------------------------------------------------
*//**
   \param[in,out]    *wpdb      Whole PDB structure
   \param[in]        *buffer    A line containing a CONECT record

   Stores the connectivity data from a CONECT record into the PDB linked
   list

-  18.02.15  Original   By: ACRM
-  05.03.15  Initialize all connected atoms to NULL and check for NULLs
             when storing the CONECTs. This deals with cases where an
             atom specified in a CONECT record has alternate occupancies
             and an alternate position is removed. Previously this led
             to core dumps. e.g. PDB code 3pnw
*/
static void StoreConectRecords(WHOLEPDB *wpdb, char *buffer)
{
   char record_type[8];
   int  i, j,
        nConect,
        atoms[5];
   PDB  *p,
        *atomsP[5];
   BOOL gotLink;

   fsscanf(buffer,"%6s%5d%5d%5d%5d%5d", 
           record_type,&atoms[0],&atoms[1],&atoms[2],&atoms[3],&atoms[4]);

   /* Find the PDB pointers for the (up to) 5 CONECT atoms              */
   nConect = 0;
   for(i=0; i<5; i++)
   {
      /* Have we run out of specified atoms?                            */
      if(atoms[i] == 0)
         break;

      atomsP[i] = NULL;                                 /* 05.03.15     */

      /* Look for this atom                                             */
      for(p=wpdb->pdb; p!=NULL; NEXT(p))
      {
         if(atoms[i] == p->atnum)
         {
            /* Found it so store and break out                          */
            atomsP[i] = p;
            nConect++;
            break;
         }
      }
   }

   /* Set the connections from atom 0                                   */
   if(atomsP[0] != NULL)                                /* 05.03.15     */
   {
      for(i=1; i<nConect; i++)
      {
         /* Look to see if we have this conect stored already           */
         gotLink = FALSE;
         for(j=0; j<atomsP[0]->nConect; j++)
         {
            if(atomsP[i] != NULL)                       /* 05.03.15     */
            {
               if(atomsP[0]->conect[j] == atomsP[i])
               {
                  gotLink = TRUE;
                  break;
               }
            }
         }
         
         /* If not then store it                                        */
         if(!gotLink && (atomsP[0]->nConect < MAXCONECT))
         {
            if(atomsP[i] != NULL)                       /* 05.03.15     */
            {
               atomsP[0]->conect[atomsP[0]->nConect] = atomsP[i];
               (atomsP[0]->nConect)++;
            }
         }
      }

      /* Set the connections in the other direction                     */
      for(i=1; i<nConect; i++)
      {
         if(atomsP[i] != NULL)                          /* 05.03.15     */
         {
            /* Look to see if we have this conect stored already        */
            gotLink = FALSE;
            for(j=0; j<atomsP[i]->nConect; j++)
            {
               if(atomsP[i]->conect[j] == atomsP[0])
               {
                  gotLink = TRUE;
                  break;
               }
            }
            
            /* If not then store it                                     */
            if(!gotLink && (atomsP[i]->nConect < MAXCONECT))
            {
               atomsP[i]->conect[atomsP[i]->nConect] = atomsP[0];
               (atomsP[i]->nConect)++;
            }
         }
      }
   }
}

#ifdef XML_SUPPORT
/************************************************************************/
/*>static STRINGLIST *blParseHeaderPDBML(FILE *fpin)
   -------------------------------------------------
*//**

   \param[in]     *fpin   File pointer
   \return                STRINGLIST with basic header information.

   Parses a PDBML file and creates HEADER and TITLE lines.

-  22.04.14 Original. By: CTP
-  07.07.14 Renamed to blParseHeaderPDBML() By: CTP
-  18.08.14 Return NULL if XML not supported. By: CTP
-  10.09.14 Use blSetPDBDateField() to set date field. By: CTP

*/
static STRINGLIST *blParseHeaderPDBML(FILE *fpin)
{
#ifndef XML_SUPPORT

   /* PDBML format not supported.                                       */
   return NULL;

#else

   /* Parse PDBML header                                                */
   xmlParserCtxtPtr ctxt;
   xmlDoc  *document;
   xmlNode *root_node = NULL, 
           *node      = NULL,
           *subnode   = NULL,
           *n         = NULL;
   int     size_t;
   char    xml_buffer[XML_BUFFER];
   xmlChar *content, *attribute;
   
   STRINGLIST *wpdb_header = NULL,
              *title_lines = NULL;

   char header_line[82]  = "",
        title_line[82]   = "",
        header_field[41] = "",
        pdb_field[5]     = "",
        date_field[10]   = "",
        title_field[71]  = "";
        
   int cut_from = 0, 
       cut_to   = 0,
       nlines   = 0,
       i        = 0;

   
   /* Generate Document From Filehandle                                 */
   size_t = fread(xml_buffer, 1, XML_BUFFER, fpin);
   ctxt = xmlCreatePushParserCtxt(NULL, NULL, xml_buffer, size_t, "file");
   while ((size_t = fread(xml_buffer, 1, XML_BUFFER, fpin)) > 0) 
   {
      xmlParseChunk(ctxt, xml_buffer, size_t, 0);
   }
   xmlParseChunk(ctxt, xml_buffer, 0, 1);
   document = ctxt->myDoc;
   xmlFreeParserCtxt(ctxt);

   if(document == NULL){ return(NULL); } /*        failed to parse file */


   /* Parse Document Tree                                               */
   root_node = xmlDocGetRootElement(document);
   for(node = root_node->children; node; node = node->next)
   {
      if(node->type != XML_ELEMENT_NODE){ continue; }

      /* get header                                                     */
      if(!strcmp("struct_keywordsCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("pdbx_keywords",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);
               strncpy(header_field,(char *) content,40);
               xmlFree(content);
            }
         }
      }

      /* get date                                                       */
      if(!strcmp("database_PDB_revCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("date_original",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);
               
               /* convert date format                                   */
               blSetPDBDateField(date_field, (char *)content);
               
               xmlFree(content);
            }
         }
      }

      /* get pdb code                                                   */
      if(!strcmp("entryCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            if(strcmp("entry",(char *) subnode->name)){ continue; }
            attribute = xmlGetProp(subnode,(xmlChar *) "id");
            strncpy(pdb_field,(char *)attribute,4);
            pdb_field[4] = '\0';
            xmlFree(attribute);
         }
      }

      /* get title                                                      */
      if(!strcmp("structCategory",(char *)node->name))
      {
         for(subnode = node->children; subnode; subnode = subnode->next)
         {
            for(n=subnode->children; n; n = n->next)
            {
               if(strcmp("title",(char *) n->name)){ continue; }
               content = xmlNodeGetContent(n);

               /* Get title lines as STRINGLIST                         */
               for(i=0; i<strlen((char *) content); i++)
               {
                  if(content[i] == ' ' ) cut_to = i;
                  if(i == strlen((char *) content) - 1) cut_to = i+1;

                  /* split and store title line                         */
                  if( (i && !((i - cut_from)%70)) || 
                      i == strlen((char *) content)-1 )
                  {
                     nlines++;
                     cut_to = (cut_from == cut_to) ? i : cut_to;
                     strncpy(title_field,
                             (char *) &content[cut_from],
                             cut_to - cut_from);
                     title_field[cut_to - cut_from] = '\0';
                     PADMINTERM(title_field,70);
                     cut_from = cut_to;
                     i        = cut_to;

                     if(nlines == 1)
                     {
                        sprintf(title_line, "TITLE     %s\n",
                                title_field);
                        title_lines = blStoreString(NULL,title_line);
                     }
                     else
                     {
                        sprintf(title_line, "TITLE   %2d%s\n", nlines,
                                title_field);
                        blStoreString(title_lines,title_line);
                     }
                  }
               }
               xmlFree(content);
            }
         }
      }
   }

   /* Free document                                                     */
   xmlFreeDoc(document);

   /* Cleanup xml parser                                                */
   xmlCleanupParser();

   /* Create Header Line                                                */
   if(!strlen(header_field))
   {
      strcpy(header_field,"Converted from PDBML");
   }
   sprintf(header_line, "HEADER    %-40s%9s   %4s              \n",
           header_field, date_field, pdb_field);
   
   /* Make Stringlist                                                   */
   wpdb_header = blStoreString(wpdb_header, header_line);
   wpdb_header->next = title_lines;
   
   return(wpdb_header);

#endif
}

/************************************************************************/
/*>static BOOL blSetPDBDateField(char *pdb_date, char *pdbml_date)
   ---------------------------------------------------------------
*//**

   \param[out]    *pdb_date      PDB date string   'dd-MTH-yy'
   \param[in]     *pdbml_date    PDBML date string 'yyyy-mm-dd'
   \return                       Success?

   Convert pdbml date format to pdb date format.

-  10.09.14 Original. By: CTP

*/
static BOOL blSetPDBDateField(char *pdb_date, char *pdbml_date)
{
   char month_letter[12][4] = {"JAN","FEB","MAR","APR","MAY","JUN",
                               "JUL","AUG","SEP","OCT","NOV","DEC"};
   int day   = 0,
       month = 0,
       year  = 0,
       items = 0;
   
   /* parse pdbml date                                                  */
   items = sscanf(pdbml_date, "%4d-%2d-%2d", &year, &month, &day);

   /* error check                                                       */
   if(items != 3 || 
      year == 0 || month == 0 || day == 0 || 
      day   < 1 || day > 31   ||
      month < 1 || month > 12 ||
      year  < 1900)
   {
      /* conversion failed                                              */
      strncpy(pdb_date, "         ", 10);
      return FALSE;
   }
   
   /* set pdb date                                                      */
   sprintf(pdb_date, "%02d-%3s-%02d",
           day, month_letter[month - 1], year % 100);

   return TRUE;
}
#endif


/*************************************************************************

   Program:    ProFit
   File:       main.c
   
   Version:    V3.1
   Date:       31.03.09
   Function:   Protein Fitting program. Main routines.
   
   Copyright:  SciTech Software / UCL 1992-2009
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain.

   It may not be copied or made available to third parties, but may be
   freely used by non-profit-making organisations who have obtained it
   directly from the author or by FTP.

   You are requested to send EMail to the author to say that you are 
   using this code so that you may be informed of future updates.

   The code may not be made available on other FTP sites without express
   permission from the author.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If
   someone else breaks this code, the author doesn't want to be blamed
   for code that does not work! You may not distribute any
   modifications, but are encouraged to send them to the author so
   that they may be incorporated into future versions of the code.

   Such modifications become the property of Dr. Andrew C.R. Martin and
   SciTech Software though their origin will be acknowledged.

   The code may not be sold commercially or used for commercial purposes
   without prior permission from the author.
   
**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V0.1  25.09.92  Original
   V0.2  02.10.92  Added CURSES support
   V0.3  07.10.92  Added Amiga windows and paging support
   V0.4  09.10.92  Added N&W alignment support & fixed bug in multi-zones
   V0.5  08.10.93  Various tidying for Unix & chaned for booklib 
   V0.6  05.01.94  Reads HELPDIR and DATADIR environment variables under
                   Unix
   V0.7  24.11.94  Uses ReadPDBAtoms()
   V0.8  17.07.95  Stripped out the windowing support (never used!). Will
                   add a Tcl/Tk interface later!
                   Multiple chains now work correctly.
   V1.0  18.07.95  Insert codes now work.
                   First official release (at last!).
   V1.1  20.07.95  Added WEIGHT command support and translation vector
                   output from MATRIX command
   V1.1a 22.07.95  Stop crash on writing when not fitted
   V1.2  25.07.95  Added GAPPEN command
                   Added chain label printing in Status
                   Added optional filename parameter to RESIDUE command
   V1.3  31.07.95  Fixed bug in fitting.c where end of zone=end of chain
                   other than the last chain
   V1.4  14.08.95  Fixed bug in fitting.c with RESIDUE command was not
                   printing RMS for last residue
   V1.5  21.08.95  Fixed bug in NWAlign.c (when last zone not at end of
                   chain) and Bioplib library bug in align()
   V1.5a 02.10.95  Added printing of CofG information with MATRIX command
   V1.5b 15.11.95  NWAlign.c: Prints normalised score
   V1.6  20.11.95  Added READALIGNMENT command and code
   V1.6a 21.11.95  Fixed a couple of warnings under gcc
   V1.6b 22.11.95  Modified code in SetNWZones() such that deletions at
                   the same position in both sequences don't cause a 
                   problem.
   V1.6c 13.12.95  Further fix to double deletions. Added info when
                   zones mismatch in fitting.c
   V1.6d 24.01.96  Fixed bug in status command when printing atom names
                   containing spaces.
   V1.6e 31.05.96  Added BVAL command and code
   V1.6f 13.06.96  Added BWEIGHT command and code
   V1.6g 18.06.96  Replaced MODE_* with ZONE_MODE_* and use FindZonePDB()
                   from bioplib rather than own version
   V1.7  23.07.96  Supports atom wildcards. Some comment tidying.
   V1.7b 11.11.96  Added REF keyword to end of BVAL command allowing
                   only the reference structure to be considered
   V1.7c 18.11.96  Added IGNOREMISSING option
   V1.7d 20.12.96  Added NFITTED command
   V1.7e 27.06.97  Allows WRITE and RESIDUE to output to a pipe
   V1.7f 03.07.97  Added break into CreateFitArrays() to fix core dump
                   on bad multiple-occupancy PDB files
   V1.7g 06.05.98  Rewrite of NWAlign/SetNWZones()
   V1.8  07.05.98  Skipped for release
   V2.0  01.03.01  Now supports multiple structure fitting and iterative
                   zone updating
   V2.1  28.03.01  Parameter for ITERATE and added CENTRE command
   V2.2  20.02.01  Fixed some Bioplib problems related to raw atom names
                   and multiple occupancies
   V2.3  01.12.04  Skipped for release
   V2.4  03.06.05  Some significant changes in Bioplib for better handling
                   of raw atom names with multiple occupancies
   V2.5  07.06.05  Another bug in bioplib
   V2.5.1 10.06.05 Bug in bioplib PDB2Seq() with CA-only chains
   V2.5.2 14.10.05 Modified doReadPDB() to cope with corrupt partial
                   occupancies like 1zeh/B13
   V2.5.2 06.07.06 Fixed raw atom name when ILE-CD changed to ILE-CD1
   V2.5.4 29.06.07 Bioplib change for Mac compilation
   V-.-.- 14.03.08 Added DELZONE command.
   V-.-.- 17.03.08 Added DELRZONE command.
   V-.-.- 19.03.08 Added SETCENTRE/SETCENTER command.
   V-.-.- 27.03.08 Added SCRIPT command.
   V-.-.- 02.04.08 Added HEADER command.
   V-.-.- 03.04.08 Added PAIRDIST command.
   V-.-.- 07.04.08 Added DISTCUTOFF command.
   V-.-.- 22.04.08 Added handling of lowercase chain and inserts.
   V-.-.- 29.04.08 Added calls to ReadWholePDB() for reading PDB structures.
   V-.-.- 01.05.08 Added support for reading low occupancy atoms.
   V-.-.- 02.05.08 Removed InitPDBCoordRecord() and IsPDBCoordRecord().
                   Removed ReadPDBHeader(), WritePDBHeader() and 
                   WritePDBFooter().
   V2.5.9 22.05.08 Cleaned-up comments. Added version number. Removed unused 
                   functions, ShowOverlap() and ShowAllOverlap().
   V2.6   12.06.08 Modified ReadStructure() to use updated read functions in 
                   bioplib.
   V2.6   11.08.08 ShowRMS() turned-off for SetRMSAtoms() and SetRMSZone()
                   when in quiet mode. Added MULTIREF command.
   V2.6   16.10.08 Added ALLVSALL, SETREF, ORDERFIT and TRIMZONES commands.
   V2.6   20.10.08 Cleaned-up code.
   V2.6   23.10.08 Added gWtAverage flag to allow old use of old 
                   weighting system.
   V2.6   27.10.08 Added NOFIT command.
   V3.0   06.11.08 Fixed bug with SETREF.
   V3.0   07.11.08 Added option for STATUS to output to a file.
   V3.0   14.11.08 Added GNU Readline Library support.
   V3.0   11.02.09 Changed DELZONE and DELRZONE commands to use ALL.
   V3.0   16.02.09 Removed Invert3x3Matrix() and Transpose3x3Matrix().
   V3.0   13.03.09 Fixed bug with ATOMS and RATOMS commands. Set 
                   ReadStructure() to clear fitted coordinates when 
                   loading a new mobile structure.
   V3.1   31.03.09 Updated version number.

*************************************************************************/
#define MAIN 
#include "ProFit.h"


/************************************************************************/
/*>void logo(void)
   ---------------
   Displays the ProFit logo.
   
   25.09.92 Original
   01.10.92 Added version
   10.10.93 Version 0.5
   05.01.94 Version 0.6
   24.11.94 Version 0.7
   17.07.95 Version 0.8; removed screen() version
   18.07.95 Version 1.0
   20.07.95 Version 1.1
   22.07.95 Version 1.2
   31.07.95 Version 1.3
   14.08.95 Version 1.4
   21.08.95 Version 1.5
   20.11.95 Version 1.6
   23.07.96 Version 1.7
   07.05.98 Version 1.8
   15.02.01 Version 2.0
   20.03.01 Version 2.1
   20.12.01 Version 2.2
   01.12.04 Version 2.3
   03.06.05 Version 2.4
   07.06.05 Version 2.5
   10.06.05 Version 2.5.1
   14.10.05 Version 2.5.2
   06.07.06 Version 2.5.3
   29.06.07 Version 2.5.4
   28.03.08 Version Working_Copy By: CTP
   22.05.08 Version 2.6
   04.06.08 Version Working Copy
   03.11.08 Version 3.0
   31.03.09 Version 3.1
*/
void logo(void)
{
   printf("\n                  PPPPP                FFFFFF ii   tt\n");
   printf("                  PP  PP               FF          tt\n");
   printf("                  PP  PP rrrrr   oooo  FF     ii  ttttt\n");
   printf("                  PPPPP  rr  rr oo  oo FFFF   ii   tt\n");
   printf("                  PP     rr     oo  oo FF     ii   tt\n");
   printf("                  PP     rr     oo  oo FF     ii   tt\n");
   printf("                  PP     rr      oooo  FF      ii   ttt\n\n");
   printf("                      Protein Least Squares Fitting\n\n");
   printf("                               Version 3.1\n\n");
   /*printf("                                   3.0\n\n");*/
   printf("      Copyright (c) Dr. Andrew C.R. Martin, SciTech Software \
1992-2009\n");
   printf("              Copyright (c) Dr. Craig T. Porter, UCL \
2008-2009\n\n");
}


/************************************************************************/
/*>int main(int  argc, char **argv)
   --------------------------------
   The main program. Checks parameters, sets some defaults, reads files if
   specified, initialises command parser and starts the command loop.
   
   25.09.92 Original
   29.09.92 Added gCurrentMode, Changed to call ReadStructure()
   02.10.92 Calls Cleanup() as it should have done!
   05.10.92 Added AMIGA_WINDOWS support
   06.10.92 Changed argc check for icon start
   08.10.93 Changed version number
            Chaned CURSES calls for book library
   17.07.95 Removed windowing support
   19.07.95 Handles -h flag
   20.07.95 Initialise gDoWeights
   13.06.96 Changed initialisation of gDoWeights since it now has 3 states
   18.06.96 Replaced MODE_* with ZONE_MODE_*
   01.03.01 Added -x handling
   19.03.08 Initialise gCZoneList[] By: CTP
   27.03.08 Changed call to DoCommandLoop()
   02.04.08 Initialise gRefHeader, etc... for PDB headers and footers.
            Added call to InitPDBCoordRecord().
   15.04.08 Added option to run script from command line.
   02.05.08 Removed call to InitPDBCoordRecord(). Removed gRefHeader and 
            gRefFooter
   04.06.08 Added gMatchSymAtoms and call to SetSymmetricalAtomPAirs()
   21.08.08 Added gMultiVsRef and gRotateRefit flags and calculation of 42 
            degree rotation matrix for refitting structure to avoid local 
            minimum.
   07.11.08 Added gMultiRef.
   04.02.09 gRotateRefit replaced by ROTATE_REFIT #define
*/
int main(int  argc, char **argv)
{
   BOOL CmdError = FALSE;
   BOOL xmasFormat = FALSE;
   BOOL runscript = FALSE;
   char filename[MAXSTRLEN];
   int  i;
   
   logo();

   argc--;
   argv++;

   /* See if a flag has been given                                      */
   while(argc>0 && argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 'h':
         gHetAtoms = TRUE;
         break;
      case 'x':
#ifdef USE_XMAS
         xmasFormat = TRUE;
#else
         fprintf(stderr,"Error: XMAS support not compiled in this \
version\n");
         return(1);
#endif
         break;
      case 'f':
        if(runscript)
          {
            CmdError = TRUE;
            break;
          }
        runscript = TRUE;
        filename[0] = '\0';
        strcpy(filename,argv[1]);
        argc--;
        argv++;
        break;
      default:
         CmdError = TRUE;
         break;
      }
      argc--;
      argv++;
   }
   

   if(CmdError || (argc>0 && argc!=2))
   {
      printf("\n\nSyntax: ProFit [-h] [-f <scriptfile.txt>]");
      printf(" [<reference.pdb> <mobile.pdb>]\n");
      printf("        -h Include HETATM records when reading PDB \
files\n");
      printf("        -f Run script\n\n");
      
      exit(1);
   }
   
   /* Initialise various things                                         */
   strcpy(gFitAtoms[0],"*");
   strcpy(gRMSAtoms[0],"*");
   gRefFilename[0] = '\0';
   gCurrentMode     = ZONE_MODE_RESNUM;
   gUserRMSAtoms    = FALSE;
   gUserRMSZone     = FALSE;
   gUserFitZone     = FALSE;
   gFitted          = FALSE;
   gDoWeights       = WEIGHT_NONE;
   gReadHeader      = FALSE;
   gMatchSymAtoms   = FALSE;
   gMultiVsRef      = FALSE;
   gTwistAngle      = 42.0;
   CalculateRotationMatrix(gTwistAngle, gRotMatTwist);

   gWtAverage       = TRUE;
   gMultiRef        = 0;

   for(i=0; i<MAXSTRUC; i++)
   {
      gMobPDB[i]         = NULL;
      gMobWPDB[i]        = NULL;
      gFitPDB[i]         = NULL;
      gMobCoor[i]        = NULL;
      gMobFilename[i][0] = '\0';
      gZoneList[i]       = NULL;
      gRZoneList[i]      = NULL;
      gCZoneList[i]      = NULL;
   }
   gLimit[0] = gLimit[1] = (-1);

   
   if(argc==2)
   {
      /* Filenames have been specified, so open and read files          */
      ReadStructure(STRUC_REFERENCE, argv[0], 0, xmasFormat);
      ReadStructure(STRUC_MOBILE,    argv[1], 0, xmasFormat);
      gMultiCount = 1;
   }


   InitParser();
   SetSymmetricalAtomPAirs();

   if(runscript)
      RunScript(filename);
   else 
      DoCommandLoop(stdin);
   
   Cleanup();

   return(0);
}


/************************************************************************/
/*>void SetRMSAtoms(char *command)
   -------------------------------
   Responds to a command to set the atoms for calculation of RMS

   28.09.92 Original.
   29.09.92 Added gUserRMSAtoms
   30.09.92 Added NOT option; pad gRMSAtoms[]
   17.07.95 Replaced calls to screen() with printf()
   19.07.95 Added parameter to ShowRMS()
   25.07.95 Added another parameter to ShowRMS()
   23.07.96 Changed to call PADMINTERM() rather than padterm()
            Checks for legal wildcard specifications
            Improved out-of-bounds error checking
   03.04.08 Added parameter to ShowRMS() By: CTP
   11.08.08 ShowRMS() turned-off when in quiet mode.
   13.03.09 Added temporary fix for atom type recognition.

*/
void SetRMSAtoms(char *command)
{
   int   comptr,
         atmnum,
         atmpos,
         comstart,
         strucnum;
   char  *cmnd;
   
   /* Put to upper case                                                 */
   UPPER(command);
   
   /* Assume this is not a NOT selection                                */
   gNOTRMSAtoms = FALSE;
   cmnd        = command;
   
   /* Set NOT flag if required and step over the symbol                 */
   if(command[0] == '~' || command[0] == '^')
   {
      gNOTRMSAtoms = TRUE;
      cmnd++;
   }

   /* Blank all rmsatoms                                                */
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      for(atmpos=0; atmpos<8; atmpos++)
         gRMSAtoms[atmnum][atmpos] = '\0';
   
   atmnum=0;
   atmpos=0;
   comstart=0;
   for(comptr=0;comptr<strlen(cmnd);comptr++)
   {
      if(cmnd[comptr]==',')
      {
         comstart=comptr+1;
         gRMSAtoms[atmnum][atmpos]   = '\0';  /* Terminate string       */
         PADMINTERM(gRMSAtoms[atmnum],4);
         if(++atmnum >= NUMTYPES)
         {
            if(!gQuiet)
            {
               printf("   Warning==> Too many atoms in specification. \
Only %d used\n",NUMTYPES);
            }
            break;
         }
         atmpos=0;
      }
      else
      {
         if(atmpos >= MAXATSPEC)
         {
            TERMAT(cmnd+comptr, ',');
            if(!gQuiet)
            {
               printf("   Warning==> Atom name specification too long \
(%s)\n", cmnd+comstart);
               printf("              Using all atoms.\n");
            }
            gRMSAtoms[0][0] = '*';
            gRMSAtoms[0][1] = '\0';
            break;
         }
         gRMSAtoms[atmnum][atmpos++] = cmnd[comptr];
      }
   }
   
   /* Terminate last one                                                */
   gRMSAtoms[atmnum][atmpos] = '\0';
   PADMINTERM(gRMSAtoms[atmnum],4);

   /* See if a * was specified not in first position; move it if so.
      Also check for legal atom wildcard specifications
   */
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
   {
      if(!LegalAtomSpec(gRMSAtoms[atmnum]))
      {
         if(!gQuiet)
         {
            printf("   Warning==> Illegal atom specification (%s). Using \
all atoms.\n",
                   gRMSAtoms[atmnum]);
         }
         gRMSAtoms[0][0] = '*';
         gRMSAtoms[0][1] = '\0';
         break;
      }
      
      if(gRMSAtoms[atmnum][0] == '*')
      {
         if(gNOTRMSAtoms)
         {
            if(!gQuiet)
            {
               printf("   Warning==> NOT option ignored.\n");
            }
            gNOTRMSAtoms = FALSE;
         }
         gRMSAtoms[0][0] = '*';
         gRMSAtoms[0][1] = '\0';
         break;
      }
   }
   
   gUserRMSAtoms = TRUE;
   

   /* Bugfix - Atom Recognition By: CTP */
   /* Temporary fix for atom type recognition. Removed trailing spaces
      for atom type strings. Need to check with ACRM why strings were
      padded. */
   if(TRUE)
   {
      int i = 0;
      for(i=0;i<NUMTYPES;i++)
      {
         if(gRMSAtoms[i][0] == '\0') break;
         KILLTRAILSPACES(gRMSAtoms[i]);
      }
   }
   /* End of Bugfix */


   if(!gQuiet)
   {
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         ShowRMS(FALSE,NULL,strucnum,FALSE,FALSE);
      }
   }

   return;
}


/************************************************************************/
/*>void SetRMSZone(char *command)
   ------------------------------
   Responds to a command to set the zone for RMS calculation

   29.09.92 Original. 
   17.07.95 Replaced calls to screen() with printf()
   18.07.95 Added initialisation of inserts in zones
   19.07.95 Added parameter to ShowRMS()
            * on its own equivalent to CLEAR
   25.07.95 Added another parameter to ShowRMS()
   18.06.96 Replaced MODE_* with ZONE_MODE_*
   01.02.01 Added multiple structures
   28.02.01 Error message from ParseZone() when trying to specify multiple
            zones for multi-structure now here instead of in ParseZone()
   03.04.08 Added parameter to ShowRMS() By: CTP
   11.04.08 Added warning for overlapping zones.
   22.04.08 Added handling of lowercase chain and inserts.
   11.08.08 ShowRMS() turned-off when in quiet mode.
   29.10.08 Fixed bug that gave segmentation fault with CLEAR.
*/
void SetRMSZone(char *command)
{
   int   start1,  stop1,
         start2,  stop2,
         SeqZone, strucnum;
   char  chain1,  chain2,
         startinsert1, stopinsert1,
         startinsert2, stopinsert2;
   ZONE  *z;
   int   warned = 0;

   
   /* See if this is clearing the zones                                 */
   if(!upstrncmp(command,"CLEAR",5) || !strcmp(command,"*"))
   {
      gUserRMSZone = FALSE;
      
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         if(gRZoneList[strucnum]!=NULL)
         {
            FREELIST(gRZoneList[strucnum],ZONE);
            gRZoneList[strucnum] = NULL;
         }
         
         if(!gQuiet)
         {
            ShowRMS(FALSE,NULL,strucnum,FALSE,FALSE);
         }
      }
      return;
   }
   
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      SeqZone = ParseZone(command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1, 
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2,
                          strucnum);

      if((SeqZone == (-2)) && (!warned))
      {
         printf("   Error==> You cannot specify zones for each \
structure when performing\n");
         printf("            multiple structure fitting.\n");
         warned = 1;
      }
      
      if(SeqZone > -1)
      {
         /* If the user has not already specified an RMS zone, blank
            the current list.
         */
         if(!gUserRMSZone)
         {
            if(gRZoneList[strucnum]!=NULL)
            {
               FREELIST(gRZoneList[strucnum],ZONE);
               gRZoneList[strucnum] = NULL;
            }
         }
         
         /* Allocate an entry in zone list                              */
         if(gRZoneList[strucnum])
         {
            /* Move to end of zone list                                 */
            z=gRZoneList[strucnum];
            LAST(z);
            ALLOCNEXT(z,ZONE);
         }
         else
         {
            INIT(gRZoneList[strucnum],ZONE);
            z = gRZoneList[strucnum];
         }
         
         if(z==NULL)
         {
            printf("   Error==> No memory for zone!\n");
         }
         else
         {
            /* Add this zone to the zone list                           */
            z->chain1       = chain1;
            z->start1       = start1;
            z->startinsert1 = startinsert1;
            z->stop1        = stop1;
            z->stopinsert1  = stopinsert1;
            z->chain2       = chain2;
            z->start2       = start2;
            z->startinsert2 = startinsert2;
            z->stop2        = stop2;
            z->stopinsert2  = stopinsert2;
            z->mode         = SeqZone?ZONE_MODE_SEQUENTIAL:gCurrentMode;

            /* CTP: Check for overlap                                   */
            if(CheckOverlap(z,gRZoneList[strucnum],strucnum) > 1)
              printf("   Warning: New zone overlaps existing zone.\n");
            else if(CheckOverlap(z,gRZoneList[strucnum],strucnum) == -1)
              printf("   Error: Failed to find new zone.\n");
         }
         
         gUserRMSZone = TRUE;
         if(!gQuiet)
         {
            ShowRMS(FALSE,NULL,strucnum,FALSE,FALSE);
         }
      }
   }

   return;
}


/************************************************************************/
/*>void Cleanup(void)
   ------------------
   Clean up screen and any malloc'd memory

   28.09.92 Original
   29.09.92 Added sequences
   01.10.92 Freed strparam and *Coor memory.
   05.10.92 Added libraries for AMIGA_WINDOWS
   08.10.93 Changed CURSES calls for book library
   17.07.95 Removed CURSES clean up and Amiga Windows cleanup
   20.07.95 Frees gWeights
   12.01.01 gMobPDB[] now an array
   01.02.01 Various other arrays now freed for multi-structure support
   28.03.08 Added gCZoneList[] By: CTP
   02.04.08 Added Header and Footer storage.
   02.04.08 Removed Header and Footer storage.
*/
void Cleanup(void)
{
   int   i;
   
   for(i=0; i<MAXSTRUC; i++)
   {
      if(gMobPDB[i])  FREELIST(gMobPDB[i], PDB);
      gMobPDB[i]     = NULL;
   }

   if(gRefPDB)          FREELIST(gRefPDB, PDB);
   gRefPDB     = NULL;
   
   for(i=0; i<MAXSTRUC; i++)
   {
      if(gFitPDB[i])          FREELIST(gFitPDB[i], PDB);
      gFitPDB[i]     = NULL;

      if(gMobCoor[i])         free(gMobCoor[i]);
      gMobCoor[i]    = NULL;

      if(gMobSeq[i])          free(gMobSeq[i]);
      gMobSeq[i]     = NULL;

      if(gZoneList[i])        FREELIST(gZoneList[i], ZONE);
      gZoneList[i]   = NULL;

      if(gRZoneList[i])       FREELIST(gRZoneList[i], ZONE);
      gRZoneList[i]  = NULL;

      if(gCZoneList[i])       FREELIST(gCZoneList[i], ZONE);
      gCZoneList[i]  = NULL;
   }

   if(gRefSeq)          free(gRefSeq);
   gRefSeq     = NULL;
   
   if(gRefCoor)         free(gRefCoor);
   gRefCoor    = NULL;
   
   if(gWeights)         free(gWeights);
   gWeights    = NULL;
   
   for(i=0; i<MAXSTRPARAM; i++)
      free(gStrParam[i]);
   
   Help("Dummy","CLOSE");
}


/************************************************************************/
/*>void Die(char *message)
   -----------------------
   Program death. Writes message and cleans up before exit.

   25.09.92 Original
   17.07.95 Replaced calls to screen() with puts()
*/
void Die(char *message)
{
   puts(message);
   
   Cleanup();
   
   exit(0);
}


/************************************************************************/
/*>BOOL ReadStructure(int structure, char *filename, int strucnum,
                      int xmasFormat)
   ---------------------------------------------------------------
   Reads one of the 2 PDB structures.

   28.09.92 Original
   30.09.92 Added coordinate array allocation & check for inserts
   08.10.93 Modified for new ReadPDB()
   24.11.94 Changed to call ReadPDBAtoms() rather than ReadPDB()
   17.07.95 Replaced calls to screen() with printf()
   18.07.95 Uses fopen() rather than OpenRead()
            Uses fclose() rather than CloseFile()
   19.07.95 Selects whether to read HETATM records depending on
            gHetAtoms
   24.01.96 Improved error messages when no atoms read since ReadPDB()
            can distinguish between no memory and no atoms.
   03.07.97 Warning message if finds multiple occupancies
            Changed to CD1 rather than CD for ILE
   12.01.01 gMobPDB[] now an array
   01.02.01 Added strucnum param; returns BOOL
   20.02.01 Frees all previously loaded structures if we had a multiple
            set loaded previously.
   01.03.01 Added xmasFormat parameter
   06.07.08 Correction of CD -> CD1 had an extra space in the raw atom 
            name (i.e. it was " CD1 " rather than " CD1") By: CTP
   02.04.08 CTP added call to ReadPDBHeader()
   29.04.08 CTP added calls to ReadWholePDB() for reading pdb structures.
            Removed call to ReadPDBHeader().
   01.05.08 CTP added support for reading low occupancy atoms.
   02.05.08 CTP header and trailer now stored within wholepdb datatype.
   12.06.08 Changed call from doReadWholePDB() to ReadWholePDB() and 
            ReadWholePDBAtom(). 
            Added calls to ReadPDBOccRank() and  ReadPDBAtomsOccRank() 
            for reading partial occupancies.
   13.03.09 Set function to clear fitted coordinates when loading a new
            mobile structure. Required for the NOFIT command which 
            (unlike the FIT command) does not update any existing fitted
            coordinates.
*/
BOOL ReadStructure(int  structure,
                   char *filename,
                   int  strucnum,
                   int  xmasFormat)
{
   file  *fp;
   int   natoms;
   PDB   *p;
   int   i;
   static int sLastStrucNum = (-1);

   gPDBPartialOcc = FALSE;
   
   /* Open the file                                                     */
   if((fp = fopen(filename,"r"))==NULL)
   {
      printf("   Unable to open file %s\n",filename);
      return(FALSE);
   }

   /* Handle appropriate pdb linked list and give message               */
   if(structure==STRUC_REFERENCE)
   {
      printf("   Reading reference structure...\n");
      strcpy(gRefFilename,filename);

      /* CTP: If there is something in the current PDB list, free it    */
      if(gRefWPDB)
      {
         /* Free current WHOLEPDB.                                      */
         FreeWholePDB(gRefWPDB);
         gRefWPDB   = NULL;
         gRefPDB    = NULL;
         natoms     = 0;
      }
      /* CTP: Free current PDB (Needed if using XMAS format)            */
      if(gRefPDB)
      {
         FREELIST(gRefPDB,PDB);
         gRefPDB = NULL;
         gRefWPDB   = NULL;
         natoms     = 0;
      }
   
      /* Read the structure                                             */
      if(xmasFormat)
      {
         if(gHetAtoms)
            gRefPDB = ReadXMAS(fp, &natoms);
         else
            gRefPDB = ReadXMASAtoms(fp, &natoms);
      }
      else
      {
         if(gHetAtoms)
            gRefWPDB = ReadWholePDB(fp);
         else
            gRefWPDB = ReadWholePDBAtoms(fp);
         
         if(gOccRank > 1)
         {
            FREELIST(gRefWPDB->pdb, PDB);
            gRefWPDB->pdb = NULL;
            rewind(fp);
            
            if(gHetAtoms)
               gRefWPDB->pdb = ReadPDBOccRank(fp,&natoms,gOccRank);
            else 
               gRefWPDB->pdb = ReadPDBAtomsOccRank(fp,&natoms,gOccRank);
         }
         
         gRefWPDB->pdb = RemoveAlternates(gRefWPDB->pdb);
         natoms     = gRefWPDB->natoms;
         gRefPDB    = gRefWPDB->pdb;
      }
      
      if(gRefPDB==NULL)
      {
         if(!natoms)
         {
            printf("   Error==> No atoms read from reference PDB \
file!\n");
         }
         else
         {
            printf("   Error==> No memory to read reference PDB file!\n");
         }
            
         gRefFilename[0] = '\0';
         fclose(fp);
         return(FALSE);
      }

      /* 03.07.97 Warning about multiple-occupancies                    */
      if(gPDBPartialOcc && !gQuiet)
      {
         printf("   Warning==> Reference set contains multiple occupancy \
atoms.\n");
         if(gOccRank == 1)
           printf("              Only the first atom will be \
considered.\n");
      }      

      /* Allocate coordinate array                                      */
      if(gRefCoor) free(gRefCoor);
      if((gRefCoor = (COOR *)malloc(natoms * sizeof(COOR))) == NULL)
         printf("   Error==> Unable to allocate reference coordinate \
memory!\n");

      /* Convert to sequence                                            */
      if(gRefSeq != NULL) free(gRefSeq);
      if((gRefSeq = PDB2Seq(gRefPDB))==NULL)
         printf("   Error==> Unable to read sequence for reference \
structure!\n");
      
      /* Check for inserts                                              */
      for(p=gRefPDB; p!=NULL; NEXT(p))
      {
         if(p->insert[0] != ' ')
         {
            if(!gQuiet)
            {
               printf("   Warning==> Reference protein contains \
insertions.\n");
            }
            break;
         }
      }
      
      /* Fix ILE CD to CD1 (03.07.97 was other way round)               */
      for(p=gRefPDB; p!=NULL; NEXT(p))
      {
         /* 28.02.01 added ->atnam_raw                                  */
         if(!strncmp(p->resnam,"ILE ",4) && !strncmp(p->atnam,"CD  ",4))
         {
            strcpy(p->atnam,"CD1 ");
            strcpy(p->atnam_raw," CD1"); /* 06.07.06 Removed extra space*/
         }
      }
   }
   else if(structure==STRUC_MOBILE)
   {
      printf("   Reading mobile structure...\n");
      strcpy(gMobFilename[strucnum],filename);
      
      /* If there is something in the current list, free it             */
      for(i=strucnum; i<=sLastStrucNum; i++)
      {
         if(gMobWPDB[i])
         {
            FreeWholePDB(gMobWPDB[i]);
            gMobWPDB[i]   = NULL;
            gMobPDB[i]    = NULL;
            natoms        = 0;
         }       
         if(gMobPDB[i])
         {
            FREELIST(gMobPDB[i],PDB);
            gMobPDB[i]    = NULL;
            gMobWPDB[i]   = NULL;
            natoms        = 0;
         }       
      }
      sLastStrucNum = strucnum;
   

      /* Read the structure                                             */
      if(xmasFormat)
      {
         if(gHetAtoms)
            gMobPDB[strucnum] = ReadXMAS(fp, &natoms);
         else
            gMobPDB[strucnum] = ReadXMASAtoms(fp, &natoms);
      }
      else
      {
         if(gHetAtoms)
            gMobWPDB[strucnum]   = ReadWholePDB(fp);
         else 
            gMobWPDB[strucnum]   = ReadWholePDBAtoms(fp);
         
         if(gOccRank > 1)
         {
            FREELIST(gMobWPDB[strucnum]->pdb, PDB);
            gMobWPDB[strucnum]->pdb = NULL;
            rewind(fp);
            
            if(gHetAtoms)
               gMobWPDB[strucnum]->pdb = 
                  ReadPDBOccRank(fp,&natoms,gOccRank);
            else 
               gMobWPDB[strucnum]->pdb = 
                  ReadPDBAtomsOccRank(fp,&natoms,gOccRank);
         }
         
         gMobWPDB[strucnum]->pdb = 
            RemoveAlternates(gMobWPDB[strucnum]->pdb);
         natoms               = gMobWPDB[strucnum]->natoms;
         gMobPDB[strucnum]    = gMobWPDB[strucnum]->pdb;
      }
      
      
      if(gMobPDB[strucnum]==NULL)
      {
         if(!natoms)
         {
            printf("   Error==> No atoms read from mobile PDB file!\n");
         }
         else
         {
            printf("   Error==> No memory to read mobile PDB file!\n");
         }
         
         gMobFilename[strucnum][0] = '\0';
         fclose(fp);
         return(FALSE);
      }

      /* 03.07.97 Warning about multiple-occupancies                    */
      if(gPDBPartialOcc && !gQuiet)
      {
         printf("   Warning==> Mobile set contains multiple occupancy \
atoms.\n");
         if(gOccRank == 1)
            printf("              Only the first atom will be \
considered.\n");
      }      

      /* Allocate coordinate array                                      */
      if(gMobCoor[strucnum]) free(gMobCoor[strucnum]);
      if((gMobCoor[strucnum] = (COOR *)malloc(natoms * sizeof(COOR))) 
         == NULL)
         printf("   Error==> Unable to allocate mobile coordinate \
memory!\n");
      
      /* Convert to sequence                                            */
      if(gMobSeq[strucnum] != NULL) free(gMobSeq[strucnum]);
      if((gMobSeq[strucnum] = PDB2Seq(gMobPDB[strucnum]))==NULL)
         printf("   Error==> Unable to read sequence for reference \
structure!\n");
      
      /* Check for inserts                                              */
      for(p=gMobPDB[strucnum]; p!=NULL; NEXT(p))
      {
         if(p->insert[0] != ' ')
         {
            if(!gQuiet)
            {
               printf("   Warning==> Mobile protein contains \
insertions.\n");
            }
            break;
         }
      }
      
      /* Fix ILE CD to CD1 (03.07.97 was other way round)               */
      for(p=gMobPDB[strucnum]; p!=NULL; NEXT(p))
      {
         /* 28.02.01 added ->atnam_raw                                  */
         if(!strncmp(p->resnam,"ILE ",4) && !strncmp(p->atnam,"CD  ",4))
         {
            strcpy(p->atnam,"CD1 ");
            strcpy(p->atnam_raw," CD1"); /* 06.07.06 Removed extra space*/
         }
      }

      /* Free Fitted Structure */
      if(gFitPDB[strucnum] != NULL) FREELIST(gFitPDB[strucnum], PDB);
      gFitPDB[strucnum] = NULL;

   }
   else
   {
      printf("   ReadStructure(): Internal Error!\n");
      return(FALSE);
   }

   gUserRMSAtoms = FALSE;
   gUserRMSZone  = FALSE;
   gFitted       = FALSE;

   /* and close file                                                    */
   fclose(fp);

   return(TRUE);
}


/************************************************************************/
/*>int InitParser(void)
   --------------------
   Initialise the command parser.
   
   25.09.92 Original
   09.10.92 Changed ALIGN to 0 parameters
   19.07.95 Added RESIDUE, HETATOMS & NOHETATOMS
   20.07.95 Added WEIGHT/NOWEIGHT
   22.07.95 Added GAPPEN
   20.11.95 Added READALIGNMENT
   31.05.96 Added BVALUE
   13.06.96 Added BWEIGHT
   11.11.96 BVALUE command now takes 1 or 2 params
   18.11.96 Added (NO)IGNOREMISSING
   20.12.96 Added NFITTED
   12.01.01 Added ITERATE
   01.02.01 Added MULTI, QUIET and MWRITE
   20.02.01 Added LIMIT
   01.03.01 Added optional 'xmas' parameter to REF, MOB, MULTI if
            XMAS support is compiled in
   20.03.01 Added CENTER/CENTRE commands
   28.03.01 Added optional 'reference' parameter to WRITE
   14.03.08 Added DELZONE command. By: CTP
   17.03.08 Added DELRZONE command.
   19.03.08 Added SETCENTRE/SETCENTER command.
   27.03.08 Added SCRIPT command.
   02.04.08 Added HEADER command.
   03.04.08 Added PAIRDIST command.
   07.04.08 Added DISTCUTOFF command.
   15.04.08 Added RETURN command
   04.06.08 Added OCCRANK, BZONE and SYMMATOMS commands.
            Removed RETURN command. 
   16.07.08 GAPPEN command now takes 1 or 2 params
   18.07.08 Altered ALIGN command to take string parameters.
            Added PRINTALIGN command.
   07.08.08 Added SETMULTI command.
   16.10.08 Added ALLVSALL, SETREF, ORDERFIT and TRIMZONES commands.
   23.10.08 Added WTAVERAGE
   27.10.08 Added NOFIT.
   07.11.08 Changed STATUS to take 0 or 1 string parameters.

*/
int InitParser(void)
{
   int   i;
   
   /* Initialise returned string array                                  */
   for(i=0; i<MAXSTRPARAM; i++)
      gStrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char));

   /* Construct the gKeyWords                                           */
#ifdef USE_XMAS
   MAKEMKEY(gKeyWords[0],  "REFERENCE",       STRING,1,2);
   MAKEMKEY(gKeyWords[1],  "MOBILE",          STRING,1,2);
#else
   MAKEMKEY(gKeyWords[0],  "REFERENCE",       STRING,1,1);
   MAKEMKEY(gKeyWords[1],  "MOBILE",          STRING,1,1);
#endif
   MAKEMKEY(gKeyWords[2],  "FIT",             NUMBER,0,0);
   MAKEMKEY(gKeyWords[3],  "ATOMS",           STRING,1,1);
   MAKEMKEY(gKeyWords[4],  "ZONE",            STRING,1,1);
   MAKEMKEY(gKeyWords[5],  "GRAPHIC",         NUMBER,0,0);
   MAKEMKEY(gKeyWords[6],  "ALIGN",           STRING,0,2);
   MAKEMKEY(gKeyWords[7],  "RATOMS",          STRING,1,1);
   MAKEMKEY(gKeyWords[8],  "RZONE",           STRING,1,1);
   MAKEMKEY(gKeyWords[9],  "WRITE",           STRING,1,2);
   MAKEMKEY(gKeyWords[10], "MATRIX",          NUMBER,0,0);
   MAKEMKEY(gKeyWords[11], "STATUS",          STRING,0,1);
   MAKEMKEY(gKeyWords[12], "QUIT",            NUMBER,0,0);
   MAKEMKEY(gKeyWords[13], "NUMBER",          STRING,1,1);
   MAKEMKEY(gKeyWords[14], "RMS",             NUMBER,0,0);
   MAKEMKEY(gKeyWords[15], "RESIDUE",         STRING,0,1);
   MAKEMKEY(gKeyWords[16], "HETATOMS",        NUMBER,0,0);
   MAKEMKEY(gKeyWords[17], "NOHETATOMS",      NUMBER,0,0);
   MAKEMKEY(gKeyWords[18], "WEIGHT",          NUMBER,0,0);
   MAKEMKEY(gKeyWords[19], "NOWEIGHT",        NUMBER,0,0);
   MAKEMKEY(gKeyWords[20], "GAPPEN",          NUMBER,1,2);
   MAKEMKEY(gKeyWords[21], "READALIGNMENT",   STRING,1,1);
   MAKEMKEY(gKeyWords[22], "BVALUE",          STRING,1,2);
   MAKEMKEY(gKeyWords[23], "BWEIGHT",         NUMBER,0,0);
   MAKEMKEY(gKeyWords[24], "IGNOREMISSING",   NUMBER,0,0);
   MAKEMKEY(gKeyWords[25], "NOIGNOREMISSING", NUMBER,0,0);
   MAKEMKEY(gKeyWords[26], "NFITTED",         NUMBER,0,0);
   MAKEMKEY(gKeyWords[27], "ITERATE",         STRING,0,1);
#ifdef USE_XMAS
   MAKEMKEY(gKeyWords[28], "MULTI",           STRING,1,2);
#else
   MAKEMKEY(gKeyWords[28], "MULTI",           STRING,1,1);
#endif
   MAKEMKEY(gKeyWords[29], "QUIET",           STRING,0,1);
   MAKEMKEY(gKeyWords[30], "MWRITE",          STRING,0,1);
   MAKEMKEY(gKeyWords[31], "LIMIT",           STRING,1,2);
   MAKEMKEY(gKeyWords[32], "CENTRE",          STRING,0,1);
   MAKEMKEY(gKeyWords[33], "CENTER",          STRING,0,1);
   MAKEMKEY(gKeyWords[34], "DELZONE",         STRING,1,1);
   MAKEMKEY(gKeyWords[35], "DELRZONE",        STRING,1,1);
   MAKEMKEY(gKeyWords[36], "SETCENTRE",       STRING,1,1);
   MAKEMKEY(gKeyWords[37], "SETCENTER",       STRING,1,1);
   MAKEMKEY(gKeyWords[38], "SCRIPT",          STRING,1,1);
   MAKEMKEY(gKeyWords[39], "HEADER",          STRING,0,1);
   MAKEMKEY(gKeyWords[40], "PAIRDIST",        STRING,0,1);
   MAKEMKEY(gKeyWords[41], "DISTCUTOFF",      STRING,0,1);
   MAKEMKEY(gKeyWords[42], "OCCRANK",         NUMBER,1,1);
   MAKEMKEY(gKeyWords[43], "BZONE",           NUMBER,0,0);
   MAKEMKEY(gKeyWords[44], "SYMMATOMS",       STRING,0,2);
   MAKEMKEY(gKeyWords[45], "PRINTALIGN",      STRING,0,2);
   MAKEMKEY(gKeyWords[46], "MULTREF",         STRING,0,1);
   MAKEMKEY(gKeyWords[47], "ALLVSALL",        STRING,0,1);
   MAKEMKEY(gKeyWords[48], "SETREF",          NUMBER,0,1);
   MAKEMKEY(gKeyWords[49], "ORDERFIT",        STRING,0,0);
   MAKEMKEY(gKeyWords[50], "TRIMZONES",       STRING,0,0);
   MAKEMKEY(gKeyWords[51], "WTAVERAGE",       STRING,0,1);
   MAKEMKEY(gKeyWords[52], "NOFIT",           STRING,0,0);

   return(0);
}


/************************************************************************/
/*>int DoCommandLoop(FILE *script)
   -------------------------------
   Sit and wait for commands to be processed by the parser. Also checks 
   for help commands.
   
   25.09.92 Original
   28.09.92 Added ReadStructure(), OS commands
   02.10.92 Changed gets() to GetKybdString()
   17.07.95 Replaced calls to screen() with printf()
            Back to fgets() rather than GetKybdString()
   19.07.95 Added ByRes parameter to ShowRMS() and added RESIDUE, HETATOMS
            and NOHETATOMS
   20.07.95 Added WEIGHT/NOWEIGHT
   25.07.95 Changed to use mparse() and RESIDUE command may now take
            a filename for output. Added filename parameter to ShowRMS()
   20.11.95 Added READALIGNMENT handling
   21.11.95 Corrected `defaut:' to `default:'
   31.05.96 Added BVALUE handling
   13.06.96 Added BWEIGHT handling
   11.11.96 Handles second parameter to BVALUE
   18.11.96 Added IGNOREMISSING handling
   20.12.96 Added NFITTED handling
   15.01.01 Added ITER handling
   01.02.01 Added MULTI, QUIET, MWRITE handling
   20.02.01 Added LIMIT handling
   01.03.01 Added xmas parameters to REF, MOB, MULTI
   20.03.01 Added numeric parameter for ITERATE
            Added CENTRE/CENTER
   28.03.01 Added REFERENCE parameter to WRITE command
   14.03.08 Added DELZONE command. By: CTP
   17.03.08 Added DELRZONE command.
   19.03.08 Added SETCENTRE/SETCENTER command.
   27.03.08 Added SCRIPT command. Added parameter - DoCommandLoop() now
            takes input from FILE (either stdin or from script file).
            Added comment marker # 
            Prompt only displayed when input from stdin.
   02.04.08 Added HEADER command.
   03.04.08 Added PAIRDIST command. Added parameter to ShowRMS()
   07.04.08 Added DISTCUTOFF command.
   16.04.08 Quitting anywhere exits from ProFit.
   04.06.08 Added SYMMATOMS command.
   16.07.08 GAPPEN now takes one or two parameters.
   18.07.08 Altered ALIGN command to call AlignWrapper().
            Added PRINTALIGN command.
   24.07.08 Set RMS, RESIDUE and PAIRDIST commands to cycle through all 
            mobile structures.
   07.08.08 Added MULTIREF command. 
   16.10.08 Added ALLVSALL, SETREF, ORDERFIT and TRIMZONES commands.
   23.10.08 Added WTAVERAGE command
   27.10.08 Added NOFIT command.
   07.11.08 Changed STATUS to take optional output filename.
   14.11.08 Added GNU Readline Library support.
   14.01.09 Added Output to file for PRINTALIGN.

*/
int DoCommandLoop(FILE *script)
{
   int   nletters,
         NParam,
         i;
   char  comline[MAXBUFF],
         *ptr;
   
   for(;;)
   {
#ifdef READLINE_SUPPORT
      if(script == stdin)
      {
         /* Get user input using readline                               */
         char *line_read;
         line_read = readline("ProFit>");

         if(line_read && *line_read)
            add_history (line_read);
     
         strcpy(comline, line_read);
         TERMINATE(comline);
         
         if(line_read)
         {
            free(line_read);
            line_read =(char *)NULL;
         }
      }
      else 
      {
         /* Read script using fgets                                     */
         if(!fgets(comline, MAXBUFF, script))
            return(0);
         TERMINATE(comline);
      }
#else
      /* Get all input using fgets                                      */
      if(script == stdin)
         stdprompt("ProFit");
      
      if(!fgets(comline, MAXBUFF, script))
         return(0);
      TERMINATE(comline);
#endif
     
      /* Any line which starts with a $ is passed to the OS             */
      if(comline[0] == '$')
      {
         ptr = comline+1;
         system(ptr);
         continue;
      }

      /* Any line which starts with a # is echoed to the screen         */
      ptr = KillLeadSpaces(comline);
      if(ptr[0] == '#')
      {
         if(!gQuiet)
            printf("   %s\n",comline);
         continue;
      }

      /* We need to check HELP outside the main parser as it has an 
         optional parameter (N.B. mparse() could now handle this)
      */
      if(match(comline,"HELP",&nletters))
      {
         if(nletters == 4)
         {
            DoHelp(comline,HELPFILE);
            continue;
         }
      }
      
      /* Main parser                                                    */
      switch(mparse(comline,NCOMM,gKeyWords,gNumParam,gStrParam,&NParam))
      {
      case PARSE_ERRC:
         printf("   Unrecognised keyword: %s\n",comline);
         break;
      case PARSE_ERRP:
         printf("   Invalid parameters: %s\n",comline);
         break;
      case PARSE_COMMENT:
         break;
      case 0:  /* REFERENCE                                             */
         gMultiCount = 1;
         if(NParam==1)
            ReadStructure(STRUC_REFERENCE,gStrParam[0],0,FALSE);
         else
            ReadStructure(STRUC_REFERENCE,gStrParam[1],0,TRUE);
         break;
      case 1:  /* MOBILE                                                */
         gMultiCount = 1;
         if(NParam==1)
            ReadStructure(STRUC_MOBILE,gStrParam[0],0,FALSE);
         else
            ReadStructure(STRUC_MOBILE,gStrParam[1],0,TRUE);
         break;
      case 2:  /* FIT                                                   */
         FitStructures();
         break;
      case 3:  /* ATOMS                                                 */
         if(gIterate)
         {
            printf("   Warning==> You cannot change the atoms when \
ITERATE is set.\n");
            printf("              Command ignored\n");
         }
         else
         {
            SetFitAtoms(gStrParam[0]);
         }
         break;
      case 4:  /* ZONE                                                  */
         SetFitZone(gStrParam[0], -1);
         break;
      case 5:  /* GRAPHIC                                               */
         GraphicAlign();
         break;
      case 6:  /* ALIGN                                                 */
        /* CTP: Change to new alignment code and new command            */
        /*
         for(i=0; i< gMultiCount; i++)
            NWAlign(i);
         break;
        */
         if(NParam==0)
         {
            AlignmentWrapper(-1,NULL,FALSE);
         }
         else if(NParam==1) 
         {
            AlignmentWrapper(-1,gStrParam[0],FALSE);
         }
         else if(NParam==2 && !upstrncmp(gStrParam[1],"APPEND",6))
         {
            AlignmentWrapper(-1,gStrParam[0],TRUE);
         }
         else 
         {
            printf("   Error: Command format is ");
            printf("ALIGN [[WHOLE|*]|zone [APPEND]]\n");
         }
         break;
      case 7:  /* RATOMS                                                */
         SetRMSAtoms(gStrParam[0]);
         break;
      case 8:  /* RZONE                                                 */
         SetRMSZone(gStrParam[0]);
         break;   
      case 9:  /* WRITE                                                 */
         if(NParam == 2)
         {
            if(!upstrncmp(gStrParam[0],"REF", 3))
            {
               WriteCoordinates(gStrParam[1], -1);
            }
            else
            {
               printf("   Error==> Invalid qualifier for WRITE: %s\n",
                      gStrParam[0]);
            }
         }
         else
         {
            WriteCoordinates(gStrParam[0], 0);
         }
         break;
      case 10: /* MATRIX                                                */
         ShowMatrix();
         break;
      case 11: /* STATUS                                                */
         ShowStatus(NParam?gStrParam[0]:NULL);
         break;
      case 12: /* QUIT                                                  */
         Cleanup();
         exit(0);
         break;
      case 13: /* NUMBER                                                */
         SetZoneStatus(gStrParam[0]);
         break;
      case 14: /* RMS                                                   */
         if(gFitted)
         {
            if(gMultiCount == 1)
            {
               ShowRMS(FALSE,NULL,0,FALSE,FALSE);
            }
            else 
            {
               int i;
               for(i=0;i<gMultiCount;i++)
               {
                  ShowRMS(FALSE,NULL,i,FALSE,FALSE);
               }
            }
         }
         else 
         {
            printf("   Warning: Structures have not been fitted.\n");
         }
         break;
      case 15: /* RESIDUE                                               */
         if(gFitted)
         {
            if(gMultiCount == 1)
            {
               ShowRMS(TRUE,(NParam?gStrParam[0]:NULL),0,FALSE,FALSE);
            }
            else 
            {
               int i;
               for(i=0;i<gMultiCount;i++)
                  ShowRMS(TRUE,(NParam?gStrParam[0]:NULL),i,FALSE,FALSE);
            }
         }
         else 
         {
            printf("   Warning: Structures have not been fitted.\n");
         }
         break;
      case 16: /* HETATOMS                                              */
         gHetAtoms = TRUE;
         printf("   Hetatoms will be read with future MOBILE or \
REFERENCE commands\n");
         break;
      case 17: /* NOHETATOMS                                            */
         gHetAtoms = FALSE;
         printf("   Hetatoms will be ignored with future MOBILE or \
REFERENCE commands\n");
         break;
      case 18: /* WEIGHT                                                */
         gDoWeights = WEIGHT_BVAL;
         break;
      case 19: /* NOWEIGHT                                              */
         gDoWeights = WEIGHT_NONE;
         break;
      case 20: /* GAPPEN                                                */
         gGapPen = (int)gNumParam[0];
         if(NParam == 2)
         {
            gGapPenExt = (int)gNumParam[1];
         }
         break;
      case 21: /* READALIGNMENT                                         */
         ReadAlignment(gStrParam[0]);
         break;
      case 22: /* BVALUE                                                */
         if((!isdigit(gStrParam[0][0])) && 
            (gStrParam[0][0] != '-')    &&
            (gStrParam[0][0] != '+')    &&
            (gStrParam[0][0] != '.'))
         {
            gUseBVal = 0;
            printf("   Atoms will be included regardless of B-value\n");
         }
         else
         {
            if(NParam == 2)
            {
               if(!upstrncmp(gStrParam[1],"REF",3))
               {
                  gUseBVal = 2;
               }
               else if(!upstrncmp(gStrParam[1],"MOB",3))
               {
                  gUseBVal = 3;
               }
               else
               {
                  printf("   %s is not a valid parameter to BVALUE. \
Command ignored.\n",gStrParam[1]);
                  break;
               }
            }
            else
            {
               gUseBVal = 1;
            }
            sscanf(gStrParam[0], "%lf", &gBValue);
         }
         break;
      case 23: /* BWEIGHT                                               */
         gDoWeights = WEIGHT_INVBVAL;
         break;
      case 24: /* IGNOREMISSING                                         */
         gIgnoreMissing = TRUE;
         break;
      case 25: /* NOIGNOREMISSING                                       */
         gIgnoreMissing = FALSE;
         break;
      case 26: /* NFITTED                                               */
         ShowNFitted();
         break;
      case 27: /* ITERATE                                               */
         if((NParam == 1) && !upstrncmp(gStrParam[0], "OFF",3))
         {
            gIterate = FALSE;
         }
         else
         {
            int ok = TRUE;

            /* If a parameter specified use this as the cutoff for adding
               or removing equivalenced pairs
            */
            if(NParam == 1)
            {
               if(!sscanf(gStrParam[0], "%lf", &gMaxEquivDistSq))
               {
                  printf("   Error==> Couldn't read cutoff from \
parameter\n");
                  break;
               }
               gMaxEquivDistSq *= gMaxEquivDistSq;
            }
            
            /* CTP: Removed this check as we now handle multiple chains */
            /* Check for numbers of chains                              */
/***
            if(countchar(gRefSeq,'*') > 0)
            {
               printf("   Error==> Structures must have only one chain \
for iterative zones\n");
               ok = FALSE;
            }
            for(i=0; i<gMultiCount; i++)
            {
               if(countchar(gMobSeq[i],'*') > 0)
               {
                  printf("   Error==> Structures must have only one \
chain for iterative zones\n");
                  ok = FALSE;
                  break;
               }
            }
***/

            if(ok)
            {
               gIterate = TRUE;
               SetFitAtoms("CA");
               if(!gQuiet)
               {
                  printf("   Info==> Setting atom selection to CA \
only\n");
               }
            }
         }
         break;
      case 28: /* MULTI                                                 */
         if(NParam==1)
            ReadMulti(gStrParam[0], FALSE);
         else
            ReadMulti(gStrParam[0], TRUE);
         break;
      case 29: /* QUIET                                                 */
         if((NParam == 1) && !upstrncmp(gStrParam[0], "OFF",3))
         {
            gQuiet = FALSE;
         }
         else
         {
            gQuiet = TRUE;
         }
         break;
      case 30: /* MWRITE                                                */
         if(NParam == 1)
         {
            WriteMulti(gStrParam[0]);
         }
         else
         {
            WriteMulti("fit");
         }
         break;
      case 31: /* LIMIT                                                 */
         if((NParam == 1) && (!upstrncmp(gStrParam[0], "OFF",3)))
         {
            gLimit[0] = gLimit[1] = (-1);
            break;
         }
         else if(NParam == 2)
         {
            if(sscanf(gStrParam[0], "%d", &(gLimit[0])) &&
               sscanf(gStrParam[1], "%d", &(gLimit[1])))
            {
               break;
            }
         }

         /* Error message                                               */
         printf("   Invalid parameters: %s\n",comline);
         break;
      case 32: /* CENTER & CENTRE                                       */
      case 33:
         if((NParam == 1) && (!upstrncmp(gStrParam[0], "OFF",3)))
         {
            gCentre = FALSE;
         }
         else
         {
            gCentre = TRUE;
         }
         break;
      case 34: /* DELZONE                                               */
         DelFitZone(gStrParam[0], -1);
         break;
      case 35: /* DELRZONE                                              */
         DelRMSZone(gStrParam[0], -1);
         break;
      case 36: /* SETCENTRE                                             */
      case 37: /* SETCENTER                                             */
         SetCentreResidue(gStrParam[0]);
         break;
      case 38: /* SCRIPT                                                */
         RunScript(gStrParam[0]);
         break;
      case 39: /* HEADER                                                */
         if((NParam == 1) && !upstrncmp(gStrParam[0], "OFF",3))
         {
            gReadHeader = FALSE;
         }
         else 
         {
            gReadHeader = TRUE;
         }
         break;
      case 40: /* PAIRDIST                                              */
         if(gFitted)
         {
            if(gMultiCount == 1)
            {
               ShowRMS(TRUE,(NParam?gStrParam[0]:NULL),0,FALSE,TRUE);
            }
            else 
            {
               int i = 0;
               for(i = 0; i<gMultiCount; i++)
               {
                  printf("\n   Mobile Structure: %d\n",i+1);
                  ShowRMS(TRUE,(NParam?gStrParam[0]:NULL),i,FALSE,TRUE);
               }
            }
         }
         else 
         {
            printf("   Warning: Structures have not been fitted.\n");
         }
         
         break;
      case 41: /* DISTCUTOFF                                            */
         if((NParam == 1) && !upstrncmp(gStrParam[0], "OFF",3))
         {
            gUseDistCutoff = FALSE;
         }
         else if((NParam == 1) && !upstrncmp(gStrParam[0], "ON",2))
         {
            gUseDistCutoff = TRUE;
         }
         else
         {
            sscanf(gStrParam[0], "%lf", &gDistCutoff);
            gUseDistCutoff = TRUE;          
         }
         if(!gQuiet)
         {
            if(gUseDistCutoff)
            {
               printf("   Atom pairs will be discarded if their \
interatomic distance is > %.2f\n", gDistCutoff); 
               if(!gDistCutoff)
               {
                  printf("   Warning: No atom pairs will be \
included when the cutoff is set to zero.\n");
               }
            }
            else
            {
               printf("   Atom pairs will be included regardless\
 of interatomic distance\n");
            }
         }
         break;
      case 42: /* OCCRANK                                               */
         if((int)gNumParam[0] >= 1)
            gOccRank = (int)gNumParam[0];
         else 
            printf(   "Error: Occupancy rank must be >= 1\n");
         break;
      case 43: /* BZONE                                                 */
         SetZoneFromBValCol();
         break;
      case 44: /* SYMMATOMS                                             */
         /* Set Status                                                  */
         if(NParam)
         {
            if(!upstrncmp(gStrParam[0],"ON", 2) || 
               !upstrncmp(gStrParam[1],"ON", 2) ||
               !upstrncmp(gStrParam[0],"ALL",3))
               gMatchSymAtoms = TRUE;
            
            if(!upstrncmp(gStrParam[0],"OFF",3) ||
               (!upstrncmp(gStrParam[0],"ALL",3) &&
                !upstrncmp(gStrParam[1],"OFF",3)))
               gMatchSymAtoms = FALSE;
            
            for(i=0; i < SYMM_ATM_PAIRS; i++)
            {
               if(!upstrncmp(gStrParam[0],gSymType[i][0],3) ||
                  !upstrncmp(gStrParam[0],"ALL",3))
               {
                  if(!upstrncmp(gStrParam[1],"OFF",3))
                     gSymType[i][3][0] = FALSE;
                  else
                     gSymType[i][3][0] = TRUE;
               }
            }
         }
         
         /* Print Status */
         if(gMatchSymAtoms)
            printf("   Match Symmetric Atoms is ON\n");
         else 
            printf("   Match Symmetric Atoms is OFF\n");
         
         printf("   Atom pairs matched:\n");
         
         for(i=0; i<SYMM_ATM_PAIRS; i++)
         {
            printf("     %s%s -%s ",gSymType[i][0],
                   gSymType[i][1],gSymType[i][2]);
            
            if(gSymType[i][3][0] == TRUE)
               printf("  ON\n");
            else 
               printf(" OFF\n");
         }
         
         break;
      case 45: /* PRINTALIGN                                            */
         /* Error if no user fit zones                                  */
         if(!gUserFitZone)
         {
            printf("   Error: No user-defined zones found.\n");
            break;
         }
         
         /* Parameter FASTA prints a FASTA-formatted pairwise alignment */
         if(NParam>=1 && !upstrncmp(gStrParam[0],"FASTA",5) && 
            strlen(gStrParam[0]) == 5)
         {
            AlignmentFromZones((NParam == 2)?gStrParam[1]:NULL, TRUE);
         }
         else if(NParam>=1 && !upstrncmp(gStrParam[0],"PIR",3) && 
                 strlen(gStrParam[0]) == 3)
         {
            AlignmentFromZones_PIR((NParam == 2)?gStrParam[1]:NULL);
         }
         else 
         {
            AlignmentFromZones(NParam?gStrParam[0]:NULL, FALSE);    
         } 
         break;
      case 46: /* MULTREF                                               */
         if(NParam==1 && !upstrncmp(gStrParam[0],"OFF",3))
         {
            gMultiVsRef = FALSE;
            /* Reset reference to current multistructure reference */
            SetMobileToReference(gMultiRef);
            printf("   Multi: RMS, RESIDUE, PARDIST and MATRIX \
commands\n");
            printf("          compare with mobile structure: %d.\n",
                   gMultiRef+1);
            
         }
         else 
         {
            gMultiVsRef = TRUE;
            printf("   Multi: RMS, RESIDUE, PARDIST and MATRIX \
commands\n");
            printf("          compare with Averaged Reference \
structure.\n");
          } 
         break;
      case 47: /* ALLVSALL                                              */
         if(gMultiCount > 1 && !gIterate)
         {
            AllVsAllRMS((NParam?gStrParam[0]:NULL),
                        (NParam?TRUE:!gQuiet),FALSE);
         }
         else
         {
            /* Error Messages                                           */
            if(gMultiCount < 2)
            {
               printf("   Error: All vs all comparison ");
               printf("can only be used with MULTI.\n");
            }
            if(gIterate)
            {
               printf("   Error: All vs all comparison ");
               printf("cannot be used with iterative zones.\n");
            }
         }
        break;
      case 48: /* SETREF                                                */
         if(gMultiCount > 1)
         {
            if(NParam==0)
            {
               if(!gIterate)
               {
                  /* Automated set                                      */
                  AllVsAllRMS(NULL,FALSE,TRUE);
               }
               else
               {
                  /* Error: Can only use with MULTI                     */
                  printf("   Error: Automated reference selection ");
                  printf("cannot be used with iterative zones.\n");       
                  
               }
            }
            else 
            {
               if((int)gNumParam[0] >  0 && 
                  (int)gNumParam[0] <= gMultiCount)
               {
                  /* Manual set                                         */
                  SetMobileToReference((int)gNumParam[0] - 1);
                  gMultiRef = (int)gNumParam[0] - 1;
                  if(!gQuiet)
                     printf("   Reference set to mobile %d\n", 
                            (int)gNumParam[0]);
                  
               }
               else 
               {
                  /* Structure number out of range                      */
                  printf("   Error: ");
                  printf("Number entered must be between 1 and %d\n",
                         gMultiCount);
               }
            }
         }
         else 
         {
            /* Error: Can only use with MULTI */
            printf("   Error: ");
            printf("SETREF can only be used with MULTI.\n");        
         }
         
         break;
      case 49: /* ORDERFIT                                              */
         FitStructuresWrapper();
         break;
      case 50: /* TRIMZONES                                             */
         if(gMultiCount > 1)
         {
            TrimZones();
         }
         else
         {
            printf("   Error: ");
            printf("TRIMZONES can only be used with MULTI.\n");
         }
         break;
      case 51: /* WTAVERAGE                                             */
         
         if(NParam==1 && !upstrncmp(gStrParam[0],"OFF",3))
         {
            gWtAverage = FALSE;
         }
         else if((NParam==1 && !upstrncmp(gStrParam[0],"ON",3)) ||
                 (NParam==0))
         {
            gWtAverage = TRUE;
         }
         break;
      case 52: /* NOFIT                                                 */
         NoFitStructures();
         break;
      default:
         break;
      }
   }
   return(0);
}


/************************************************************************/
/*>void SetFitAtoms(char *command)
   -------------------------------
   This splits up a list of atoms names and inserts them in the 
   fitatoms list.
   
   28.09.92 Framework
   29.09.92 Original
   30.09.92 Added NOT option; pad fitatoms[]
   17.07.95 Replaced calls to screen() with printf()
   23.07.96 Changed to call PADMINTERM() rather than padterm()
            Checks for legal wildcard specifications
            Improved out-of-bounds error checking
   13.03.09 Added temporary fix for atom type recognition.
*/
void SetFitAtoms(char *command)
{
   int   comptr,
         comstart,
         atmnum,
         atmpos;
   char  *cmnd;
   
   /* Put to upper case                                                 */
   UPPER(command);
   
   /* Assume this is not a NOT selection                                */
   gNOTFitAtoms = FALSE;
   cmnd = command;
   
   /* Set NOT flag if required and step over the symbol                 */
   if(command[0] == '~' || command[0] == '^')
   {
      gNOTFitAtoms = TRUE;
      cmnd++;
   }

   /* Blank all fitatoms                                                */
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
      for(atmpos=0; atmpos<8; atmpos++)
         gFitAtoms[atmnum][atmpos] = '\0';
   
   atmnum=0;
   atmpos=0;
   comstart=0;
   for(comptr=0;comptr<strlen(cmnd);comptr++)
   {
      if(cmnd[comptr]==',')
      {
         comstart=comptr+1;
         gFitAtoms[atmnum][atmpos]   = '\0';  /* Terminate string       */
         PADMINTERM(gFitAtoms[atmnum],4);
         if(++atmnum >= NUMTYPES)
         {
            if(!gQuiet)
            {
               printf("   Warning==> Too many atoms in specification. \
Only %d used\n",NUMTYPES);
            }
            break;
         }
         atmpos=0;
      }
      else
      {
         if(atmpos >= MAXATSPEC)
         {
            TERMAT(cmnd+comptr, ',');
            if(!gQuiet)
            {
               printf("   Warning==> Atom name specification too long \
(%s)\n", cmnd+comstart);
               printf("              Using all atoms.\n");
            }
            gFitAtoms[0][0] = '*';
            gFitAtoms[0][1] = '\0';
            break;
         }

         gFitAtoms[atmnum][atmpos++] = cmnd[comptr];
      }
   }
   
   /* Terminate last one                                                */
   gFitAtoms[atmnum][atmpos] = '\0';
   PADMINTERM(gFitAtoms[atmnum],4);

   /* See if a * was specified not in first position; move it if so.
      Also check for legal atom wildcard specifications
   */
   for(atmnum=0; atmnum<NUMTYPES; atmnum++)
   {
      if(!LegalAtomSpec(gFitAtoms[atmnum]))
      {
         if(!gQuiet)
         {
            printf("   Warning==> Illegal atom specification (%s). Using \
all atoms\n",
                   gFitAtoms[atmnum]);
         }
         gFitAtoms[0][0] = '*';
         gFitAtoms[0][1] = '\0';
         break;
      }
      
      if(gFitAtoms[atmnum][0] == '*')
      {
         if(gNOTFitAtoms)
         {
            if(!gQuiet)
            {
               printf("   Warning==> NOT option ignored.\n");
            }
            gNOTFitAtoms = FALSE;
         }
         gFitAtoms[0][0] = '*';
         gFitAtoms[0][1] = '\0';
         break;
      }
   }
   
   gFitted = FALSE;
   

   /* Bugfix - Atom Recognition By: CTP */
   /* Temporary fix for atom type recognition. Removed trailing spaces
      for atom type strings. Need to check with ACRM why strings were
      padded. */
   if(TRUE)
   {
      int i = 0;
      for(i=0;i<NUMTYPES;i++)
      {
         if(gFitAtoms[i][0] == '\0') break;
         KILLTRAILSPACES(gFitAtoms[i]);
      }
   }
   /* End of Bugfix */


   return;
}


/************************************************************************/
/*>void SetFitZone(char *command, int strucnum)
   --------------------------------------------
   Processes commands to specify the zone for fitting.
   
   28.09.92 Original
   29.09.92 Added gCurrentMode
   01.10.92 Added gUserFitZone flag
   17.07.95 Replaced calls to screen() with printf()
   18.07.95 Added initialisation of inserts in zones
   19.07.95 * on its own equivalent to CLEAR
   18.06.96 Replaced MODE_* with ZONE_MODE_*
   01.02.01 Added multi-structure support incl. strucnum param
   28.02.01 Error message from ParseZone() when trying to specify multiple
            zones for multi-structure now here instead of in ParseZone()
   11.04.08 Added warning for overlapping zones. By: CTP
   22.04.08 Added handling of lowercase chain and inserts.
*/
void SetFitZone(char *command, int strucnum)
{
   int   start1,  stop1,
         start2,  stop2,
         SeqZone, maxstruc,
         snum;
   char  chain1,  chain2,
         startinsert1, stopinsert1,
         startinsert2, stopinsert2;
   ZONE  *z;
   int   warned = 0;

   snum = strucnum;
   
   /* See if this is clearing the zones                                 */
   if(!gUserFitZone || !upstrncmp(command,"CLEAR",5) || 
      !strcmp(command,"*"))
   {
      if(strucnum > (-1))
      {
         if(gZoneList[snum] != NULL)
         {
            FREELIST(gZoneList[snum],ZONE);
            gZoneList[snum] = NULL;
         }
      }
      else
      {
         for(snum=0; snum<gMultiCount; snum++)
         {
            if(gZoneList[snum] != NULL)
            {
               FREELIST(gZoneList[snum],ZONE);
               gZoneList[snum] = NULL;
            }
         }
      }
      gUserFitZone = FALSE;
      if(!upstrncmp(command,"CLEAR",5) || !strcmp(command,"*"))  
         return;
   }

   /* If strucnum is -1 then set range so we will do all structures     */
   if(strucnum == (-1))
   {
      snum = 0;
      maxstruc = gMultiCount;
   }
   else  /* Just do the one structure specified                         */
   {
      maxstruc = strucnum+1;
   }

   for(; snum<maxstruc; snum++)
   {
      SeqZone = ParseZone(command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1,
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2, snum);
      if((SeqZone == (-2)) && (!warned))
      {
         printf("   Error==> You cannot specify zones for each \
structure when performing\n");
         printf("            multiple structure fitting.\n");
         warned = 1;
      }
      
      if(SeqZone > -1)
      {
         /* Allocate an entry in zone list                              */
         if(gZoneList[snum])
         {
            /* Move to end of zone list                                 */
            z=gZoneList[snum];
            LAST(z);
            ALLOCNEXT(z,ZONE);
         }
         else
         {
            INIT(gZoneList[snum],ZONE);
            z = gZoneList[snum];
         }
         
         if(z==NULL)
         {
            printf("   Error==> No memory for zone!\n");
         }
         else
         {
            /* Add this zone to the zone list                           */
            z->chain1       = chain1;
            z->start1       = start1;
            z->startinsert1 = startinsert1;
            z->stop1        = stop1;
            z->stopinsert1  = stopinsert1;
            z->chain2       = chain2;
            z->start2       = start2;
            z->startinsert2 = startinsert2;
            z->stop2        = stop2;
            z->stopinsert2  = stopinsert2;
            z->mode         = SeqZone ? 
                              ZONE_MODE_SEQUENTIAL : gCurrentMode;

            /* CTP: Check for overlap                                   */
            if(CheckOverlap(z,gZoneList[snum],snum) > 1)
              printf("   Warning: New zone overlaps existing zone.\n");
            else if(CheckOverlap(z,gZoneList[snum],snum) == -1)
              printf("   Error: Failed to find new zone.\n");
            
         }
      }
   }
   
   gFitted      = FALSE;
   gUserFitZone = TRUE;
   
   return;
}

/************************************************************************/
/*>void DelFitZone(char *command, int strucnum)
   --------------------------------------------
   Processes commands to delete a zone for fitting.
   
   14.03.08 Original based on SetFitZone() By: CTP
   22.04.08 Added handling of lowercase chain and inserts.
   11.02.09 Changed CLEAR to ALL.
*/
void DelFitZone(char *command, int strucnum)
{
   int  start1,  stop1;
   int  start2,  stop2;
   int  SeqZone, maxstruc;
   int  snum;
   char chain1,  chain2;
   char startinsert1, stopinsert1;
   char startinsert2, stopinsert2;
   ZONE *z;
   int  warned = 0;

   snum = strucnum;
  
   /* Clear All Zones                                                   */
   if(!gUserFitZone || !upstrncmp(command,"ALL",3) || 
      !strcmp(command,"*"))
   {
      if(strucnum > (-1))
      {
         if(gZoneList[snum] != NULL)
         {
            FREELIST(gZoneList[snum],ZONE);
            gZoneList[snum] = NULL;
         }
      }
      else
      {
         for(snum=0; snum<gMultiCount; snum++)
         {
            if(gZoneList[snum] != NULL)
            {
               FREELIST(gZoneList[snum],ZONE);
               gZoneList[snum] = NULL;
            }
         }
      }
      gUserFitZone = FALSE;
      if(!upstrncmp(command,"ALL",3) || !strcmp(command,"*"))
      {
         printf("   All zones cleared.\n");
         return;
      }
   }
   
   /* If strucnum is -1 then set range so we will do all structures     */
   if(strucnum == (-1))
   {
      snum = 0;
      maxstruc = gMultiCount;
   }
   else  /* Just do the one structure specified                         */
   {
      maxstruc = strucnum+1;
   }

   /* Main Loop: Find and Delete Zone                                   */
   for(;snum<maxstruc; snum++)
   {
      /* Parse Zone                                                     */
      SeqZone = ParseZone(command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1,
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2, snum);
      
      if((SeqZone == (-2)) && (!warned))
      {
         printf("   Error==> You cannot specify zones for each \
structure when performing\n");
         printf("            multiple structure fitting.\n");
         warned = 1;
      }

      /* Compare with Zone list                                         */
      z=gZoneList[snum]; 

      while(z!=NULL)
      {
         /* Compare with current zone                                   */
         if((z->chain1       == chain1) &&
            (z->start1       == start1) &&
            (z->startinsert1 == startinsert1) &&
            (z->stop1        == stop1) &&
            (z->stopinsert1  == stopinsert1) &&
            (z->chain2       == chain2) &&
            (z->start2       == start2) &&
            (z->startinsert2 == startinsert2) &&
            (z->stop2        == stop2) &&
            (z->stopinsert2  == stopinsert2) &&
            (z->mode         == SeqZone ? 
             ZONE_MODE_SEQUENTIAL : gCurrentMode))
         {
            /* Delete Zone and Return                                   */
            DELETE(gZoneList[snum],z,ZONE);
            gFitted = FALSE;
            if(gZoneList[snum]!=NULL) 
               gUserFitZone = TRUE;
            else 
               gUserFitZone = FALSE;
            printf("   Zone Deleted\n");
            return;
         }
         NEXT(z);     
      }
   }
   /* Zone Not Found                                                    */
   printf("   No matching zone found.\n");
   return;
}


/************************************************************************/
/*>void DelRMSZone(char *command, int strucnum)
   --------------------------------------------
   Processes commands to delete a user-defined RMSd comparison zone.
   
   17.03.08 Original based on DelFitZone() By: CTP
   22.04.08 Added handling of lowercase chain and inserts.
   11.02.09 Changed CLEAR to ALL.
*/
void DelRMSZone(char *command, int strucnum)
{
   int  start1,  stop1;
   int  start2,  stop2;
   int  SeqZone, maxstruc;
   int  snum;
   char chain1,  chain2;
   char startinsert1, stopinsert1;
   char startinsert2, stopinsert2;
   ZONE *z;
   int  warned = 0;
   
   snum = strucnum;

   /* Clear All Zones                                                   */
   if(!upstrncmp(command,"ALL",3) || !strcmp(command,"*"))
   {
      if(strucnum > (-1))
      {
         if(gRZoneList[snum] != NULL)
         {
            FREELIST(gRZoneList[snum],ZONE);
            gRZoneList[snum] = NULL;
         }
      }
      else
      {
         for(snum=0; snum<gMultiCount; snum++)
         {
            if(gRZoneList[snum] != NULL)
            {
               FREELIST(gRZoneList[snum],ZONE);
               gRZoneList[snum] = NULL;
            }
         }
      }
      gUserRMSZone = TRUE;
      printf("   All zones cleared.\n");
      return;
   }

   /* If strucnum is -1 then set range so we will do all structures     */
   if(strucnum == (-1))
   {
      snum = 0;
      maxstruc = gMultiCount;
   }
   else  /* Just do the one structure specified                         */
   {
      maxstruc = strucnum+1;
   }
   
   /* Main Loop: Find and Delete RMSd Calc Zone                         */
   for(;snum<maxstruc; snum++)
   {
      /* Parse Zone                                                     */
      SeqZone = ParseZone(command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1,
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2, snum);
      
      if((SeqZone == (-2)) && (!warned))
      {
         printf("   Error==> You cannot specify zones for each \
structure when performing\n");
         printf("            multiple structure fitting.\n");
         warned = 1;
      }
      
      /* Compare with Zone list                                         */
      z=gRZoneList[snum]; 
      
      while(z!=NULL)
      {
         /* Compare to current zone                                     */
         if((z->chain1       == chain1) &&
            (z->start1       == start1) &&
            (z->startinsert1 == startinsert1) &&
            (z->stop1        == stop1) &&
            (z->stopinsert1  == stopinsert1) &&
            (z->chain2       == chain2) &&
            (z->start2       == start2) &&
            (z->startinsert2 == startinsert2) &&
            (z->stop2        == stop2) &&
            (z->stopinsert2  == stopinsert2) &&
            (z->mode         == SeqZone ?
             ZONE_MODE_SEQUENTIAL : gCurrentMode))
         {
            /* Delete Zone and Return                                   */
            DELETE(gRZoneList[snum],z,ZONE);
            gUserRMSZone = TRUE;
            printf("   Zone Deleted\n");
            return;
         }
         NEXT(z);     
      }
   }
   /* Zone Not Found                                                    */
   printf("   No matching zone found.\n");
   return;
}


/************************************************************************/
/*>void SetZoneStatus(char *status)
   --------------------------------
   Processes the command to set the mode for interpretation of residue
   numbers.
   
   28.09.92 Framework
   29.09.92 Original
   08.10.93 Modified only to toupper lowercase letters
   17.07.95 Replaced calls to screen() with printf()
   18.06.96 Replaced MODE_* with ZONE_MODE_*
*/
void SetZoneStatus(char *status)
{
   char ch;

   ch = status[0];
   if(islower(ch)) ch = toupper(ch);

   if(ch == 'R')
      gCurrentMode = ZONE_MODE_RESNUM;
   else if(ch == 'S') 
      gCurrentMode = ZONE_MODE_SEQUENTIAL;
   else
      printf("   Error==> Invalid numbering mode. Must be RESIDUE or \
SEQUENTIAL\n");

   return;
}


/************************************************************************/
/*>int GetResSpec(char *resspec, int *resnum, char *chain, char *insert)
   ---------------------------------------------------------------------
   Extracts residue number and chain name from a string. If the string is
   blank, the residue number will be specified as -999 and the chain name
   will be unmodified.

   29.09.92 Original
   17.07.95 Returns insert code. Uses ParseResSpec()
   20.02.01 Now makes its own copy of the string with an \ characters
            removed
   20.02.01 -999 for start or end of structure rather than -1
   22.04.08 Set to use ParseResSpecNoUpper() allowing function to deal 
            with lowercase chain names and inserts.  By: CTP
*/
int GetResSpec(char *resspec, int *resnum, char *chain, char *insert)
{
   char  *ptr, *buffer;
   int   i,
         retval = 0;

   if((buffer = (char *)malloc(strlen(resspec)+1))==NULL)
   {
      printf("   Error==> no memory for copy of residue spec\n");
      return(2);
   }
   
   /* Move pointer over any spaces                                      */
   for(ptr=resspec; *ptr==' '||*ptr=='\t'; ptr++) ;

   /* Copy the resspec string skipping \ characters                     */
   for(i=0; *ptr; ptr++)
   {
      if(*ptr != '\\')
         buffer[i++] = *ptr;
   }
   buffer[i] = '\0';

   /* Assume blank insert code                                          */
   *insert = ' ';
   
   /* If the zone is blank, it will be all residues                     */
   if(*buffer == '\0')
   {
      *resnum = -999;
      free(buffer);
      return(0);
   }

   if(!ParseResSpecNoUpper(buffer, chain, resnum, insert))
      retval = 1;              /* Indicates an error                    */

   free(buffer);
   return(retval);
}


/************************************************************************/
/*>int ParseZone(char *zonespec, 
                 int *start1, int *stop1, char *chain1, 
                 char *startinsert1, char *stopinsert1,
                 int *start2, int *stop2, char *chain2,
                 char *startinsert2, char *stopinsert2,
                 int strucnum)
   ----------------------------------------------------
   Returns: 0 for numeric zones,
            1 for sequence specified zones,
           -1 for error.
           -2 for error in multi-zone

   Sorts out zone specification. The specified residue ranges will be
   returned in start1/stop1/chain1 and start2/stop2/chain2. The first or 
   last residue present will be specified by -999. Zones may be specified 
   by number or by sequence.

   28.09.92 Original
   29.09.92 Added Sequence zones and return flag
   17.07.95 Replaced calls to screen() with printf()
            Added inserts
   01.02.01 Added strucnum parameter
   20.02.01 Allow negative residue numbers to be escaped with a \ 
            Code realises this - is not a range spcifier.
   20.02.01 -999 for start or end of structure rather than -1
   28.02.01 -2 return for multi-zone error - moved message out to calling
            routine
   22.04.08 Modified to alow use of lower case chain names. By: CTP
   16.05.08 Case sensitivity set after zonespec split into zone1 and zone2

*/
int  ParseZone(char *zonespec,
               int  *start1, 
               int  *stop1,
               char *chain1,
               char *startinsert1,
               char *stopinsert1,
               int  *start2, 
               int  *stop2,
               char *chain2,
               char *startinsert2,
               char *stopinsert2,
               int  strucnum)
{
   char *zone1,
        *zone2,
        *ptr,
        *dash = NULL,
        buffer[80];
   int  retval = 0;

   /* Blank the chain and insert names                                  */
   *chain1       = *chain2       = ' ';
   *startinsert1 = *startinsert2 = ' ';
   *stopinsert1  = *stopinsert2  = ' ';
   
   /* First split into 2 parts if required
      17.07.95 Added calls to KILLLEADSPACES()
   */
   KILLLEADSPACES(zone1, zonespec);
   zone2 = zone1;

   if((ptr=strchr(zonespec,':'))!=NULL)
   {
      /* 20.02.01 We don't allow this type of zone spec when we have 
         multiple structures
      */
      if(gMultiCount > 1)
      {
         return(-2);
      }
      
      KILLLEADSPACES(zone2, ptr+1);
      *ptr  = '\0';
   }

   /* Convert to uppercase unless a full stop is used as a separator. */
   /* A full stop is used as a separator with either a numeric or lower 
      case chain names. */
   if(strchr(zone1,'.')==NULL)
     UPPER(zone1);
   if(strchr(zone2,'.')==NULL)
     UPPER(zone2);
   
/*
*** Do zone1 first                                                     ***
*/

   /* See if there is a dash representing a range, skipping any - sign
      escaped with a \ (\- is used to represent a negative residue
      number)
   */
   dash=FindDash(zone1);

   /* If there's a * in the zone spec. it's all residues                */
   if((ptr=strchr(zone1,'*'))!=NULL)
   {
      *start1 = -999;
      *stop1  = -999;
      
      /* If * was not the first char in the zone spec., then the first
         character was the chain name.
      */
      if(ptr != zone1)
         *chain1 = *zone1;
   }
   else if(dash != NULL)                  /* - indicates a numeric zone */
   {
      /* Make a temporary copy                                          */
      strcpy(buffer,zone1);
      /* Terminate copy at the -                                        */
      *(FindDash(buffer)) = '\0';
      /* Read the start of the zone                                     */
      if(GetResSpec(buffer,start1,chain1,startinsert1))
      {
         printf("   Invalid zone specification: ");
         printf("%s\n",buffer);
         return(-1);
      }

      /* Read end of zone                                               */
      if(GetResSpec(dash+1,stop1,chain1,stopinsert1))
      {
         printf("   Invalid zone specification: ");
         printf("%s\n",dash+1);
         return(-1);
      }
   }
   else                            /* It's a sequence specified zone    */
   {
      if(FindSeq(zone1,gRefSeq,start1,stop1,chain1))
      {
         printf("   Sequence zone specification not found: ");
         printf("%s\n",zone1);
         return(-1);
      }
      else
      {
         retval = 1;
      }
   }

/*
*** Do zone2 if different from zone1                                   ***
*/
   if((zone1 == zone2) && retval==0)
   {
      *start2       = *start1;
      *stop2        = *stop1;
      *chain2       = *chain1;
      *startinsert2 = *startinsert1;
      *stopinsert2  = *stopinsert1;
   }
   else
   {
      dash=FindDash(zone2);

      /* If there's a * in the zone spec. it's all residues             */
      if((ptr=strchr(zone2,'*'))!=NULL)
      {
         *start2 = -999;
         *stop2  = -999;
         
         /* If * was not the first char in the zone spec., then the first
            character was the chain name.
         */
         if(ptr != zone2)
            *chain2 = *zone2;
      }
      else if(dash!=NULL)                    /* - shows a numeric zone  */
      {
         /* Make a temporary copy                                       */
         strcpy(buffer,zone2);
         /* Terminate copy at the -                                     */
         *(FindDash(buffer)) = '\0';
         /* Read the start of the zone                                  */
         if(GetResSpec(buffer,start2,chain2,startinsert2))
         {
            printf("   Invalid zone specification: ");
            printf("%s\n",buffer);
            return(-1);
         }
         
         /* Read end of zone                                            */
         if(GetResSpec(dash+1,stop2,chain2,stopinsert2))
         {
            printf("   Invalid zone specification: ");
            printf("%s\n",dash+1);
            return(-1);
         }
      }
      else                           /* It's a sequence specified zone  */
      {
         if(FindSeq(zone2,gMobSeq[strucnum],start2,stop2,chain2))
         {
            printf("   Sequence zone specification not found: ");
            printf("%s\n",zone2);
            return(-1);
         }
         else
         {
            retval = 1;
         }
      }
   }
   
   return(retval);
}


/************************************************************************/
/*>int FindSeq(char *zonespec, char *sequence, int *start, int *stop, 
               int *chain)
   ------------------------------------------------------------------
   Returns: 0 Sequence found
            1 Sequence not found

   Finds the start and stop on the basis of sequence.
            
   29.09.92 Original
   30.09.92 +1 Correction to length
*/
int FindSeq(char *zonespec,
            char *sequence, 
            int  *start, 
            int  *stop, 
            char *chain)
{
   char  zoneseq[40],
         *ptr;
   int   length,
         occurence,
         noccur,
         j;
   
   /* Extract the sequence part                                         */
   strcpy(zoneseq,zonespec);
   ptr = strchr(zoneseq,',');
   if(ptr!=NULL) *ptr = '\0';
   ptr = strchr(zoneseq,'/');
   if(ptr!=NULL) *ptr = '\0';
   
   /* Extract the length                                                */
   ptr = strchr(zonespec,',');
   if(ptr!=NULL)    sscanf(ptr+1,"%d",&length);
   else             length = strlen(zoneseq);
   
   /* Extract the occurence number                                      */
   ptr = strchr(zonespec,'/');
   if(ptr!=NULL)    sscanf(ptr+1,"%d",&occurence);
   else             occurence = 1;
   
   /* Now search for the occurence'th occurence of zoneseq              */
   noccur = 0;
   for(j=0; j<strlen(sequence)-length+1; j++)
   {
      if(!strncmp(zoneseq,sequence+j,strlen(zoneseq)))   /* 30.09.92    */
      {
         if(++noccur == occurence)
         {
            *start = j+1;
            *stop  = j+length;
            *chain = ' ';
            return(0);
         }
      }
   }

   *start = -2;
   *stop  = -2;
   *chain = 'X';

   return(1);
}


/************************************************************************/
/*>void ShowMatrix(void)
   ---------------------
   Displays the rotation matrix.
   
   28.09.92 Framework
   30.09.92 Original
   17.07.95 Replaced calls to screen() with printf()
   02.10.95 Added printing of CofGs
   20.02.01 gMobCofG now an array
   25.07.08 Modified to give translation and rotation vector to 
            gMobCofG[0] instead of the averaged reference structure, 
            gRefCofG. By: CTP
   11.11.08 Modified to give translation + rotation to gMultiRef.
   16.02.09 Modified rotation matrix inversion.
*/
void ShowMatrix(void)
{
   int   i,j,
         strucnum;

   /* Return if not fitted                                              */
   if(!gFitted)
   {
      printf("   Warning==> Structures have not yet been fitted.\n");
   }

   if(gMultiCount == 1)
   {
      /* Deal with single mobile structure.                             */
      /* ==================================                             */
      
      printf("   Reference CofG...\n");
      printf("   %8.4f %8.4f %8.4f\n",gRefCofG.x,
                                      gRefCofG.y,
                                      gRefCofG.z);
      printf("   Mobile CofG...\n");
      printf("   %8.4f %8.4f %8.4f\n",gMobCofG[0].x,
                                      gMobCofG[0].y,
                                      gMobCofG[0].z);
      
      printf("   Rotation matrix...\n");
      for(i=0; i<3; i++)
      {
         printf("   %8.4f %8.4f %8.4f\n",gRotMat[0][i][0],
                                         gRotMat[0][i][1],
                                         gRotMat[0][i][2]);
      } 

      printf("   Translation vector (between CofGs)...\n");
      printf("   %8.4f %8.4f %8.4f\n",gRefCofG.x - gMobCofG[0].x,
                                      gRefCofG.y - gMobCofG[0].y,
                                      gRefCofG.z - gMobCofG[0].z);
   }
   else
   {
      /* Deal with multiple mobile stuctures.                           */
      /* ====================================                           */
      REAL InvRotMat[3][3], ModRotMat[3][3];

      /* Invert rotation matrix gRotMat[gMultiRef]                      */
      /* The inverse of a rotation matrix is its transpose so copy 
         gRotMat[gMultiRef] to InvRotMat while transposing elements.    */
      for(i=0;i<3;i++)
      {
         for(j=0;j<3;j++)
         {
            InvRotMat[j][i] = gRotMat[gMultiRef][i][j];
         }
      }
            
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         printf("   Structure %d CofG...\n", strucnum+1);
         
         printf("   %8.4f %8.4f %8.4f\n",gMobCofG[strucnum].x,
                                         gMobCofG[strucnum].y,
                                         gMobCofG[strucnum].z);
      }
           
      printf("   Rotation matrix...\n");
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         /* Calculate rotation matrix                                   */
         MatMult33_33(gRotMat[strucnum],InvRotMat,ModRotMat);  
         printf("   Structure %d:\n", strucnum+1);

         if(gMultiVsRef)
         {
            for(i=0; i<3; i++)
            {
               printf("   %8.4f %8.4f %8.4f\n",
                      gRotMat[strucnum][i][0],
                      gRotMat[strucnum][i][1],
                      gRotMat[strucnum][i][2]);
            }
         }
         else 
         {
            for(i=0; i<3; i++)
            {
               printf("   %8.4f %8.4f %8.4f\n",
                      ModRotMat[i][0],
                      ModRotMat[i][1],
                      ModRotMat[i][2]);
            }
         }
         
      }

      printf("   Translation vector (between CofGs)...\n");
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         printf("   Structure %d:\n", strucnum+1);
         
         if(gMultiVsRef)
         {
            printf("   %8.4f %8.4f %8.4f\n",
                   gRefCofG.x - gMobCofG[strucnum].x,
                   gRefCofG.y - gMobCofG[strucnum].y,
                   gRefCofG.z - gMobCofG[strucnum].z);
         }
         else 
         {
            printf("   %8.4f %8.4f %8.4f\n",
                   gMobCofG[gMultiRef].x - gMobCofG[strucnum].x,
                   gMobCofG[gMultiRef].y - gMobCofG[strucnum].y,
                   gMobCofG[gMultiRef].z - gMobCofG[strucnum].z);
         }
      }
   }

   return;
}


/************************************************************************/
/*>void ShowStatus(char *filename)
   -------------------------------
   Shows the current program status.
   
   28.09.92 Framework
   29.09.92 Original
   30.09.92 Modified for padding of atom names
   01.10.92 Modification to All zone printing
   17.07.95 Replaced calls to screen() with printf()
   18.07.95 Prints zones with chain first and with inserts
            Prints *'s rather than -1 for all residues
            Added calls to FormatZone()
   20.07.95 Added WEIGHTS
   25.07.95 Added code to print chain labels
   24.01.96 Fixed bug in printing atom names containing spaces
   31.06.96 Added BValue cutoff 
   13.06.96 Modified B-value weighting message
   18.06.96 Replaced MODE_* with ZONE_MODE_*
   11.11.96 BValue cutoff message reflects new REF and MOB 
            parameters
   12.01.01 gMobPDB[] now an array
            Added iterate mode printing
   01.02.01 Added printing of multi structure data
   20.02.01 -999 for start or end of structure rather than -1
   28.03.01 Reports CENTRE mode
   20.03.08 Reports fit centred on residue. By: CTP
   07.04.08 Distance cutoff included.
   10.06.08 Trying auto-setting zones to current numbering mode. 
   01.08.08 Auto-setting zones feature disabled. 
   05.09.08 Added reporting of gap penalties.
   07.11.08 Altered option to output to a file.
            Marked reference for multi fitting.
*/
void ShowStatus(char *filename)
{
   char  buffer[240],
         atm[8],
         *chains;
   int   i, j, strucnum;
   ZONE  *z;
   FILE   *fp = stdout;

   /* Convert zones to current numbering mode                           */
   /* Feature disabled - 01.08.08
      This will put all zones into the current residue numbering scheme.
      However, this is slow and may not be what users want, so has been
      disabled. It hasn't been tested for some while, so re-enable at your
      own risk!
   */
/***
   SortAllZones();
   ConvertAllZones(gCurrentMode);
***/

   /* Open output file/pipe                                             */
   if(filename)
   {
      if((fp=OpenOrPipe(filename))==NULL)
      {
         printf("   Warning: Unable to open output file\n");
         fp = stdout;
      }
   }

   fprintf(fp,"   Reference structure:        %s\n",
           (gRefFilename[0]?gRefFilename:"Undefined"));

   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      /* Mark reference                                                 */
      if(gMultiCount > 1 && strucnum == gMultiRef)
         fprintf(fp,"  >");
      else 
         fprintf(fp,"   ");
      if(gMultiCount == 1)
      {
         fprintf(fp,"Mobile structure:           %s\n",
                 (gMobFilename[strucnum][0] ?
                  gMobFilename[strucnum]:"Undefined"));
      }
      else
      {
         fprintf(fp,"Mobile structure:%4d       %s\n",
                 strucnum + 1,
                 (gMobFilename[strucnum][0] ?
                  gMobFilename[strucnum]:"Undefined"));
      }
   }

   fprintf(fp,"   HETATM records are:         %s\n",
           (gHetAtoms?"Included":"Ignored"));

   fprintf(fp,"   Align gap penalty:          %2d\n",gGapPen);
   fprintf(fp,"   Align gap extend penalty:   %2d\n",gGapPenExt);

   fprintf(fp,"   Fitting will be:            %s\n",
           (gDoWeights==WEIGHT_BVAL?"Weighted by B-value":
            (gDoWeights==WEIGHT_INVBVAL?"Weighted by 1/B-value":
             "Normal (unweighted)")));

   if(((chains = GetPDBChainLabels(gRefPDB)) != NULL) &&
      ((strlen(chains) > 1) || (chains[0] != ' ')))
      fprintf(fp,"   Reference structure Chains: %s\n",chains);
   if(chains != NULL)
      free(chains);
   
   if(((chains = GetPDBChainLabels(gMobPDB[0])) != NULL) &&
      ((strlen(chains) > 1) || (chains[0] != ' ')))
      fprintf(fp,"   Mobile structure Chains:    %s\n",chains);
   if(chains != NULL)
      free(chains);

   fprintf(fp,"   Current numbering mode:     ");
   if(gCurrentMode == ZONE_MODE_RESNUM)
      fprintf(fp,"Residue\n");
   else if(gCurrentMode == ZONE_MODE_SEQUENTIAL)
      fprintf(fp,"Sequential\n");
   
   fprintf(fp,"   Iterative zone updating:    ");
   fprintf(fp,"%s\n", ((gIterate)?"On":"Off"));

   fprintf(fp,"   Atoms being fitted:         ");
   if(gFitAtoms[0][0] == '*')
   {
      fprintf(fp,"All\n");
   }
   else
   {
      if(gNOTFitAtoms) fprintf(fp,"NOT ");
      
      for(i=0;i<NUMTYPES;i++)
      {
         if(gFitAtoms[i][0] == '\0') break;
         strcpy(atm,gFitAtoms[i]);
         /* Remove trailing spaces                                      */
         /* 24.01.96 Was terminating at the first space. This broke
            atom names containing spaces (e.g. "N A" in heme groups)
         */
         for(j=strlen(atm)-1; j>=0; j--)
         {
            if(atm[j] == ' ')
               atm[j] = '\0';
            else
               break;
         }

         if(i)
            sprintf(buffer, ", %s", atm);
         else
            sprintf(buffer, "%s",   atm);
         fprintf(fp,buffer);
      }
      fprintf(fp,"\n");
   }
   if(gUseBVal)
   {
      fprintf(fp,"   Atoms will be discarded if their B-value is > %.2f \
in %s structure\n", gBValue, 
              (gUseBVal==2?"the reference":
               (gUseBVal==3?"the mobile":
                "either")));
   }
   else
   {
      fprintf(fp,"   Atoms will be included regardless of B-value\n");
   }
      
   if(gUseDistCutoff)
   {
      fprintf(fp,"   Atom pairs will be discarded if their interatomic \
distance is > %.2f\n", gDistCutoff); 
   }
   else
   {
      fprintf(fp,"   Atom pairs will be included regardless of \
interatomic distance\n");
   }
   
   fprintf(fp,"   Zones being fitted:         ");

   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      if(gMultiCount > 2)
         fprintf(fp,"\n   (Mobile Structure: %d)   ", strucnum+1);

      if(gZoneList[strucnum]==NULL || !gUserFitZone)
      {
         fprintf(fp,"All\n");
      }
      else
      {
         int overlap = 0;

         fprintf(fp,"\n");
         for(z=gZoneList[strucnum]; z!=NULL; NEXT(z))
         {
            char zone1[64],
                 zone2[64];
            
            FormatZone(zone1, z->chain1, 
                       z->start1, z->startinsert1, 
                       z->stop1,  z->stopinsert1);
            
            FormatZone(zone2, z->chain2, 
                       z->start2, z->startinsert2, 
                       z->stop2,  z->stopinsert2);

            fprintf(fp,"      %-16s with %-16s %s",
                   zone1, zone2,
                   ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                    :"(Sequential numbering)"));

            /* CTP: Check for Overlap                                   */
            if(CheckOverlap(z,gZoneList[strucnum],strucnum) > 1)
            {
               overlap++;
               fprintf(fp,"*\n");
            }
            else 
            {
               fprintf(fp,"\n");
            }
         }
         if(overlap)
            fprintf(fp,"%44s * Overlapping zones.\n","");
      }
   
      /* CTP: Display centre for fitting                                */
      if(gCZoneList[strucnum] != NULL)
      {
         char res1[64], res2[64];
         
         sprintf(res1,"%c%d%c",
                 gCZoneList[strucnum]->chain1,
                 gCZoneList[strucnum]->start1,
                 gCZoneList[strucnum]->startinsert1);
         
         sprintf(res2,"%c%d%c",
                 gCZoneList[strucnum]->chain2,
                 gCZoneList[strucnum]->start2,
                 gCZoneList[strucnum]->startinsert2);
          
         fprintf(fp,"   Fit centred on residues:\n");
         fprintf(fp,"                %-6s   to           %-6s %s\n",
                 res1, res2,
                 ((gCZoneList[strucnum]->mode == ZONE_MODE_RESNUM)?
                  "(Residue numbering)":"(Sequential numbering)"));
      }
   }

   if(gFitted)  /* Only display this when its definitely valid          */
   {
      fprintf(fp,"   Atoms for RMS calculation:  ");
   
      if(gUserRMSAtoms)
      {
         if(gRMSAtoms[0][0] == '*')
         {
            fprintf(fp,"All\n");
         }
         else
         {
            if(gNOTRMSAtoms) fprintf(fp,"NOT ");

            for(i=0;i<NUMTYPES;i++)
            {
               if(gRMSAtoms[i][0] == '\0') break;
               strcpy(atm,gRMSAtoms[i]);
               if(strchr(atm,' ')) *strchr(atm,' ') = '\0';
               
               if(i)
                  sprintf(buffer, ", %s", atm);
               else
                  sprintf(buffer, "%s",   atm);
               fprintf(fp,buffer);
            }
            fprintf(fp,"\n");
         }
      }
      else
      {
         if(gFitAtoms[0][0] == '*')
         {
            fprintf(fp,"All\n");
         }
         else
         {
            if(gNOTFitAtoms) fprintf(fp,"NOT ");

            for(i=0;i<NUMTYPES;i++)
            {
               if(gFitAtoms[i][0] == '\0') break;
               strcpy(atm,gFitAtoms[i]);
               if(strchr(atm,' ')) *strchr(atm,' ') = '\0';

               if(i)
                  sprintf(buffer, ", %s", atm);
               else
                  sprintf(buffer, "%s",   atm);
               fprintf(fp,buffer);
            }
            fprintf(fp,"\n");
         }
      }
      
      fprintf(fp,"   Zones for RMS calculation:  ");
   
      if(gUserRMSZone)
      {
         for(strucnum=0; strucnum<gMultiCount; strucnum++)
         {
            int overlap = 0;
            if(gMultiCount > 2)
               fprintf(fp,"\n   (Mobile Structure: %d)   ", strucnum+1);

            if(gRZoneList[strucnum]==NULL) fprintf(fp,"All\n");
            else                fprintf(fp,"\n");
            
            for(z=gRZoneList[strucnum]; z!=NULL; NEXT(z))
            {
               char zone1[64],
                    zone2[64];
               
               FormatZone(zone1, z->chain1, 
                          z->start1, z->startinsert1, 
                          z->stop1,  z->stopinsert1);
               
               FormatZone(zone2, z->chain2, 
                          z->start2, z->startinsert2, 
                          z->stop2,  z->stopinsert2);
               
               fprintf(fp,"      %-16s with %-16s %s",
                      zone1, zone2,
                      ((z->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
                       :"(Sequential numbering)"));
               if(CheckOverlap(z,gRZoneList[strucnum],strucnum) > 1)
                 {
                   overlap++;
                   fprintf(fp,"*\n");
                 }
               else 
                 {
                   fprintf(fp,"\n");
                 }
            }
            if(overlap)
              fprintf(fp,"%44s * Overlapping zones.\n","");
         }
      }
      else
      {
         for(strucnum=0; strucnum<gMultiCount; strucnum++)
         {
            if(gMultiCount > 2)
               fprintf(fp,"\n   (Mobile Structure: %d)   ", strucnum+1);
            if(gZoneList[0]==NULL || !gUserFitZone)
            {
               fprintf(fp,"All\n");
            }
            else
            {
               int overlap = 0;
               fprintf(fp,"\n");
               
               for(z=gZoneList[strucnum]; z!=NULL; NEXT(z))
               {
                  char zone1[64],
                       zone2[64];
                  
                  FormatZone(zone1, z->chain1, 
                             z->start1, z->startinsert1, 
                             z->stop1,  z->stopinsert1);
                  
                  FormatZone(zone2, z->chain2, 
                             z->start2, z->startinsert2, 
                             z->stop2,  z->stopinsert2);
                  
                  fprintf(fp,"      %-16s with %-16s %s",
                         zone1, zone2,
                         ((z->mode == ZONE_MODE_RESNUM)?
                          "(Residue numbering)"
                          :"(Sequential numbering)"));

                  /* CTP: Check for overlap                             */
                  if(CheckOverlap(z,gZoneList[strucnum],strucnum) > 1)
                  {
                     overlap++;
                     fprintf(fp,"*\n");
                  }
                  else 
                  {
                     fprintf(fp,"\n");
                  }
               }
               if(overlap)
                  fprintf(fp,"%44s * Overlapping zones.\n","");
            }
         }
      }
   }
   
   if(gCentre)
   {
      fprintf(fp,"\n   Coordinates written centred on: ORIGIN\n\n");
   }
   else
   {
      fprintf(fp,"\n   Coordinates written centred on: REFERENCE \
SET\n\n");
   }
   
   fprintf(fp,"   Reference sequence:         ");
   if(gRefSeq)
   {
      for(i=0, j=0; i<strlen(gRefSeq); i++)
      {
         buffer[j++] = gRefSeq[i];
         
         if(j>59 || i==strlen(gRefSeq)-1)
         {
            buffer[j] = '\0';
            fprintf(fp,"\n   ");
            fprintf(fp,buffer);
            j=0;
         }
      }
      fprintf(fp,"\n");
   }
   else
   {
      fprintf(fp,"Undefined\n");
   }

   /* CTP                                                               */
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      if(gMultiCount == 1)
         fprintf(fp,"   Mobile sequence:            ");
      else
         fprintf(fp,"   Mobile sequence:%4d        ", strucnum+1);
      
      if(gMobSeq[strucnum])
      {
         for(i=0, j=0; i<strlen(gMobSeq[strucnum]); i++)
         {
            buffer[j++] = gMobSeq[strucnum][i];
            
            if(j>59 || i==strlen(gMobSeq[strucnum])-1)
            {
               buffer[j] = '\0';
               fprintf(fp,"\n   ");
               fprintf(fp,buffer);
               j=0;
            }
         }
         fprintf(fp,"\n");
      }
      else
      {
         fprintf(fp,"Undefined\n");
      }
   }
   
   /* Close output file/pipe                                            */
   if(fp != stdout)
      CloseOrPipe(fp);

   return;
}


/************************************************************************/
/*>void stdprompt(char *string)
   ----------------------------
   Issues a prompt to stdout providing stdin is a tty

   18.07.95 Original    By: ACRM
*/
#include <unistd.h>
void stdprompt(char *string)
{
#if (unix || __unix__)
   if(!isatty(0))
      return;
#endif
   
   printf("%s> ",string);
   fflush(stdout);
}


/************************************************************************/
/*>void FormatZone(char *zone, char chain, int start, char startinsert, 
                   int stop,  char stopinsert)
   --------------------------------------------------------------------
   Formats a zone specification accounting for -999 values which represent
   all residues

   18.07.95 Original    By: ACRM
   20.02.01 -999 for start or end of structure rather than -1
*/
void FormatZone(char *zone, char chain, int start, char startinsert, 
                int stop,  char stopinsert)
{
   char part1[16],
        part2[16];
   
   if((start == (-999)) && (stop  == (-999)))
   {
      if(chain == ' ')
         sprintf(zone,"All residues");
      else
         sprintf(zone,"Chain %c",chain);
   }
   else
   {
      if(start==(-999))
      {
         sprintf(part1,"%c*",    chain);
         sprintf(part2,"%c%d%c", chain, stop, stopinsert);
      }
      else if(stop==(-999))
      {
         sprintf(part1,"%c%d%c", chain, start, startinsert);
         sprintf(part2,"%c*",    chain);
      }
      else
      {
         sprintf(part1,"%c%d%c", chain, start, startinsert);
         sprintf(part2,"%c%d%c", chain, stop, stopinsert);
      }

      sprintf(zone,"%-6s to %-6s", part1, part2);
   }
}


/************************************************************************/
/*>void ReadMulti(char *filename, BOOL xmasFormat)
   -----------------------------------------------
   Reads in multiple structures as given by the MULTI command

   01.02.01 Original   By: ACRM
   01.03.01 Added xmasFormat parameter
   08.01.09 CTP modified error check for maximum number of structures.
*/
void ReadMulti(char *filename, BOOL xmasFormat)
{
   FILE *fof = NULL;
   BOOL GotRef = FALSE;
   char buffer[MAXBUFF], *ch;

   gMultiCount = 0;

   if((fof=fopen(filename,"r"))==NULL)
   {
      printf("   Error==> Can't open list of files: %s\n",
             filename);
   }
   else
   {
      while(fgets(buffer, MAXBUFF, fof))
      {
         TERMINATE(buffer);
         KILLTRAILSPACES(buffer);
         KILLLEADSPACES(ch, buffer);
         
         if((*ch != '#') && (*ch != '!') && (strlen(ch)))
         {
            if(!GotRef)
            {
               ReadStructure(STRUC_REFERENCE, ch, 0, xmasFormat);
               GotRef=TRUE;
            }
            
            /* CTP                                                      */
            if(gMultiCount < MAXSTRUC)
            {
               if(ReadStructure(STRUC_MOBILE,ch,gMultiCount,xmasFormat))
               {
                  gMultiCount++;
               }
            }
            else 
            {
               printf("   Error==> Maximum structure count (%d) \
exceeded. Increase MAXSTRUC\n", MAXSTRUC);
               printf("            Skipped structure: %s\n",ch);
            }
         }      
      }
   }
}


/************************************************************************/
/*>void WriteCoordinates(char *filename, int strucnum)
   ---------------------------------------------------
   Writes a coordinate file.

   28.09.92 Framework
   01.10.92 Original
   17.07.95 Replaced calls to screen() with printf()
   18.07.95 Uses fopen() rather than OpenWrite()
            Uses fclose() rather than CloseFile()
   21.07.95 Corrected logic of test for un-fitted PDB
   27.06.97 Changed call to fopen() to OpenOrPipe
   11.01.01 gFitPDB now an array
   01.02.01 Added strucnum parameter
   28.03.01 Added code to support strucnum < 0 to write the centred
            reference set and centering of coordinates
   02.04.08 added calls to WritePDBHeader() and WritePDBFooter() By: CTP
   02.05.08 changed WritePDBHeader() and WritePDBFooter() calls to 
            calls to WriteWholePDBHeader() and WriteWholePDBTrailer()
*/
void WriteCoordinates(char *filename, int strucnum)
{
   FILE  *fp;
   PDB   *pdb,
         *pdbc;
   VEC3F CofG;
   WHOLEPDB *wpdb;

   /* Point pdb to the coordinate set of interest                       */
   if(strucnum < 0)  /* use the reference set                           */
   {
      pdb  = gRefPDB;
      wpdb = gRefWPDB;
   }
   else              /* use a mobile set                                */
   {
      pdb  = gFitPDB[strucnum];
      wpdb = gMobWPDB[strucnum];
   }

   CofG.x = -1.0 * gRefCofG.x;
   CofG.y = -1.0 * gRefCofG.y;
   CofG.z = -1.0 * gRefCofG.z;
   
   /* If centering, then make a copy of the coordinates and translate
      to the origin
   */
   if(gCentre && (pdb!=NULL))
   {
      if((pdbc = DupePDB(pdb))==NULL)
      {
         printf("   Error==> No memory for translating the reference \
coordinates\n");
         printf("            Centred reference coordinates not \
written.\n");
         return;
      }

      /* Point pdb to this copy of the coordinates                      */
      pdb = pdbc;
      TranslatePDB(pdb, CofG);
   }
   
   if(!gFitted || pdb == NULL)
   {
      printf("   Error==> Fitting has not yet been performed.\n");
   }
   else
   {
      if((fp=OpenOrPipe(filename))==NULL)
      {
         printf("   Error==> Enable to open file for writing.\n");
      }
      else
      {
         printf("   Writing coordinates...\n");
         
         if(gReadHeader)
            WriteWholePDBHeader(fp, wpdb);

         WritePDB(fp, pdb);

         if(gReadHeader)
            WriteWholePDBTrailer(fp, wpdb);
         
         CloseOrPipe(fp);
      }
   }

   /* If centering, free the copy of the reference set                  */
   if(gCentre && pdb!=NULL)
   {
      FREELIST(pdb, PDB);
   }

   return;
}


/************************************************************************/
/*>void WriteMulti(char *ext)
   --------------------------
   Writes a set of fitted files from multi-fitting

   01.02.01 Original   By: ACRM
   20.02.01 Added more defensive checks on string lengths and handle case
            where input filename didn't contain a .
*/
void WriteMulti(char *ext)
{
   char filename[MAXBUFF];
   int  i,
        j;

   for(i=0; i<gMultiCount; i++)
   {
      strncpy(filename, gMobFilename[i], MAXBUFF);
      /* Work from the end of string to strip off the extension         */
      j = strlen(filename)-1;


      /* Step back until we find a . delimiting the extension           */
      while((filename[j] != '.') && (j>=0))
      {
         /* Break out if we find a / \ ] or : since we are in the path  */
         if(filename[j] == '/'  || 
            filename[j] == '\\' ||
            filename[j] == ']'  ||
            filename[j] == ':')
         {
            j=0;
            break;
         }
         j--;
      }
      
      /* There was a . in the filename, so truncate the string there    */
      if(j)
      {
         filename[j] = '\0';
      }

      /* Report error if filename too long                              */
      if((strlen(filename) + strlen(ext) + 1) >= MAXBUFF)
      {
         printf("   Error==> Filename too long to add new extension\n");
         printf("            Fitted file not written: %s\n", 
                gMobFilename[i]);
         return;
      }

      /* Append extension                                               */
      strcat(filename,".");
      strcat(filename, ext);

      /* Write the file                                                 */
      WriteCoordinates(filename, i);
   }
}


/************************************************************************/
/*>char *FindDash(char *buffer)
   ----------------------------
   Find a dash representing a range, skipping any - sign escaped with a 
   \ (\- is used to represent a negative residue number)

   20.02.01 Original   By: ACRM
*/
char *FindDash(char *buffer)
{
   char *dash = NULL;

   dash=strchr(buffer,'-');
   while((dash != NULL) && (dash > buffer) && (*(dash-1) == '\\'))
   {
      buffer = dash+1;
      dash=strchr(buffer,'-');
   }

   return(dash);
}


/************************************************************************/
/*>PDB *ReadXMAS(FILE *fp, int *natoms)
   ------------------------------------
   Like ReadPDB() but reads from an XMAS file instead of a PDB file

   01.03.01 Original   By: ACRM
*/
PDB *ReadXMAS(FILE *fp, int *natoms)
{
   return(doReadXMAS(fp, natoms, TRUE));
}


/************************************************************************/
/*>PDB *ReadXMASAtoms(FILE *fp, int *natoms)
   -----------------------------------------
   Like ReadPDBAtoms() but reads from an XMAS file instead of a PDB file

   01.03.01 Original   By: ACRM
*/
PDB *ReadXMASAtoms(FILE *fp, int *natoms)
{
   return(doReadXMAS(fp, natoms, FALSE));
}


/************************************************************************/
#ifdef USE_XMAS
#  define _LIBACRM_H          /* Stops libacrm.h being included        */
#  define FDWRAP_MAX_BUFF 1024
   typedef struct
   {
      int fd;
      char buffer[FDWRAP_MAX_BUFF];
      int  buffpos,
           maxbuff,
           eof,
           socket;
   }
   FDWRAP;

#  include "xmas.h"
#endif


/************************************************************************/
/*>PDB *doReadXMAS(FILE *fp, int *natoms, int readhet)
   ---------------------------------------------------
   Reads an XMAS file into a PDB linked list

   01.03.01 Original   By: ACRM
   15.03.01 Removed special call to ReadXmasData() - no longer needed as
            XMAS column index is now stored in the XMAS structure
*/
PDB *doReadXMAS(FILE *fp, int *natoms, int readhet)
{
   PDB  *pdb = NULL;

#ifdef USE_XMAS
   XMAS *xmas = NULL;
   PDB  *p;
   char  atnum[16], 
         atnam[16], 
         x[16], 
         y[16], 
         z[16], 
         occup[16], 
         bval[16], 
         resnam[16], 
         resnum[16], 
         chain[16], 
         type[16];

   *natoms = 0;
   
   /* Read the XMAS header                                              */
   if((xmas = ReadXmasHeader(fp))==NULL)
   {
      printf("   Error==> Couldn't read XMAS header: %s\n", gXMASError);
      return(NULL);
   }

   /* Read in the XMAS data                                             */
   if(!CacheXmasData(xmas))
   {
      printf("   Error==> Couldn't read XMAS data: %s\n", gXMASError);
      FreeXmasData(xmas);
      return(NULL);
   }
   
   /* Check the data contains atom records                              */
   if(!DoesXmasContain(xmas, "atoms"))
   {
      fprintf(stderr,"   Error==> XMAS file does not have ATOM \
records!\n");
      FreeXmasData(xmas);
      return(NULL);
   }
   
   /* All looks OK, read the ATOM data into a PDB linked list           */
   while(ReadXmasData(xmas, "atoms",
                      "atnum, atnam, x, y, z, occup, bval, resnam, \
                       resnum, chain, type",
                      atnum, atnam, x, y, z, occup, 
                      bval, resnam, resnum, chain, type))
   {
      if((!strcmp(type, "ATOM")) ||
         (!strcmp(type, "HETATM") && readhet))
      {
         
         /* Allocate memory in the PDB linked list                      */
         if(pdb==NULL)
         {
            INIT(pdb, PDB);
            p=pdb;
         }
         else
         {
            ALLOCNEXT(p, PDB);
         }
         if(p==NULL)
         {
            FREELIST(pdb, PDB);
            *natoms = 0;
            return(NULL);
         }
         
         sscanf(atnum, "%d", &(p->atnum));

         DEDOTIFY(atnam);
         strcat(atnam, " ");
         strcpy(p->atnam_raw, atnam);
         strcpy(p->atnam, FixAtomName(atnam));
         p->atnam[4] = '\0';
         sscanf(x, "%lf", &(p->x));
         sscanf(y, "%lf", &(p->y));
         sscanf(z, "%lf", &(p->z));
         sscanf(occup, "%lf", &(p->occ));
         sscanf(bval,  "%lf", &(p->bval));
         strncpy(p->resnam, resnam, 3);
         strcat(p->resnam, " ");
         
         DEDOTIFY(resnum);
         p->insert[0] = resnum[4];
         p->insert[1] = '\0';
         resnum[4] = '\0';
         sscanf(resnum,  "%d", &(p->resnum));
         
         p->chain[0] = chain[0];
         p->chain[1] = '\0';
         
         strcpy(p->junk, type);
         PADCHARMINTERM(p->junk, ' ', 6);

         (*natoms)++;
      }
   }

   /* Free the XMAS data                                                */
   FreeXmasData(xmas);
   
#endif   

   return(pdb);
}


/************************************************************************/
/*>void SetCentreResidue(char *command)
   ------------------------------------
   Sets a residue (a single-residue zone in a one zone list) as the 
   centre for fitting.
   As the function sets-up a list of zones, future developments of ProFit
   could be set up to define the centre for fitting around a list of 
   zones. (For instance, around an active site.)

   19.03.08 Original based on SetRMSZone()  By: CTP
   22.04.08 Added handling of lowercase chain and inserts.
*/
void SetCentreResidue(char *command)
{
   int   start1,  stop1,
         start2,  stop2,
         SeqZone, strucnum;
   char  chain1,  chain2,
         startinsert1, stopinsert1,
         startinsert2, stopinsert2;
   ZONE  *z;
   int   warned = 0;

   char *zone1, *zone2, *ptr;
   char zone_command[20];
   
   /* See if this is clearing the zones                                 */
   if(!upstrncmp(command,"CLEAR",5) || !strcmp(command,"*"))
   {
      for(strucnum=0; strucnum<gMultiCount; strucnum++)
      {
         if(gCZoneList[strucnum]!=NULL)
         {
            FREELIST(gRZoneList[strucnum],ZONE);
            gCZoneList[strucnum] = NULL;
         }
      }

      gFitted = FALSE;
      return;
   }

   /* Parse command here                                                */
   KILLLEADSPACES(zone1, command);
   zone2 = zone1;
   if((ptr=strchr(command,':'))!=NULL)
   {
      if(gMultiCount > 1)
      {
         printf("   Warning: Cannot define residue using ':' \
notation.\n");
         return;
      }
       
      KILLLEADSPACES(zone2, ptr+1);
      *ptr  = '\0';
   }
    
   strcpy(zone_command,zone1);  
   strcat(zone_command,"-");
   strcat(zone_command,zone1);

   if(zone1 != zone2)
   {
      strcat(zone_command,":");  
      strcat(zone_command,zone2);  
      strcat(zone_command,"-");
      strcat(zone_command,zone2);
   }

   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      SeqZone = ParseZone(zone_command, &start1, &stop1, &chain1, 
                          &startinsert1, &stopinsert1, 
                          &start2, &stop2, &chain2,
                          &startinsert2, &stopinsert2,
                          strucnum);

      if((SeqZone == (-2)) && (!warned))
      {
         printf("   Error==> You cannot specify zones for each \
structure when performing\n");
         printf("            multiple structure fitting.\n");
         warned = 1;
      }
      
      if(SeqZone > -1)
      {
         /* Blank Current List                                          */
         if(gCZoneList[strucnum]!=NULL)
         {
            FREELIST(gCZoneList[strucnum],ZONE);
            gCZoneList[strucnum] = NULL;
         }

         /* Allocate entry                                              */
         INIT(gCZoneList[strucnum],ZONE);
         z = gCZoneList[strucnum];
         
         if(z==NULL)
         {
            printf("   Error==> No memory for zone!\n");
         }
         else
         {
            /* Add this zone to the zone list                           */
            z->chain1       = chain1;
            z->start1       = start1;
            z->startinsert1 = startinsert1;
            z->stop1        = stop1;
            z->stopinsert1  = stopinsert1;
            z->chain2       = chain2;
            z->start2       = start2;
            z->startinsert2 = startinsert2;
            z->stop2        = stop2;
            z->stopinsert2  = stopinsert2;
            z->mode         = SeqZone?ZONE_MODE_SEQUENTIAL:gCurrentMode;
         }
         
         gFitted = FALSE;
      }
   }

   return;
}


/************************************************************************/
/*>int RunScript(char *command)
   ----------------------------
   Run a script file. 
   Opens file, runs DoCommandLoop() then closes file. Also limits number
   of scripts within scripts.

   27.03.08 Original By: CTP
   15.04.08 Quitting in a script now exits from ProFit instead of 
            returning to DoCommandLoop().
*/
int RunScript(char *command)
{
  const  int max_recur_level = 1000;
  static int recur_level = 0;
  FILE   *script = NULL;
  char   filename[MAXSTRLEN];


  /* Only allow limited level of scripts running scripts                */
  if(recur_level >= max_recur_level)
  {
     printf("   Error==> The maximum number of nested scripts is %d.\n",
            max_recur_level);
     return(1);
  }
  
  /* Copy global variable 'command' to local 'filename'                 */
  strcpy(filename, command);
  
  /* Open script file                                                   */
  if((script = fopen(filename,"r")) == NULL)
  {
     printf("   Error==> Failed to open: '%s'\n",filename);
     return(1);
  }
  
  /* Run Script File                                                    */
  recur_level++;
  if(!gQuiet)
     printf("\n   Starting script: '%s'\n",filename);
  
  if(!DoCommandLoop(script))
  {
     if(!gQuiet)
        printf("   Finished script: '%s'\n\n",filename);
  }
  else 
  {
     printf("   Error in script: '%s'\n\n",filename);
  }
  
  fclose(script);
  recur_level--;
  
  return(0);
}



/************************************************************************/
/*>int ConvertResidueToSequential(PDB *pdb, ZONE *res_zone, 
                                  ZONE *seq_zone)
   --------------------------------------------------------
   Converts residue-numbered zone into sequential-numbered zone.

   08.04.08 Original By: CTP
   10.06.08 Set chain id to space for sequential-numbered zone.
*/
int ConvertResidueToSequential(ZONE *input_zone, int strucnum)
{
   PDB  *p            = NULL;
   ZONE *output_zone  = NULL; 
   int  residue_count = 0;
   int  prev_residue  = 0;
   char prev_insert   = ' ';
   char prev_chain    = ' ';

   /* Zero start and stop for output_zone                               */
   INIT(output_zone,ZONE);
   output_zone->start1 = 0;
   output_zone->start2 = 0;
   output_zone->stop1  = 0;
   output_zone->stop2  = 0;
   
   /* Reference Residue                                                 */
   residue_count = 0;
   prev_residue  = 0;
   prev_insert   = ' ';
   prev_chain    = ' ';
   
   for(p=gRefPDB; p!=NULL; NEXT(p))
   {
      /* Update residue count                                           */
      if((p->resnum    != prev_residue) ||
         (p->insert[0] != prev_insert)  ||
         (p->chain[0]  != prev_chain))
      {
         residue_count++;
      }
      
      /* Check for Start                                                */
      /* Start residue number == -999                                   */
      if(((input_zone->chain1 == ' ') || 
          (input_zone->chain1 == p->chain[0])) &&
         (input_zone->start1 == -999 && output_zone->start1 == 0))
      {
         output_zone->start1 = residue_count;
      }
      /* Start residue number defined                                   */
      if(((input_zone->chain1 == ' ') || 
          (input_zone->chain1 == p->chain[0])) &&
         (input_zone->start1 == p->resnum) &&
         (input_zone->startinsert1 == p->insert[0]) &&
         (output_zone->start1 == 0))
      {
         output_zone->start1 = residue_count;
      }
      
      /* Check for Finish                                               */
      /* Stop residue number == -999                                    */
      if(((input_zone->chain1 == ' ') || 
          (input_zone->chain1 == p->chain[0])) &&
         (input_zone->stop1 == -999))
      {
         output_zone->stop1 = residue_count;
      }
      /* Stop residue number defined                                    */
      if(((input_zone->chain1 == ' ') || 
          (input_zone->chain1 == p->chain[0])) &&
         (input_zone->stop1 == p->resnum) &&
         (input_zone->stopinsert1 == p->insert[0]) &&
         (output_zone->stop1 == 0))
      {
         output_zone->stop1 = residue_count;
      }
      
      /* Update previous residue                                        */
      prev_residue = p->resnum;
      prev_insert  = p->insert[0];
      prev_chain   = p->chain[0];
   }
   
   /* Mobile Residue                                                    */
   residue_count = 0;
   prev_residue  = 0;
   prev_insert   = ' ';
   prev_chain    = ' ';
   
   for(p=gMobPDB[strucnum]; p!=NULL; NEXT(p))
   {
      /* Update residue count                                           */
      if((p->resnum    != prev_residue) ||
         (p->insert[0] != prev_insert)  ||
         (p->chain[0]  != prev_chain))
      {
         residue_count++;
      }
      
      /* Check for Start                                                */
      /* Start residue number == -999                                   */
      if(((input_zone->chain2 == ' ') || 
          (input_zone->chain2 == p->chain[0])) &&
         (input_zone->start2 == -999 && output_zone->start2 == 0))
      {
         output_zone->start2 = residue_count;
      }
      /* Start residue number defined                                   */
      if(((input_zone->chain2 == ' ') || 
          (input_zone->chain2 == p->chain[0])) &&
         (input_zone->start2 == p->resnum) &&
         (input_zone->startinsert2 == p->insert[0]) &&
         (output_zone->start2 == 0))
      {
         output_zone->start2 = residue_count;
      }
      
      /* Check for Finish                                               */
      /* Stop residue number == -999                                    */
      if(((input_zone->chain2 == ' ') || 
          (input_zone->chain2 == p->chain[0])) &&
         (input_zone->stop2 == -999))
      {
         output_zone->stop2 = residue_count;
      }
      /* Stop residue number defined                                    */
      if(((input_zone->chain2 == ' ') || 
          (input_zone->chain2 == p->chain[0])) &&
         (input_zone->stop2 == p->resnum) &&
         (input_zone->stopinsert2 == p->insert[0]) &&
         (output_zone->stop2 == 0))
      {
         output_zone->stop2 = residue_count;
      }
      
      /* Update previous residue                                        */
      prev_residue = p->resnum;
      prev_insert  = p->insert[0];
      prev_chain   = p->chain[0];
   }
   
   /* Start and Stop Assigned?                                          */
   if(!(output_zone->start1) || !(output_zone->stop1) ||
      !(output_zone->start2) || !(output_zone->stop2))
   {
      /* ZONE NOT FOUND                                                 */
      return(1);
   }
   
   /* Set Sequential Numbering                                          */
   input_zone->start1 = output_zone->start1;
   input_zone->stop1  = output_zone->stop1;
   input_zone->start2 = output_zone->start2;
   input_zone->stop2  = output_zone->stop2;
   input_zone->chain1 = ' ';
   input_zone->chain2 = ' ';
   input_zone->mode   = ZONE_MODE_SEQUENTIAL;
   
   /* Cleanup                                                           */
   FREELIST(output_zone,ZONE);
   return(0);
}


/************************************************************************/
/*>int CheckOverlap(ZONE *inputtest, ZONE *inputlist, int strucnum)
   ----------------------------------------------------------------
   Function returns number of pairs of zones overlapping or -1 if the 
   function cannot find the zone inputtest or the list of zones 
   inputlist. 

   This function is called by: 
      ShowStatus() when flagging overlapping zones.
      SetFitZone() and SetRMSZone() when checking new user defined zones
         for overlaps with existing zones.

   The function duplicates and converts residue-numbered zones into 
   sequentially-numbered zones for comparison.

   11.04.08 Original By: CTP
   14.04.08 Modified to make zonelist a list of two zones, A and B.
*/
int CheckOverlap(ZONE *inputtest, ZONE *inputlist, int strucnum)
{
  ZONE *zonelist = NULL; /* New list                                    */
  ZONE *za       = NULL; /* Pointers to zones                           */
  ZONE *zb       = NULL;
  ZONE *zc       = NULL;
  int  overlap   = 0;    /* Pairs overlapping zones                     */

  /* Check input                                                        */
  if(!inputtest || !inputlist)
  {
     /* No Zone or Zone List Specified                                  */
     return(-1);
  }

  /* Duplicate test zone and put at head of list of zones               */
  INIT(zonelist,ZONE);
  zonelist->chain1       = inputtest->chain1;
  zonelist->start1       = inputtest->start1;
  zonelist->startinsert1 = inputtest->startinsert1;
  zonelist->stop1        = inputtest->stop1;
  zonelist->stopinsert1  = inputtest->stopinsert1;
  zonelist->chain2       = inputtest->chain2;
  zonelist->start2       = inputtest->start2;
  zonelist->startinsert2 = inputtest->startinsert2;
  zonelist->stop2        = inputtest->stop2;
  zonelist->stopinsert2  = inputtest->stopinsert2;
  zonelist->mode         = inputtest->mode;
  
  /* Convert to sequential zone                                         */
  if(zonelist->mode == ZONE_MODE_RESNUM)
  {
     if(ConvertResidueToSequential(zonelist, strucnum))
     {
        /* Could not convert test zone                                  */
        return(-1);
     }
  }
  
  /* Set zone za and zb                                                 */
  za = zb = zonelist; 
  ALLOCNEXT(zb,ZONE); 

  /* Check zones in inputlist                                           */
  for(zc=inputlist; zc!=NULL; NEXT(zc))
  {
      zb->chain1       = zc->chain1;
      zb->start1       = zc->start1;
      zb->startinsert1 = zc->startinsert1;
      zb->stop1        = zc->stop1;
      zb->stopinsert1  = zc->stopinsert1;
      zb->chain2       = zc->chain2;
      zb->start2       = zc->start2;
      zb->startinsert2 = zc->startinsert2;
      zb->stop2        = zc->stop2;
      zb->stopinsert2  = zc->stopinsert2;
      zb->mode         = zc->mode;

      /* Convert residue zones to sequential zones                      */
      if(zb->mode == ZONE_MODE_RESNUM)
         ConvertResidueToSequential(zb, strucnum);

      /* Test for overlap between zone za and zone zb                   */
      /* Ignore zones that can't be converted                           */
      if((((za->start1 <= zb->start1) && (za->stop1 >= zb->start1)) ||
          ((zb->start1 <= za->start1) && (zb->stop1 >= za->start1)) ||
          ((za->start2 <= zb->start2) && (za->stop2 >= zb->start2)) ||
          ((zb->start2 <= za->start2) && (zb->stop2 >= za->start2))) &&
         (zb->mode != ZONE_MODE_RESNUM))
         overlap++;
  }

  /* Cleanup                                                            */
  FREELIST(zonelist,ZONE);

  return(overlap);
}


/************************************************************************/
/*>void SetZoneFromBValCol(void)
   -----------------------------
  Function sets zones using markers set in the temperature factor column 
  of a PDB file. 

  Zones are marked by stretches of positive whole numbers or integers in 
  the temperature factor column. Zeros are ignored. If both the Reference 
  structure and the mobile structure(s) are marked then zones are assigned 
  between corresponding stretches of labeled residues. If only the 
  reference structure is marked then the same residue numbers for both the
  reference structure and the mobile structure zones.

  12.05.08 Original By: CTP
  16.05.08 Corrected bug with assigning zones. Added assignment from 
           Reference Structure only if there are no markers set in the 
           mobile structure. 
  22.05.08 Removed Remove-Duplicate-Zones section - not needed.
*/
void SetZoneFromBValCol(void)
{
   int  snum         = 0;
   int  label        = 0;
   int  maxbval      = 0; 
   int  ref_res      = 0;
   int  mob_res      = 0;
   PDB  *pdbref      = NULL;
   PDB  *pdbmob      = NULL;
   PDB  *currpdbref  = NULL;
   PDB  *currpdbmob  = NULL;
   ZONE *zonelist    = NULL;
   ZONE *z           = NULL;
   ZONE *zn          = NULL;  
   BOOL converged    = TRUE;
   BOOL ref_assigned = FALSE;
   BOOL mob_assigned = FALSE;
   BOOL ref_only     = FALSE;
   BOOL all_zero     = TRUE;
   
   if(!gQuiet)
      printf("   Assigning zones from Temperature Factor column.\n");

   /* Find highest BVal + Error check                                   */
   for(pdbref = gRefPDB; pdbref!= NULL; NEXT(pdbref))
   {
      if( (int)(pdbref->bval) > maxbval)
         maxbval = (int)(pdbref->bval);
      
      /* Return if column not integers/whole numbers                    */
      if(pdbref->bval != (int)(pdbref->bval))
      {
         printf("   ERROR: Integer or whole numbers");
         printf(" needed in Temperature Factor column\n");
         return;
      }
      
      /* Return if column is less than zero                             */
      if((int)(pdbref->bval) < 0)
      {
         printf("   ERROR: Temperature Factors cannot be negative.\n");
         return;
      }
   }
   
   /* Return if no values in column                                     */
   if(!maxbval)
   {
      printf("   ERROR: No values found in Temperature Factor column\n");
      return; 
   }
   
   /* Check format of Mobile Structures                                 */
   for(snum = 0; snum < gMultiCount; snum++)
   {
      for(pdbmob = gMobPDB[snum]; pdbmob != NULL; NEXT(pdbmob))
      {
         /* Check Format                                                */
         if(pdbmob->bval != (int)(pdbmob->bval))
            ref_only = TRUE;
         
         if(pdbmob->bval != 0.0)
            all_zero = FALSE;
      }
   }

   if(all_zero)
      ref_only = TRUE;
   
   /* Clear Existing Zones                                              */
   for(snum = 0; snum < gMultiCount; snum++)
   {
      if(gZoneList[snum] != NULL)
      {
         FREELIST(gZoneList[snum],ZONE);
         gZoneList[snum] = NULL;
      }
   }
   gUserFitZone = FALSE;
   
   /* Cycle through mobile structures                                   */
   for(snum = 0; snum < gMultiCount; snum++)
   {
      /* Loop through ZoneLabels (B Value Column)                       */
      for(label = 1; label <= maxbval; label++)
      {
         /* Assign Zones                                                */
         /* Clear zonelist                                              */
         if(zonelist != NULL)
         {
            FREELIST(zonelist,ZONE);
            zonelist = NULL;
         }
         
         /* Reset Residue Count                                         */
         ref_res      = 1;
         mob_res      = 1;
         pdbref       = gRefPDB;
         pdbmob       = gMobPDB[snum];
         currpdbref   = gRefPDB;
         currpdbmob   = gMobPDB[snum];
         ref_assigned = FALSE;
         mob_assigned = FALSE;
         
         /* Loop Through Ref PDB                                        */
         for(; (pdbref!= NULL) && (pdbmob != NULL); NEXT(pdbref))
         {
            if((pdbref->resnum != currpdbref->resnum)  ||
               strcmp(pdbref->chain, currpdbref->chain)||
               strcmp(pdbref->insert,currpdbref->insert))
            {
               currpdbref =  pdbref; 
               ref_res++;
               ref_assigned = FALSE;
            }
            
            /* Find Label in Reference                                  */
            if(((int)(pdbref->bval) == label) && !ref_assigned)
            {
               if(ref_only)
               {
                  /* Assign Zone from Reference Stucture Alone          */
                  ref_assigned = TRUE;
                  
                  /* Update Zone                                        */
                  if(!zonelist)
                  {
                     INIT(zonelist,ZONE);
                     z = zonelist;
                  }
                  else 
                  {
                     z=zonelist;
                     LAST(z);
                     ALLOCNEXT(z,ZONE);
                  }
                  
                  z->chain1       = ' ';
                  z->start1       = ref_res;
                  z->startinsert1 = ' ';
                  z->stop1        = ref_res;
                  z->stopinsert1  = ' ';
                  z->chain2       = ' ';
                  z->start2       = ref_res;
                  z->startinsert2 = ' ';
                  z->stop2        = ref_res;
                  z->stopinsert2  = ' ';
                  z->mode         = ZONE_MODE_SEQUENTIAL;
               }
               else
               {
                  /* Find Label in Mobile                               */
                  for(; pdbmob != NULL && !ref_assigned; NEXT(pdbmob))
                  {
                     if((pdbmob->resnum != currpdbmob->resnum)  || 
                        strcmp(pdbmob->chain, currpdbmob->chain)||
                        strcmp(pdbmob->insert,currpdbmob->insert))
                     {
                        currpdbmob =  pdbmob; 
                        mob_res++;
                        mob_assigned = FALSE;
                     }
                     
                     if(((int)(pdbmob->bval) == label) && 
                        !ref_assigned &&
                        !mob_assigned)
                     {
                        ref_assigned = TRUE;
                        mob_assigned = TRUE;
                        
                        /* Update Zone                                  */
                        if(!zonelist)
                        {
                           INIT(zonelist,ZONE);
                           z = zonelist;
                        }
                        else 
                        {
                           z=zonelist;
                           LAST(z);
                           ALLOCNEXT(z,ZONE);
                        }
                        
                        z->chain1       = ' ';
                        z->start1       = ref_res;
                        z->startinsert1 = ' ';
                        z->stop1        = ref_res;
                        z->stopinsert1  = ' ';
                        z->chain2       = ' ';
                        z->start2       = mob_res;
                        z->startinsert2 = ' ';
                        z->stop2        = mob_res;
                        z->stopinsert2  = ' ';
                        z->mode         = ZONE_MODE_SEQUENTIAL;
                     }         
                  }
               }
            }
         }
         
         /* Merge Zones                                                 */
         if(zonelist)
         {
            do
            {
               /* Assume we have converged                              */
               converged = TRUE;
               for(z=zonelist; z!=NULL; NEXT(z))
               {
                  zn = z->next;
                  if(zn)
                  {
                     /* See if the two zones are sequential             */
                     if((zn->start1 == (z->stop1 + 1)) &&
                        (zn->start2 == (z->stop2 + 1)))
                     {
                        z->stop1 = zn->stop1;
                        z->stop2 = zn->stop2;
                        z->next = zn->next;
                        free(zn);
                        converged = FALSE;
                     }
                  }
               }
            }  while(!converged);
         }
         
         
         /* Append to global ZoneList                                   */
         if(zonelist)
         {
            if(!gZoneList[snum])
            {
               gZoneList[snum] = zonelist;
            }
            else 
            {
               z = gZoneList[snum];
               LAST(z);
               z->next = zonelist;
            }
            zonelist = NULL;
         }
         
      }  /* End of loop through labels                                  */
   }  /* End of loop through mobile structures                          */
   
   gUserFitZone = TRUE;
   return;
}


/************************************************************************/
/*>int ConvertSequentialToResidue(ZONE *input_zone, int strucnum)
   --------------------------------------------------------------
   Converts sequential-numbered zone into residue-numbered zone.

   06.06.08 Original By: CTP
*/
int ConvertSequentialToResidue(ZONE *input_zone, int strucnum)
{
   PDB  *p                = NULL;
   PDB  *p_start          = NULL;
   PDB  *q                = NULL;
   PDB  *q_start          = NULL;
   ZONE *output_zonelist  = NULL;
   ZONE *z                = NULL;
   int  residue_count = 0;
   int  prev_residue  = 0;
   char prev_insert   = ' ';
   char prev_chain    = ' ';
   int  ref_start     = 0;
   int  ref_stop      = 0;
   int  mob_start     = 0;
   int  mob_stop      = 0;
   int  i             = 0;


   /* Check input zone is sequential and not NULL                       */
   if(!input_zone || input_zone->mode == ZONE_MODE_RESNUM)
      return(0);
   
   /* Zero start and stop for output_zone                               */
   INIT(output_zonelist,ZONE);
   output_zonelist->start1 = 0;
   output_zonelist->start2 = 0;
   output_zonelist->stop1  = 0;
   output_zonelist->stop2  = 0;

   /* Find Reference Residues                                           */
   residue_count = 0;
   prev_residue  = 0;
   prev_insert   = ' ';
   prev_chain    = ' ';

   /* Find ref start residue p                                          */
   ref_start = input_zone->start1;
   ref_stop  = input_zone->stop1;

   /* Start is undefined                                                */
   if(ref_start == -999)
      ref_start = 1;
   
   for(p=gRefPDB; p!=NULL && (p_start==NULL || ref_stop == -999); NEXT(p))
   {
      /* Update residue count                                           */
      if((p->resnum    != prev_residue) || 
         (p->insert[0] != prev_insert)  ||
         (p->chain[0]  != prev_chain))
      {
         residue_count++;
      }
      
      /* Set start atom                                                 */
      if(!p_start && residue_count == ref_start)
         p_start = p;
      
      /* Update previous residue                                        */
      prev_residue = p->resnum;
      prev_insert  = p->insert[0];
      prev_chain   = p->chain[0];
   }
   
   /* Set stop if undefined                                             */
   if(ref_stop == -999)
      ref_stop = residue_count;

   /* Find Mobile Residues                                              */
   residue_count = 0;
   prev_residue  = 0;
   prev_insert   = ' ';
   prev_chain    = ' ';

   /* Find mob start residue p                                          */
   mob_start = input_zone->start2;
   mob_stop  = input_zone->stop2;
   
   /* Start is undefined                                                */
   if(mob_start == -999)
      mob_start = 1;

   for(p=gMobPDB[strucnum]; 
       p!=NULL && (q_start==NULL || mob_stop == -999); 
       NEXT(p))
   {
      /* Update residue count                                           */
      if((p->resnum    != prev_residue) ||
         (p->insert[0] != prev_insert)  ||
         (p->chain[0]  != prev_chain))
      {
         residue_count++;
      }

      /* Set start atom                                                 */
      if(!q_start && residue_count == mob_start)
         q_start = p;
      
      /* Update previous residue                                        */
      prev_residue = p->resnum;
      prev_insert  = p->insert[0];
      prev_chain   = p->chain[0];
   }

   /* Set stop if undefined                                             */
   if(mob_stop == -999)
      mob_stop = residue_count;

   /* Error Checking                                                    */
   
   /* Error: Ref Start Not Set                                          */
   if(!p_start)
   {
      printf("   Error: Reference start residue not found.\n");
      return(1);
   }
   
   /* Error: Mob Start Not Set                                          */
   if(!q_start)
   {
      printf("   Error: Mobile start residue not found.\n");
      return(1);
   }
   
   /* Error: Number of Residues Not Matched                             */
   if((ref_stop - ref_start) != (mob_stop - mob_start))
   {
      char zone1[64];
      char zone2[64];
      
      printf("   Error: Number of residues in zone does not match.\n");
      
      FormatZone(zone1, input_zone->chain1, 
                 input_zone->start1, input_zone->startinsert1, 
                 input_zone->stop1,  input_zone->stopinsert1);
         
      FormatZone(zone2, input_zone->chain2, 
                 input_zone->start2, input_zone->startinsert2, 
                 input_zone->stop2,  input_zone->stopinsert2);
    
      printf("           %-16s with %-16s %s\n",
             zone1, zone2,
             ((input_zone->mode == ZONE_MODE_RESNUM)?"(Residue numbering)"
              :"(Sequential numbering)"));

      printf("            Reference: %d, Mobile: %d\n\n", 
             (ref_stop - ref_start + 1),(mob_stop - mob_start + 1));
         
      return(1);
   }

   /* Assign the Zones                                                  */
   
   /* Set Start Residues                                                */
   p = p_start;
   q = q_start;
   
   /* Set Zone start1 resnum,insert,chain & start2 resnum,insert,chain  */
   z = output_zonelist;
   z->chain1       = p->chain[0];
   z->start1       = p->resnum;
   z->startinsert1 = p->insert[0];
   z->stop1        = p->resnum;
   z->stopinsert1  = p->insert[0];
   z->chain2       = q->chain[0];
   z->start2       = q->resnum;
   z->startinsert2 = p->insert[0];
   z->stop2        = q->resnum;
   z->stopinsert2  = q->insert[0];
   z->mode         = ZONE_MODE_RESNUM;
   z->next         = NULL;
   
   /* Find continuous zones                                             */
   for(i=ref_start; i < ref_stop; i++)
   {
      /* Skip to Next Ref                                               */
      while(p->chain[0]  == z->chain1 && 
            p->resnum    == z->stop1  &&
            p->insert[0] == z->stopinsert1)
      {
         NEXT(p);
      }
      
      /* Skip to Next Mob                                               */
      while(q->chain[0]  == z->chain2 && 
            q->resnum    == z->stop2  &&
            q->insert[0] == z->stopinsert2)
      {
         NEXT(q);
      }
      
      /* Check Prev Chain                                               */
      if(p->chain[0] != z->chain1 || q->chain[0] != z->chain2)
      {
         /* Add New Zone                                                */
         ALLOCNEXT(z,ZONE);
         
         /* Set Start + Stop Position                                   */
         z->chain1       = p->chain[0];
         z->start1       = p->resnum;
         z->startinsert1 = p->insert[0];
         z->stop1        = p->resnum;
         z->stopinsert1  = p->insert[0];
         z->chain2       = q->chain[0];
         z->start2       = q->resnum;
         z->startinsert2 = p->insert[0];
         z->stop2        = q->resnum;
         z->stopinsert2  = q->insert[0];
         z->mode         = ZONE_MODE_RESNUM;
         z->next         = NULL;
      }
      else 
      {
         /* Update Stop Position                                        */
         z->stop1        = p->resnum;
         z->stopinsert1  = p->insert[0];
         z->stop2        = q->resnum;
         z->stopinsert2  = q->insert[0];
      }
   }
   
   
   /* Splice output_zonelist into input_zone position                   */
   /* Set End of List                                                   */
   z->next = input_zone->next;
   
   /* Set Start of List                                                 */
   input_zone->chain1       = output_zonelist->chain1;
   input_zone->start1       = output_zonelist->start1;
   input_zone->startinsert1 = output_zonelist->startinsert1;
   input_zone->stop1        = output_zonelist->stop1;
   input_zone->stopinsert1  = output_zonelist->stopinsert1;
   input_zone->chain2       = output_zonelist->chain2;
   input_zone->start2       = output_zonelist->start2;
   input_zone->startinsert2 = output_zonelist->startinsert2;
   input_zone->stop2        = output_zonelist->stop2;
   input_zone->stopinsert2  = output_zonelist->stopinsert2;
   input_zone->mode         = output_zonelist->mode;
   input_zone->next         = output_zonelist->next;
   
   /* Cleanup - this is only used as a single structure and not a linked
      list so it only needs a free()
   */
   free(output_zonelist);
   
   return(0);
}


/************************************************************************/
/*>ZONE *SortZoneList(ZONE *zonelist)
   ----------------------------------
   Sorts a sequentially-numbered zone list by setting head of sorted list 
   and inserting zones into sorted list either when start of next zone is 
   greater than current zone or at end of sorted list. Residue-numbered 
   zones are appended to the end of the linked list.

   10.06.08 Original  By: CTP
*/
ZONE *SortZoneList(ZONE *zonelist)
{
   ZONE *sortlist = NULL;
   ZONE *z_curr   = NULL;
   ZONE *z_next   = NULL;
   ZONE *z_sort   = NULL;

   /* Make sortlist with blank zone at head as placeholder              */
   INIT(sortlist,ZONE);
   
   z_curr = zonelist; /* Set current zone at head of zonelist           */
   z_sort = sortlist; /* Set sort zone at head of sortlist              */
   
   
   /* Insertion Sort (sort of...)                                       */
   z_curr = zonelist; /* Set current zone at head of zonelist           */
   while(z_curr != NULL)
   {
      /* Set zone pointers                                              */
      z_sort = sortlist;
      z_next = z_curr->next;
      
      /* Append unconverted zones to end of sortlist                    */
      if(z_curr->mode == ZONE_MODE_RESNUM)
      {
         LAST(z_sort);
         z_sort->next = z_curr;
         z_curr->next = NULL;
         z_sort       = NULL;
      }
      
      /* Scan through sortlist                                          */
      while(z_sort != NULL)
      { 
         if(z_sort->next != NULL)
         {
            /* Insert curr_zone in sortlist                             */
            if((z_sort->next->start1 > z_curr->start1) ||
               (z_sort->next->mode == ZONE_MODE_RESNUM))
            {
               z_curr->next = z_sort->next;
               z_sort->next = z_curr;
               z_sort       = NULL;
            }
            else
            {
               NEXT(z_sort);
            }
         }
         else 
         {
            /* Append unmatched zone to end of sortlist                 */
            z_sort->next = z_curr;
            z_curr->next = NULL;
            z_sort       = NULL;
         }
      }
      
      /* Set current zone to next zone                                  */
      z_curr = z_next;
   }
   
   /* Return sorted list                                                */
   zonelist = sortlist->next;
   sortlist->next = NULL;
   FREELIST(sortlist,ZONE);
   
   return(zonelist);
}


/************************************************************************/
/*>int ConvertZoneList(ZONE *zonelist, int strucnum, int mode)
   -----------------------------------------------------------
   Convert zonelist to numbering mode.

   10.06.08 Original By: CTP
*/
int ConvertZoneList(ZONE *zonelist, int strucnum, int mode)
{
   ZONE *z   = NULL;
   int error = 0;
   
   for(z=zonelist; z!=NULL; NEXT(z))
   {
      if(z->mode == mode)
      {
         continue;
      }
      
      if(mode == ZONE_MODE_RESNUM)
         error += ConvertSequentialToResidue(z, strucnum);
      else
         error += ConvertResidueToSequential(z, strucnum);
   }
   
   return(error);
}


/************************************************************************/
/*>int ConvertAllZones(int mode)
   -----------------------------
   Convert all zonelists to numbering mode.

   10.06.08 Original By: CTP
*/
int ConvertAllZones(int mode)
{
   int strucnum = 0;
   int error    = 0;
   
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      /* Convert Fitting Zones                                          */
      error += ConvertZoneList(gZoneList[strucnum], strucnum, mode);
      
      /* Convert RMSd Calculation Zones                                 */
      error += ConvertZoneList(gRZoneList[strucnum], strucnum, mode);
      
      /* Convert Centre Zone                                            */
      error += ConvertZoneList(gCZoneList[strucnum], strucnum, mode);
   }
   return(error);
}


/************************************************************************/
/*> int SortAllZones(void)
    ----------------------
    Sort all zone lists.

    10.06.08 Original By: CTP
*/
int SortAllZones(void)
{
   int strucnum = 0;
   int error    = 0;
   
   for(strucnum=0; strucnum<gMultiCount; strucnum++)
   {
      /* Convert Zones to Sequential Numbering                          */
      error += ConvertZoneList(gZoneList[strucnum],  strucnum, 
                               ZONE_MODE_SEQUENTIAL);
      error += ConvertZoneList(gRZoneList[strucnum], strucnum, 
                               ZONE_MODE_SEQUENTIAL);
      
      /* Sort Zonelists                                                 */
      gZoneList[strucnum]  = SortZoneList(gZoneList[strucnum]);
      gRZoneList[strucnum] = SortZoneList(gRZoneList[strucnum]);
   }
   return(error);
}


/************************************************************************/
/*>ZONE *ChainList(PDB *pdblist)
   -----------------------------
   Takes a PDB sequence and returns the chains as a linked list of ZONEs. 
   Output is a sequentially numbered zonelist with the chain ID set.

   Unlike the ZONEs defined by other functions (such as SetFitZone()) the
   ZONEs defined by this function only cover a single chain rather than 
   indicating the equivalent regions of two chains. 

   25.06.08 Original  By: CTP
*/
ZONE *ChainList(PDB *pdblist)
{
   ZONE *chainlist, *z; 
   PDB  *p            = NULL;
   int  residue_count = 0;
   int  prev_residue  = 0;
   char prev_insert   = ' ';
   char prev_chain    = ' ';
   
   /* Return if no pdblist                                              */
   if(!pdblist)
   {
      return(NULL);
   }
   
   /* Zero start stop for chainlist                                     */
   INIT(chainlist,ZONE);
   chainlist->start1 = 1;
   chainlist->mode   = ZONE_MODE_SEQUENTIAL;
   
   /* Zero Residue Count                                                */
   residue_count = 0;
   prev_residue  = 0;
   prev_insert   = ' ';
   
   /* Set Chain ID                                                      */
   chainlist->chain1 = pdblist->chain[0];
   prev_chain        = pdblist->chain[0];
   
   /* Cycle through pdblist                                             */
   for(z = chainlist, p=pdblist; p!=NULL; NEXT(p))
   {
      /* Update residue count                                           */
      if((p->resnum    != prev_residue) ||
         (p->insert[0] != prev_insert)  ||
         (p->chain[0]  != prev_chain))
      {
         residue_count++;
      }
      
      /* If chain changes then add zone                                 */
      if(p->chain[0]  != prev_chain)
      {
         z->stop1 = residue_count - 1;
         ALLOCNEXT(z,ZONE);
         z->start1 = residue_count;
         z->chain1 = p->chain[0];
         z->mode   = ZONE_MODE_SEQUENTIAL;
      }
      
      /* Update previous residue                                        */
      prev_residue = p->resnum;
      prev_insert  = p->insert[0];
      prev_chain   = p->chain[0];
      
   }
   
   /* Set Final Stop                                                    */
   z->stop1 =  residue_count;
   
   /* Return                                                            */
   return(chainlist);
}


/************************************************************************/
/*>BOOL SequentialZones(int strucnum)
   ----------------------------------
   Checks zones in a list occur sequentially along the reference and 
   mobile protein sequences for each chain.

   Fit zones CANNOT be converted into sequence alignments if the zones 
   aren't in sequence.

   Note: The zones in the zonelist should be converted to SEQUENTIAL 
   numbering mode and sorted prior to calling this function.

   31.07.08 Original  By: CTP
*/
BOOL SequentialZones(int strucnum)
{
   ZONE *chainlist_mob, *chainlist_ref;
   ZONE *ref, *mob, *z, *p;
   
   /* Set chains                                                        */
   chainlist_ref = ChainList(gRefPDB);
   chainlist_mob = ChainList(gMobPDB[strucnum]);
 
   /* Cycle through ref                                                 */
   for(ref=chainlist_ref; ref!=NULL; NEXT(ref))
   {
      /* Cycle through mob                                              */
      for(mob=chainlist_mob; mob!=NULL; NEXT(mob))
      {
         p=NULL; /* Reset previous residue pointer                      */
         for(z=gZoneList[strucnum];z!=NULL; NEXT(z))
         {
            /* Check zone between ref and mob                           */
            if((z->start1 >= ref->start1)&&(z->stop1 <= ref->stop1)&&
               (z->start2 >= mob->start1)&&(z->stop2 <= mob->stop1)&&
               (z->mode ==  ZONE_MODE_SEQUENTIAL))
            {
               if(p==NULL)
               {
                  p=z;
                  continue;
               }
               
               if((p->stop1  >= z->start1)||
                  (p->stop2  >= z->start2))
               {
                  /* Return - Out of Sequence                           */
                  return(FALSE);
               }
               else
               {
                  p=z;
               }
            }
         }
      }
   }
   
   /* Return - In Sequence                                              */
   return(TRUE);
}


/************************************************************************/
/*>BOOL SequentialZonesWholeSeq(int strucnum)
   ------------------------------------------
   Checks zones in a list occur sequentially along the reference and 
   mobile protein sequences for whole sequence.

   Fit zones CANNOT be converted into sequence alignments if the zones 
   aren't in sequence.

   Note: The zones in the zonelist should be converted to SEQUENTIAL 
   numbering mode and sorted prior to calling this function.

   30.01.09 Original based on SequentialZones()  By: CTP
*/
BOOL SequentialZonesWholeSeq(int strucnum)
{
   ZONE *z = NULL;
   ZONE *p = NULL;
   
   
   for(z=gZoneList[strucnum];z!=NULL; NEXT(z))
   {
      /* Check zone numbering                                           */
      if(z->mode !=  ZONE_MODE_SEQUENTIAL)
         return(FALSE);
      
      /* Set p = z for first zone                                       */
      if(p==NULL)
      {
         p=z;
         continue;
      }
      
      if((p->stop1  >= z->start1)|| (p->stop2  >= z->start2))
      {
         /* Return - Out of Sequence                                    */
         return(FALSE);
      }
      else
      {
         p=z;
      }
   }
   
   /* Return - In Sequence                                              */
   return(TRUE);
}


/************************************************************************/
/*>BOOL OneToOneChains(int strucnum)
   ---------------------------------
   Function tests to see if there is a one to one match between aligned
   chains when assigning zones based on alignment. 

   20.10.08 Original  By: CTP
*/
BOOL OneToOneChains(int strucnum)
{
   ZONE *chainlist[2];
   ZONE *ref, *mob, *z;
   BOOL found    = FALSE;
   BOOL OneToOne = TRUE;   /* Return Value                              */
   char chainid  = ' ';
   int i;
   
   if(gZoneList[strucnum] == NULL) return(TRUE);
   
   chainlist[0] = ChainList(gRefPDB);
   chainlist[1] = ChainList(gMobPDB[strucnum]);
   
   /* Check each ref chain goes to one mob chain                        */
   /* Check each mob chain goes to one ref chain                        */
   
   for(i=0; i<2; i++)
   {
      for(ref = chainlist[i]; ref!=NULL; NEXT(ref))
      {
         /* Note:
            The expression [i * -1 + 1] alternates between 0 and 1
            when: i = 0, i * -1 +1 = 1 
                  i = 1, i * -1 +1 = 0  
         */
         for(mob = chainlist[i * -1 + 1]; mob!=NULL; NEXT(mob))
         {
            for(z=gZoneList[strucnum]; z!=NULL; NEXT(z))
            {
               /* Find zones in both chains                             */
               if(((i == 0) &&
                   (z->start1>=ref->start1)&&(z->stop1<=ref->stop1)&&
                   (z->start2>=mob->start1)&&(z->stop2<=mob->stop1)&&
                   (z->mode ==  ZONE_MODE_SEQUENTIAL)) ||
                  ((i == 1) &&
                   (z->start1>=mob->start1)&&(z->stop1<=mob->stop1)&&
                   (z->start2>=ref->start1)&&(z->stop2<=ref->stop1)&&
                   (z->mode ==  ZONE_MODE_SEQUENTIAL)))
               {
                  if(!found)
                  {
                     /* Remember first chain ID                         */
                     chainid = mob->chain1;
                     found = TRUE;
                  }
                  else
                  {
                     /* Return FALSE if different chain ID found        */
                     if(mob->chain1 != chainid)
                     {
                        OneToOne = FALSE;
                        goto Return;
                     }
                  }
                  
               }
            }
         }

         /* Reset found flag                                            */
         found = FALSE;
      }
   }

Return:
   /* Free memory and return                                            */
   if(chainlist[0]) FREELIST(chainlist[0],ZONE);
   if(chainlist[1]) FREELIST(chainlist[1],ZONE);
   
   return(OneToOne);
}


/************************************************************************/
/*>int EnforceOneToOneChains(int strucnum)
   ---------------------------------------
   Enforces a one-to-one alignment of chains between the reference 
   structure and the mobile structure. 

   31.07.08 Original   By: CTP
*/
int EnforceOneToOneChains(int strucnum)
{
   ZONE *chainlist_mob, *chainlist_ref;
   ZONE *zonelist, *deletelist; 
   ZONE *ref, *mob, *z, *p, *d;
   
   int best;
   int curr;
   
   /* Null Pointers                                                     */
   chainlist_mob = chainlist_ref = zonelist = deletelist = NULL;
   ref = mob = z = p = d = NULL;
   
   /* Set blank place holder                                            */
   INIT(zonelist,  ZONE);
   INIT(deletelist,ZONE);
   
   /* Copy global zones to function                                     */
   zonelist->next = gZoneList[strucnum];
   
   /* Set chains                                                        */
   chainlist_ref = ChainList(gRefPDB);
   chainlist_mob = ChainList(gMobPDB[strucnum]);
   
   /*** Pass A: Ref to Mob                                            ***/
   
   /* Cycle through ref                                                 */
   for(ref=chainlist_ref; ref!=NULL; NEXT(ref))
   {
      best = 0;
      ref->chain2 = chainlist_mob->chain1;
      ref->start2 = chainlist_mob->start1;
      ref->stop2  = chainlist_mob->stop1;
      
      /* Cycle through mob                                              */
      for(mob=chainlist_mob; mob!=NULL; NEXT(mob))
      {
         curr = 0;
         /* Count matching residues                                     */
         for(z=zonelist->next;z!=NULL; NEXT(z))
         {
            if((z->start1 >= ref->start1)&&(z->stop1 <= ref->stop1)&&
               (z->start2 >= mob->start1)&&(z->stop2 <= mob->stop1)&&
               (z->mode ==  ZONE_MODE_SEQUENTIAL))
            {
               curr += (z->stop1 - z->start1 + 1);
            }
         }
         
         /* Set current best match.                                     */
         if(curr > best)
         {
            best       = curr;
            ref->chain2 = mob->chain1;
            ref->start2 = mob->start1;
            ref->stop2  = mob->stop1;
         }
      }
   }
   
   /*** Pass B: Mob to Ref                                            ***/
   
   /* Cycle through mob                                                 */
   for(mob=chainlist_mob; mob!=NULL; NEXT(mob))
   {
      best = 0;
      mob->chain2 = chainlist_ref->chain1;
      mob->start2 = chainlist_ref->start1;
      mob->stop2  = chainlist_ref->stop1;
      

      /* Cycle through ref                                              */
      for(ref=chainlist_ref; ref!=NULL; NEXT(ref))
      {
         curr = 0;
         /* Count matching residues                                     */
         for(z=zonelist->next;z!=NULL; NEXT(z))
         {
            if((z->start1 >= ref->start1)&&(z->stop1 <= ref->stop1)&&
               (z->start2 >= mob->start1)&&(z->stop2 <= mob->stop1)&&
               (z->mode ==  ZONE_MODE_SEQUENTIAL))
            {
               curr += (z->stop1 - z->start1 + 1);
            }
         }
      
         /* Set current best match.                                     */
         if(curr > best)
         {
            best       = curr;
            mob->chain2 = ref->chain1;
            mob->start2 = ref->start1;
            mob->stop2  = ref->stop1;
         }
      }
   }

   /* Keep / Delete Structures                                          */
   /* ------------------------                                          */
   
   /* Set delete pointer                                                */
   d = deletelist;
   
   /*** Pass A Ref to Mobile                                          ***/
   /* Cycle through ref                                                 */
   for(ref=chainlist_ref; ref!=NULL; NEXT(ref))
   {
      /* Cycle through mob                                              */
      for(mob=chainlist_mob; mob!=NULL; NEXT(mob))
      {
         /* Cycle through zones                                         */
         p = zonelist; /* previous pointer                              */
         for(z=zonelist->next;z!=NULL; NEXT(z))
         {
            if(((z->start1 >= ref->start1)&&(z->stop1 <= ref->stop1))&&
               ((z->start2 <  ref->start2)||(z->stop2 >  ref->stop2))&&
               (z->mode ==  ZONE_MODE_SEQUENTIAL))
            {
               /* Remove zone from zonelist                             */
               p->next = z->next;
               
               /* Append zone to deletelist                             */
               d->next = z;
               z->next = NULL;
               d       = z;
               
               /* Set zone pointer to previous zone                     */
               z       = p;
            }
            else 
            {
               p = z;
            }
         }
      }
   }
   
   /*** Pass B Mobile to Ref                                          ***/
   /* Cycle through mob                                                 */
   for(mob=chainlist_mob; mob!=NULL; NEXT(mob))
   {
      /* Cycle through mob                                              */
      for(ref=chainlist_ref; ref!=NULL; NEXT(ref))
      {
         /* Cycle through zones                                         */
         for(z=zonelist->next;z!=NULL; NEXT(z))
         {
            if(((z->start2 >= mob->start1)&&(z->stop2 <= mob->stop1))&&
               ((z->start1 <  mob->start2)||(z->stop1 >  mob->stop2))&&
               (z->mode ==  ZONE_MODE_SEQUENTIAL))
            {
               /* Remove zone from zonelist                             */
               p->next = z->next;
               
               /* Append zone to deletelist                             */
               d->next = z;
               z->next = NULL;
               d       = z;
               
               /* Set zone pointer to previous zone                     */
               z       = p;
            }
            else 
            {
               p = z;
            }
         }
      }
   }
   
   /* Announce deleted zones                                            */
   if(!gQuiet)
   {
      if(deletelist->next != NULL)
      {
         /* Convert to Current Numbering Scheme                         */
         ConvertZoneList(deletelist, strucnum, gCurrentMode);
         
         printf("   Deleted zones:\n");
         for(d=deletelist->next; d!=NULL; NEXT(d))
         {
            char zone1[64], zone2[64];
            
            FormatZone(zone1, d->chain1, 
                       d->start1, d->startinsert1, 
                       d->stop1,  d->stopinsert1);
            
            FormatZone(zone2, d->chain2, 
                       d->start2, d->startinsert2, 
                       d->stop2,  d->stopinsert2);
            
            printf("      %-16s with %-16s %s\n",
                   zone1, zone2,
                   ((d->mode == ZONE_MODE_RESNUM)
                    ?"(Residue numbering)"
                    :"(Sequential numbering)"));
         }
      }
      else
      {
         printf("   No zones deleted\n");
      }
   }
   
   /* Reset Global Zonelist                                             */
   gZoneList[strucnum] = zonelist->next;
   zonelist->next      = NULL;
   
   /* Free Memory                                                       */
   FREELIST(zonelist  ,ZONE);
   FREELIST(deletelist,ZONE);
   
   return(0);
}


/************************************************************************/
/*>int CopyPDBListToRef(int strucnum)
   ----------------------------------
   Copies gMobPDB[strucnum] to gRefPDB. Called by SetMobileToReference().

   20.10.08 Original  By: CTP
   05.11.08 Fixed bug - Function scanned through whole PDBList when 
            allocating memory for next item in new reference list. 
*/
int CopyPDBListToRef(int strucnum)
{
   PDB    *mp = NULL;
   PDB    *rp = NULL;
   int natoms = 0;
   
   for(mp = gMobPDB[strucnum]; mp != NULL; NEXT(mp), natoms++)
   {
      if(!gRefPDB)
      {
         INIT(gRefPDB,PDB);
         rp = gRefPDB;
      }
      else 
      {
         ALLOCNEXT(rp,PDB);
      }

      if(rp)
      {
         rp->atnum  = mp->atnum;
         rp->resnum = mp->resnum;
         rp->x      = mp->x;
         rp->y      = mp->y;
         rp->z      = mp->z;
         rp->occ    = mp->occ;
         rp->bval   = mp->bval;
         rp->altpos = mp->altpos;
         
         strcpy(rp->record_type, mp->record_type);
         strcpy(rp->atnam,       mp->atnam);
         strcpy(rp->atnam_raw,   mp->atnam_raw);
         strcpy(rp->resnam,      mp->resnam);
         strcpy(rp->chain,       mp->chain);
         strcpy(rp->insert,      mp->insert);
         
         rp->next = NULL;
      }
      else
      {
         printf("   Error==> No memory for operation!\n");
         return(0);
      }
   }
  
   return(natoms);
}


/************************************************************************/
/*>int SetMobileToReference(int strucnum)
   --------------------------------------
   This function sets a mobile structure to the reference structure. This
   function is called by the SETREF command and ALLVsAllRMS(). 

   20.10.08 Original By: CTP
*/
int SetMobileToReference(int strucnum)
{
   int  natoms = 0;
   int  i = 0;
   
   /* Copy Filename                                                     */
   strcpy(gRefFilename,gMobFilename[strucnum]);
   
   /* Free current reference WHOLEPDB                                   */
   if(gRefWPDB)
   {
      STRINGLIST *s = NULL;
      
      if(gRefWPDB->header)  FreeStringList(gRefWPDB->header);
      if(gRefWPDB->trailer) FreeStringList(gRefWPDB->trailer);
      gRefWPDB->header  = NULL;
      gRefWPDB->trailer = NULL;
      
      for(s=gMobWPDB[strucnum]->header; s!=NULL; NEXT(s))
      {
         gRefWPDB->header = StoreString(gRefWPDB->header, s->string);
      }
      
      for(s=gMobWPDB[strucnum]->trailer; s!=NULL; NEXT(s))
      {
         gMobWPDB[strucnum]->trailer = StoreString(gRefWPDB->trailer, 
                                                   s->string);
      }
      
      gRefWPDB->natoms = 0;
   }
   
   /* Free current reference PDB                                        */
   if(gRefPDB)
   {
      FREELIST(gRefPDB,PDB);
      gRefPDB  = NULL;
   }
   
   /* Copy Mobile PDB List to Reference                                 */
   natoms = CopyPDBListToRef(strucnum);
   
   if(gRefWPDB)
   {
      gRefWPDB->natoms = natoms;
   }
   
   /* Allocate coordinate array                                         */
   if(gRefCoor) free(gRefCoor);
   if((gRefCoor = (COOR *)malloc(natoms * sizeof(COOR))) == NULL)
      printf("   Error==> Unable to allocate reference coordinate \
memory!\n");
   
   /* Convert sequence                                                  */
   if(gRefSeq != NULL) free(gRefSeq);
   if((gRefSeq = PDB2Seq(gRefPDB))==NULL)
      printf("   Error==> Unable to read sequence for reference \
structure!\n");
   
   /* Renumber old ref zones with new mobile zones.                     */
   for(i=0;i<gMultiCount && gUserFitZone;i++)
   {
      ZONE *za, *zb;
      za=gZoneList[i];
      zb=gZoneList[strucnum];
      
      for(;za != NULL && zb != NULL; NEXT(za), NEXT(zb))
      {         
         za->chain1       = zb->chain2;
         za->start1       = zb->start2;
         za->startinsert1 = zb->startinsert2;
         za->stop1        = zb->stop2;
         za->stopinsert1  = zb->stopinsert2;
      }
      
      za=gRZoneList[i];
      zb=gRZoneList[strucnum];
      
      for(;za != NULL && zb != NULL; NEXT(za), NEXT(zb))
      {         
         za->chain1       = zb->chain2;
         za->start1       = zb->start2;
         za->startinsert1 = zb->startinsert2;
         za->stop1        = zb->stop2;
         za->stopinsert1  = zb->stopinsert2;
      }
      
      za=gCZoneList[i];
      zb=gRZoneList[strucnum];
      
      for(;za != NULL && zb != NULL; NEXT(za), NEXT(zb))
      {         
         za->chain1       = zb->chain2;
         za->start1       = zb->start2;
         za->startinsert1 = zb->startinsert2;
         za->stop1        = zb->stop2;
         za->stopinsert1  = zb->stopinsert2;
      }
   }
   
   return(0);
}


/************************************************************************/
/*>int fit_order_cmp(const void *ptr_scoreA, const void *ptr_scoreB)
   -----------------------------------------------------------------
   Comparison function for quicksorting a 2D array in ascending order. 
   The first element of the array is the structure number and the second 
   element is a score (usually RMSD). The array elements are first 
   compared by the score and (if necessary) the structure number. The 
   comparison returns either -1 or 1 giving a stable sort.

   20.10.08 Original By: CTP
*/
int fit_order_cmp(const void *ptr_scoreA, const void *ptr_scoreB)
{
   const REAL *pa = *(const REAL (*)[2])ptr_scoreA;
   const REAL *pb = *(const REAL (*)[2])ptr_scoreB;

   if(pa[1] != pb[1])
   {
      /* Sort by score (ascending)                                      */
      return((pa[1] < pb[1])?-1:1);
   }
   else 
   {
      /* Sort by structure number (ascending)                           */
      return((pa[0] < pb[0])?-1:1);      
   }
}


/************************************************************************/
/*>int AllVsAllRMS(char *filename, BOOL print_tab, BOOL set_ref)
   -------------------------------------------------------------
   Runs an all vs all comparison of mobile structures. The function 
   prints the all vs all comparison as a tab-delimited list if the 
   print_tab flag is set. The function prints to stdout by default but 
   will print to a file if a filename is supplied.

   The set_ref flag sets function to select the most central mobile 
   structure as the reference structure by performing an all vs all 
   comparison and selecting the structure with the lowest overall RMSD to
   the other structures.

   20.10.08 Original By: CTP
   31.10.08 Set to detect error when printing tab-delim output.
   07.11.08 Sets gMultiRef for automatic selection. Resets to gMultiRef
            after all vs all.
   25.11.08 Added error message under AllVsAll output.
*/
int AllVsAllRMS(char *filename, BOOL print_tab, BOOL set_ref)
{
   FILE *fp = stdout;
   int mob_i, mob_j; 
   BOOL error_found = FALSE;
   
   REAL (*sortrms)[2] = NULL;
   REAL **all_rms = NULL;
   
   /* Set gMultiVsRef to TRUE                                           */
   BOOL multivsref = gMultiVsRef;
   gMultiVsRef = TRUE;
   
   /* Open Output File/Pipe                                             */
   if(filename)
   {
      if((fp=OpenOrPipe(filename))==NULL)
      {
         printf("   Warning==> Uunable to open file: %s\n", filename);
         fp = stdout;
      }
   }
   
   /* Allocate memory                                                   */
   sortrms = malloc(gMultiCount * sizeof *sortrms);
   all_rms = malloc(gMultiCount * sizeof(REAL *));
   for(mob_i = 0; mob_i < gMultiCount; mob_i++)
   {
      all_rms[mob_i] = malloc(gMultiCount * sizeof(REAL));
   }
   
   /* Feedback for user                                                 */
   if(!gQuiet && !set_ref)
      printf("   All vs All...\n");
   
   /* Trim Zones                                                        */
   TrimZones();
   
   /* All vs All RMSD                                                   */
   for(mob_i = 0;mob_i < gMultiCount; mob_i++)
   {
      for(mob_j = 0;mob_j < gMultiCount; mob_j++)
         all_rms[mob_i][mob_j] = 0.0;
   }
   
   for(mob_i = 0;mob_i < gMultiCount - 1; mob_i++)
   {
      SetMobileToReference(mob_i);
      
      /* all_rms[mob_i][mob_i] = 0.0; */
      
      for(mob_j = mob_i + 1 ;mob_j < gMultiCount; mob_j++)
      {
         REAL rmsd = 0.0;
         printf("   RMS %2d vs %2d ",mob_i+1,mob_j+1);
         rmsd = FitSingleStructure(mob_j, TRUE);
         all_rms[mob_i][mob_j] = rmsd;
         all_rms[mob_j][mob_i] = rmsd;

         /* Error Flag                                                  */
         if(rmsd == -1.0) error_found = TRUE;
      }
   }
   
   /* Print All vs All Comparison                                       */
   if(print_tab)
   {
      for(mob_i = 0; mob_i < gMultiCount; mob_i++)
      {
         fprintf(fp,"\t%d",mob_i+1);
      }  
      fprintf(fp,"\n");
      
      for(mob_i = 0; mob_i < gMultiCount; mob_i++)
      {
         fprintf(fp,"%d",mob_i+1);
         for(mob_j = 0; mob_j < gMultiCount; mob_j++)
         {
            if(all_rms[mob_i][mob_j] != -1.0)
               fprintf(fp,"\t%.3f",all_rms[mob_i][mob_j]);
            else 
               fprintf(fp,"\terror");
         }
         fprintf(fp,"\n");
      }
      
      /* Print Error Message                                            */
      if(error_found)
      {
         printf("\n   Error: Error found calculating RMSD");
         printf(" in all vs all comparison\n");
      }
   }
   
   /* Set Reference Structure                                           */
   if(set_ref)
   {
      /* Auto Set Reference Structure                                   */
      /* Find sum of RMSDs from all from all other structures           */
      for(mob_i = 0; mob_i < gMultiCount; mob_i++)
      {
         sortrms[mob_i][0] = (REAL)mob_i;
         sortrms[mob_i][1] = 0.0;
         
         for(mob_j = 0; mob_j < gMultiCount; mob_j++)
         {
            sortrms[mob_i][1] += all_rms[mob_i][mob_j];
         }
      }
      
      /* Sort Array to Set Fit Order                                    */
      qsort(sortrms,gMultiCount,sizeof(*sortrms),fit_order_cmp);
      
      /* Set Reference Structure                                        */
      SetMobileToReference((int)sortrms[0][0]);
      gMultiRef = (int)sortrms[0][0];
      
      if(!gQuiet)
      {
         printf("   Mobile structure %d used as reference.\n",
                (int)sortrms[0][0] + 1);
      }
   }
   else 
   {
      /* Reset reference to current multistructure reference            */
      SetMobileToReference(gMultiRef);
   }
   
   /* Reset gMultiVsRef                                                 */
   gMultiVsRef = multivsref;
   
   /* Close Output File/Pipe                                            */
   if(fp != stdout)
      CloseOrPipe(fp);
   
   /* Free Memory                                                       */
   free(sortrms);
   for(mob_i = 0; mob_i < gMultiCount; mob_i++)
   {
      free(all_rms[mob_i]);
   }
   free(all_rms);
   
   return(0);
}


/************************************************************************/
/*>ZONE *SetOverlappingZones(ZONE *InputA, ZONE *InputB)
   -----------------------------------------------------
   This function is called by TrimZones() and returns a zonelist with the
   sections of sequence in the reference structure that are common to both
   the zonelists InputA and InputB.

   20.10.08 Original By: CTP
   04.02.09 Altered if statement from form 'x == y == z' to form 
            'x == z && y == z' as first form generates warning on Mac.
 */
ZONE *SetOverlappingZones(ZONE *InputA, ZONE *InputB)
{
   ZONE *OutputList;
   ZONE *za, *zb, *zo;
   
   OutputList = za = zb = zo = NULL;
   
   
   for(za = InputA; za != NULL; NEXT(za))
   {
      for(zb = InputB; zb != NULL; NEXT(zb))
      {
         /* Test for overlap between zone za and zone zb                */
         if(((za->start1 <= zb->start1 && za->stop1 >= zb->start1)  ||
             (zb->start1 <= za->start1 && zb->stop1 >= za->start1)) &&
            (za->mode == ZONE_MODE_SEQUENTIAL && 
             zb->mode == ZONE_MODE_SEQUENTIAL))
         {
            /* Add to Output List                                       */
            if(!OutputList)
            {
               /* Init zonelist                                         */
               INIT(OutputList,ZONE);
               zo =  OutputList;
            }
            else 
            {
               /* Append zone                                           */
               zo =  OutputList;
               LAST(zo);
               ALLOCNEXT(zo,ZONE);
            }
            
            /* Set Start/Stop                                           */
            zo->start1 = (za->start1 <= zb->start1) ? 
                         zb->start1 : za->start1;
            zo->stop1  = (za->stop1  <= zb->stop1 ) ? 
                         za->stop1  : zb->stop1;
            
            /* Set Everything Else                                      */
            zo->startinsert1 = zo->stopinsert1 = ' ';
            zo->startinsert2 = zo->stopinsert2 = ' ';
            zo->chain1 = zo->chain2 = ' ';
            zo->start2 = zo->stop2  = 0;
            zo->mode = ZONE_MODE_SEQUENTIAL;
         }
      }
   }
   
   return(OutputList);
}


/************************************************************************/
/*>ZONE *RenumberZone(ZONE *InputZone, ZONE *OverlapZone)
   ------------------------------------------------------
   This function is called by TrimZones() and returns a zonelist with the 
   fitting zones for a mobile stucture trimed/renumbered based on the 
   sections of fitting zones common to all mobile structures. Each zone in
   the InputZone list is compared to the OverlapZone list and a renumbered
   zonelist is returned.

   20.10.08 Original By: CTP
   04.02.09 Altered if statement from form 'x == y == z' to form 
            'x == z && y == z' as first form generates warning on Mac.
*/
ZONE *RenumberZone(ZONE *InputZone, ZONE *OverlapZone)
{
   ZONE *NewZoneList = NULL;
   ZONE *zi, *zo, *zn;
   
   zi = zo = zn = NULL;
   
   for(zi = InputZone; zi != NULL; NEXT(zi))
   {
      for(zo = OverlapZone; zo != NULL; NEXT(zo))
      {
         if(((zi->start1 <= zo->start1 && zi->stop1 >= zo->start1)  ||
             (zo->start1 <= zi->start1 && zo->stop1 >= zi->start1)) &&
            (zi->mode == ZONE_MODE_SEQUENTIAL && 
             zo->mode == ZONE_MODE_SEQUENTIAL))
         {
            /* If overlap then add to New Zone List                     */
            if(!NewZoneList)
            {
               /* Init zonelist                                         */
               INIT(NewZoneList,ZONE);
               zn = NewZoneList;
            }
            else 
            {
               /* Append zone                                           */
               zn =  NewZoneList;
               LAST(zn);
               ALLOCNEXT(zn,ZONE);
            }
            
            /* Set reference residues                                   */
            zn->start1 = zo->start1;
            zn->stop1  = zo->stop1;
            
            /* Renumber mobile residues                                 */
            zn->start2 = zi->start2 + (zo->start1 - zi->start1);
            zn->stop2  = zi->stop2  - (zi->stop1  - zo->stop1 );
            
            /* Set other values                                         */
            zn->startinsert1 = zn->stopinsert1 = ' ';
            zn->startinsert2 = zn->stopinsert2 = ' ';
            zn->chain1 = zn->chain2 = ' ';
            zn->mode = ZONE_MODE_SEQUENTIAL;
         }
      }
   }
   
   return(NewZoneList);
}

/************************************************************************/
/*>int TrimZones(void)
   -------------------
   This function organizes the fitting zones so that the zones are 
   identical for all mobile structures. For zones derived by pairwise 
   alignment this meant "trimming" the ends of the alignment and ensuring 
   that gaps from one sequence are included across all sequences. This 
   allows meaningful like vs like comparisons when comparing multiple 
   structures and allows setting of mobile structures to the reference 
   structure without having to redefine zones.

   20.10.08 Original By: CTP
   31.10.08 Set error traps for no overlapping/user-defined zones.
   18.12.08 Added conversion to residue numbering then sequential
            numbering to ensure chain breaks included.
 */
int TrimZones(void)
{
   int  i; 
   ZONE *OverlapList = NULL;
   
   /* Check for User-Defined Zones                                      */
   if(!gUserFitZone)
   {
      if(!gQuiet)
         printf("   Error: No user-defined zones found.\n");
      return(0);
   }
   else 
   {
      if(!gQuiet)
         printf("   Finding common zones...\n");
   }
   
   /* Convert to sequential numbering with breaks between chains        */
   if(ConvertAllZones(ZONE_MODE_RESNUM) ||
      ConvertAllZones(ZONE_MODE_SEQUENTIAL))
   {
      printf("   Error: Could not convert zones.\n");
      return(1);
   }
   
   /* Find common residues across all structures                        */
   for(i=1;i<gMultiCount && gUserFitZone;i++)
   {
      ZONE *TempList = NULL;
      if(OverlapList)
      {
         TempList = SetOverlappingZones(OverlapList, gZoneList[i]);
         FREELIST(OverlapList,ZONE);
      }
      else 
      {
         TempList = SetOverlappingZones(gZoneList[0], gZoneList[i]);
      }
      
      OverlapList = TempList;
   }
   
   /* Error Trap: No ZONES                                              */
   if(OverlapList == NULL)
   {
      if(!gQuiet)
         printf("   Warning: No common zones found.\n");
      return(1);
   }
   
   /* Renumber Fit Zones                                                */
   for(i=0;i<gMultiCount && gUserFitZone;i++)
   {
      ZONE *TempList = NULL;
      
      TempList = RenumberZone(gZoneList[i],OverlapList);
      FREELIST(gZoneList[i],ZONE);
      gZoneList[i] = TempList;
   }
   
   /* Set Fitted Flags */
   gFitted      = FALSE;
   gUserFitZone = TRUE;
   
   return(0);
}


/************************************************************************/
/*>int FitStructuresWrapper(void)
   ------------------------------
   Wrapper function for FitStructuresInOrder(). This function finds the 
   RMSD for each mobile structure and sorts the fitting order by RMSD. 
   Although the structures are ordered by RMSD, the function could be 
   expanded to use any other score (for example, fitting by alignment 
   score).

   20.10.08 Original By: CTP
   03.11.08 Turned-off iterative updates when setting fit order.
*/
int FitStructuresWrapper(void)
{
   int mob_i; 
   REAL (*sortrms)[2] = NULL;
   BOOL iterate       = FALSE;
   
   /* Allocate memory                                                   */
   sortrms = malloc(gMultiCount * sizeof *sortrms);
   
   /* Find Sort Order                                                   */
   if(!gQuiet)
      printf("   Setting Fit Order...\n");
   
   iterate = gIterate;
   gIterate = FALSE;
   
   for(mob_i = 0;mob_i < gMultiCount; mob_i++)
   {
      printf("   Mobile %2d ",mob_i + 1);
      sortrms[mob_i][0] = (REAL)mob_i;
      sortrms[mob_i][1] = FitSingleStructure(mob_i, TRUE);
   }
   
   gIterate = iterate;
   
   /* Sort Array to Set Fit Order                                       */
   qsort(sortrms,gMultiCount,sizeof(*sortrms),fit_order_cmp);
   
   /* Fit Stuctures in order                                            */
   FitStructuresInOrder(sortrms);
   
   /* Free Memory                                                       */
   free(sortrms);
   
   return(0);
}

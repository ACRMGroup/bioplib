/************************************************************************/
/**

   \file       deprecated.h
   
   \version    v1.0
   \date       07.05.14
   \brief      Redirect calls to deprecated functions.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2014
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


   Deprecated functions for BiopLib.

   Contains the DEPRECATED() macro which prints a warning message when a 
   call is made to a deprecated function.

   Contains macros with statement expressions that act as wrappers for 
   replacement functions. Calls to deprecated functions are redirected 
   their corresponding replacement functions.



**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   
-  V1.0  07.05.14 Original By: CTP

*************************************************************************/
#ifndef _DEPRECATED_H
#define _DEPRECATED_H

/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/************************************************************************/
/* Defines and macros
*/


/************************************************************************/
/*>DEPRECATED(s, t)
   ----------------
   The DEPRECATED macro gives a warning message if a function is 
   deprecated and indicates the replacement function.
   
   The default option is to give a warning message unless an environment 
   variable, BIOPLIB_DEPRECATED_QUIET, is set.
   
   Alternatively the compile options: -D BIOPLIB_DEPRECATED_CHECK or 
   -D BIOPLIB_DEPRECATED_QUIET will set the DEPRECATED macro to ignore the
   BIOPLIB_DEPRECATED_QUIET environment variable. 
   
   -D BIOPLIB_DEPRECATED_CHECK will display the warning message. 
   -D BIOPLIB_DEPRECATED_QUIET will always silence the warning message.
   
   29.04.14 Original    By: CTP
*/
#if(defined BIOPLIB_DEPRECATED_CHECK || defined BIOPLIB_DEPRECATED_QUIET)
#  ifndef BIOPLIB_DEPRECATED_QUIET
#     define DEPRECATED(s, t)                                            \
               fprintf(stderr,                                           \
                       "This code uses %s which is now deprecated!\n"    \
                       "   Use %s instead\n", (s), (t))
#  else
#      define DEPRECATED(s, t)
#  endif
#else
#  define DEPRECATED(s, t)                                               \
   {                                                                     \
      if(!getenv("BIOPLIB_DEPRECATED_QUIET"))                            \
         fprintf(stderr,                                                 \
                 "This code uses %s which is now deprecated!\n"          \
                 "   Use %s instead\n", (s), (t));                       \
   }
#endif


/************************************************************************/
/*>PDB *FindHetatmResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list, but
  requires the residue is a HETATM record.
  Uses char for chain and insert.

  26.10.11 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindHetatmResidue(). By: CTP
  25.02.14 Added error message. By: CTP
  07.05.14 Converted to macro. By: CTP
*/
#define FindHetatmResidue(pdb, chain, resnum, insert)                    \
({                                                                       \
   char chain_a[2]  = " ",                                               \
        insert_a[2] = " ";                                               \
                                                                         \
   DEPRECATED("FindHetatmResidue()","BiopFindHetatmResidue()");          \
                                                                         \
   chain_a[0]  = chain;                                                  \
   insert_a[0] = insert;                                                 \
                                                                         \
   BiopFindHetatmResidue(pdb, chain_a, resnum, insert_a);                \
})


/************************************************************************/
/*>PDB *FindResidue(PDB *pdb, char chain, int resnum, char insert)
  ----------------------------------------------------------------
  Finds a pointer to the start of a residue in a PDB linked list.
  Uses char for string and insert.

  06.02.96 Original   By: ACRM
  24.02.14 Converted into wrapper for BiopFindResidue(). By: CTP
  25.02.14 Added error message. By: CTP
  07.05.14 Converted to macro. By: CTP

*/

#define FindResidue(pdb, chain, resnum, insert)                          \
({                                                                       \
   char chain_a[2]  = " ",                                               \
        insert_a[2] = " ";                                               \
                                                                         \
   DEPRECATED("FindResidue()","BiopFindResidue()");                      \
                                                                         \
   chain_a[0]  = chain;                                                  \
   insert_a[0] = insert;                                                 \
                                                                         \
   BiopFindResidue(pdb, chain_a, resnum, insert_a);                      \
})


/************************************************************************/
/*>BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, 
                    char stopinsert, char chain, int mode, 
                    PDB **pdb_start, PDB **pdb_stop)
   -------------------------------------------------------------
   Input:   PDB   *pdb        PDB linked list
            int   start       Resnum of start of zone
            char  startinsert Insert code for start of zone
            int   stop        Resnum of end of zone
            char  stopinsert  Insert code for end of zone
            char  chain       Chain name
            int   mode        ZONE_MODE_RESNUM:     Use PDB residue 
                                                    numbers/chain
                              ZONE_MODE_SEQUENTIAL: Use sequential 
                                                    numbering
   Output:  PDB   **pdb_start Start of zone
            PDB   **pdb_stop  End of zone
   Returns: BOOL              OK?

   Finds pointers to the start and end of a zone in a PDB linked list. The
   end is the atom *after* the specified zone

   30.09.92 Original
   17.07.95 Chain name was being ignored in specs like L* (for whole
            of light chain)
   18.08.95 Now handles inserts
   31.07.95 Fixed bug when zone end==chain end
   20.02.01 Changed to -999/-999 for beginning/end of chain rather than -1/-1
   20.03.14 Function deprecated. Converted to wrapper for blFindZonePDB() 
            By: CTP
   07.05.14 Converted to macro. By: CTP
*/
#define FindZonePDB(pdb, start, startinsert, stop, stopinsert, chain,    \
                    mode, pdb_start, pdb_stop)                           \
({                                                                       \
   char startinsert_a[2] = " ",                                          \
        stopinsert_a[2]  = " ",                                          \
        chain_a[2]       = " ";                                          \
                                                                         \
   DEPRECATED("FindZonePDB()","blFindZonePDB()");                        \
                                                                         \
   strncpy(startinsert_a,&startinsert,1);                                \
   strncpy(stopinsert_a ,&stopinsert ,1);                                \
   strncpy(chain_a      ,&chain      ,1);                                \
                                                                         \
   blFindZonePDB(pdb,                                                    \
                 start,                                                  \
                 startinsert_a,                                          \
                 stop,                                                   \
                 stopinsert_a,                                           \
                 chain_a,                                                \
                 mode,                                                   \
                 pdb_start,                                              \
                 pdb_stop);                                              \
})


/************************************************************************/
/*>BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1, 
                  int resnum2, char insert2)
   ----------------------------------------------------------
   Input:   PDB    *p         Pointer to a PDB record
            char   chain      Chain name
            int    resnum1    First residue
            char   insert1    First insert code
            int    resnum2    Second residue
            char   insert2    Second insert code
   Returns: BOOL              Is p in the range specified?

   Checks that atom stored in PDB pointer p is within the specified 
   residue range.

   N.B. This assumes ASCII coding.

   29.03.95 Original    By: ACRM
   08.02.96 Insert residues inside a zone were not handled correctly!
   18.06.96 Added to bioplib from QTree (was called InZone())
   24.02.14 Converted into wrapper for BiopInPDBZone() By: CTP
   25.02.14 Added error message. By: CTP
   07.05.14 Converted to macro. By: CTP
*/
#define InPDBZone(p, chain, resnum1, insert1, resnum2, insert2)          \
({                                                                       \
   char chain_a[2]   = " ",                                              \
        insert1_a[2] = " ",                                              \
        insert2_a[2] = " ";                                              \
                                                                         \
   DEPRECATED("InPDBZone()","BiopInPDBZone()");                          \
                                                                         \
   chain_a[0]   = chain;                                                 \
   insert1_a[0] = insert1;                                               \
   insert2_a[0] = insert2;                                               \
                                                                         \
   BiopInPDBZone(p,chain_a,resnum1,insert1_a,resnum2, insert2_a);        \
})

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#endif

/************************************************************************/
/**

   \file       InPDBZone.c
   
   \version    V1.7
   \date       07.05.14
   \brief      
   
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  30.09.92 Original
-  V1.1  16.06.93 Tidied for book. Mode now a char.
-  V1.2  18.06.96 Added InPDBZone() from QTree program
-  V1.3  19.09.96 Added InPDBZoneSpec()
-  V1.4  24.02.14 Added BiopInPDBZone() By: CTP
-  V1.5  25.02.14 Added error message for InPDBZone(). By: CTP
-  V1.6  02.03.14 Use strcmp to compare inserts.
                  Bugfix for insert residues. By: CTP
-  V1.7  07.05.14 Moved InPDBZone() to deprecated.h By: CTP

*************************************************************************/
/* Includes
*/
#include "SysDefs.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/


/************************************************************************/
/*>BOOL BiopInPDBZone(PDB *p, char *chain, int resnum1, char *insert1, 
                    int resnum2, char *insert2)
   ---------------------------------------------------------------
   Input:   PDB    *p         Pointer to a PDB record
            char   *chain     Chain name
            int    resnum1    First residue
            char   *insert1   First insert code
            int    resnum2    Second residue
            char   *insert2   Second insert code
   Returns: BOOL              Is p in the range specified?

   Checks that atom stored in PDB pointer p is within the specified 
   residue range.

   N.B. This assumes ASCII coding.

   24.02.14 Based on InPDBZone() but takes chain and inserts as stings 
            instead of chars. By: CTP
   02.03.14 Use strcmp to compare inserts.
            Fixed bug handling insert residues where start and finish of 
            zone have same residue number. By: CTP
   
*/
BOOL BiopInPDBZone(PDB *p, char *chain, int resnum1, char *insert1, 
                   int resnum2, char *insert2)
{
   if(CHAINMATCH(p->chain,chain))
   {
      
      /* If residue number is *within* the range, return TRUE           */
      if((p->resnum > resnum1) && (p->resnum < resnum2))
      {
         return(TRUE);
      }
      
      /* If the range has a single residue number, check both inserts   */
      else if((p->resnum == resnum1) && (p->resnum == resnum2))
      {
         if((strcmp(p->insert, insert1) >= 0) &&
            (strcmp(p->insert, insert2) <= 0))
            return(TRUE);
      }
      
      /* If residue number matches ends of range check insert           */
      else if(((p->resnum == resnum1) && 
               (strcmp(p->insert, insert1) >= 0)) ||
              ((p->resnum == resnum2) &&
               (strcmp(p->insert, insert2) <= 0)))
      {
         return(TRUE);
      }
   }
   
   return(FALSE);
}



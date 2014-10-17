/************************************************************************/
/**

   \file       MovePDB.c
   
   \version    V1.3
   \date       07.07.14
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
-  V1.1  01.03.94
-  V1.2  27.02.98 Removed unreachable break from switch()
-  V1.2a 06.01.11 Corrected description
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blMovePDB()
   Moves a PDB record from one linked list to another. from and to should
*/
/************************************************************************/
/* Includes
*/
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"

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
/*>BOOL blMovePDB(PDB *move, PDB **from, PDB **to)
   -----------------------------------------------
*//**

   \param[in]     *move     PDB record to be moved
   \param[in,out] **from    Start of PDB linked list containing record
   \param[in,out] **to      Start of output linked list
   \return                     Success?

   Moves a PDB record from one linked list to another. from and to should
   point to the start of the 2 lists. If the to list hasn't been started,
   to should be NULL. Returns TRUE if moved, FALSE otherwise.

-  13.05.92 Original
-  19.06.92 Changed p=*to, etc. for crappy compilers
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blMovePDB(PDB *move, PDB **from, PDB **to)
{
   PDB *p;
   BOOL ret = FALSE;
   
   if(move != NULL && *from != NULL)
   {
      /* Find the item before move in the *from list                    */
      if(move == (*from))           /* Start of list                    */
      {
         p = NULL;
      }
      else                          /* Middle of list                   */
      {
         /* Move p to item before move                                  */
         for(p = (*from); p->next && p->next != move; NEXT(p)) ;
      }
      
      /* Unlink move from *from                                         */
      if(p)          /* We're moving something in the middle of the list*/
      {
         /* Unlink move                                                 */
         p->next = move->next;
      }
      else           /* We're moving the first one in the list          */
      {
         /* If first in *from list, reset *from list                    */
         *from = move->next;
      }

      /* Add move onto the end of *to                                   */
      move->next = NULL;
      if(*to)
      {
         /* Move p to end of *to list                                   */
         for(p=(*to); p->next; NEXT(p)) ;
         /* Link in move                                                */
         p->next = move;
      }
      else
      {
         /* Initialise *to list                                         */
         *to = move;
      }
      ret = TRUE;
   }
   return(ret);
}


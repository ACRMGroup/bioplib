/************************************************************************/
/**

   \file       IndexPDB.c
   
   \version    V2.2
   \date       19.04.15
   \brief      Create an array of pointers into a PDB linked list
   
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

   IndexPDB() creates an array of pointers to each PDB record in a linked
   list. This allows random access to atoms without having to step through
   the PDB linked list.

**************************************************************************

   Usage:
   ======

   pdb.h must be included before using this routine.

\code
   PDB **indx,
       *pdb;
   int natom;

   indx = IndexPDB(pdb, &natom);
\endcode

**************************************************************************

   Revision History:
   =================
-  V1.0  19.07.90 Original
-  V1.0a 15.02.91 Corrected comments to match new standard.
-  V1.1  01.06.92 ANSIed and documented, FPU condition added
-  V2.0  24.02.94 Completely re-written. Note that the calling format
                  has changed!! NOT BACKWARDLY COMPATIBLE!
-  V2.1  07.07.14 Use bl prefix for functions By: CTP
-  V2.2  19.04.15 Added blIndexAtomNumbersPDB()  By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Manipulating the PDB linked list
   #FUNCTION  blIndexPDB()
   Creates an array of pointers to PDB from a linked list. This is used
   to allow array style access to items in the linked list:
   e.g. (indx[23])->x will give the x coordinate of the 23rd item

   # FUNCTION blIndexAtomNumbersPDB()
   Creates an array of pointers to PDB from a linked list. This is used
   to allow array style access to items in the linked list by atom
   number:
   e.g. (indx[23])->x will give the x coordinate of atom number 23

*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"

/************************************************************************/
/*>PDB **blIndexPDB(PDB *pdb, int *natom)
   --------------------------------------
*//**

   \param[in]     *pdb        Pointer to the start of a PDB linked list.
   \param[out]    *natom      Number of atoms in the PDB linked list.
   \return                    An array of pointers to the PDB records.
                              NULL if unable to allocate memory.

   Creates an array of pointers to PDB from a linked list. This is used
   to allow array style access to items in the linked list:
   e.g. (indx[23])->x will give the x coordinate of the 23rd item

-  19.07.90 Original
-  01.06.92 ANSIed and documented.
-  24.02.94 Re-written. Now allocates and returns the index.
-  07.07.14 Use bl prefix for functions By: CTP
*/
PDB **blIndexPDB(PDB *pdb, int *natom)
{
   PDB *p,
       **indx;
   int i=0;

   /* Count the number of entries                                       */
   for(p=pdb, i=0; p!=NULL; NEXT(p)) i++;
   *natom = i;

   /* Allocate memory for the index array                               */
   if((indx = (PDB **)malloc((i+1) * sizeof(PDB *)))==NULL)
      return(NULL);
   

   for(p=pdb, i=0; p!=NULL; NEXT(p))
      indx[i++] = p;

   indx[i] = NULL;
   
   return(indx);
}

/************************************************************************/
/*>PDB **blIndexAtomNumbersPDB(PDB *pdb, int *indexSize)
   -----------------------------------------------------
*//**
   \param[in]   *pdb         PDB linked list
   \param[out]  *indexSize   Index size
   \return                   malloc'd array of PDB pointers indexed
                             by atom number

   Creates an array of pointers to PDB from a linked list. This is used
   to allow array style access to items in the linked list by atom
   number:
   e.g. (indx[23])->x will give the x coordinate of atom number 23

-  19.04.15 Original   By: ACRM
*/
PDB **blIndexAtomNumbersPDB(PDB *pdb, int *indexSize)
{
   PDB *p,
       **index = NULL;
   int maxAtnum = 0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atnum > maxAtnum)
         maxAtnum = p->atnum;
   }

   *indexSize = maxAtnum + 1;

   if((index=(PDB **)calloc((*indexSize), sizeof(PDB *)))!=NULL)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         index[p->atnum] = p;
      }
   }
   
   return(index);
}



/************************************************************************/
/**

   \file       aalist.c
   
   \version    V3.1
   \date       07.07.14
   \brief      Amino acid linked lists.
   
   \copyright  (c) UCL / Dr. Andrew C.R. Martin 2006-2014
   \author     Dr. Andrew C.R. Martin
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
   #include "bioplib/aalist.h" to define AA datatype.
 
**************************************************************************

   Revision History:
   =================
-  V1.0  21.08.06 Original   By: ACRM
-  V3.0  06.11.08 Incorporated into ProFit V3 By: CTP
-  V3.0  18.02.09 Moved to bioplib. By: CTP
-  V3.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP Sequence manipulation

   #FUNCTION blInsertNextResiduesInAAList()
   Inserts a set of identical residues after the current position in
   the linked list.  The returned value is the last residue which has
   been inserted so this can be called again on the returned aa to
   insert another aa

   #FUNCTION blInsertNextResidueInAAList()
   Inserts a residues after the current position in the linked list.
   The returned value is the residue which has been inserted so this
   can be called again on the returned aa to insert another aa

   #FUNCTION blBuildSeqFromAAList()
   Converts the linked list back into a string which is malloc'd

   #FUNCTION blInsertResidueInAAListAt()
   Inserts a residue after the specified position in the
   list. Residues are numbered from 1. If the position is > length of
   the list then the residue will be added at the end. If the position
   is zero, it will be at the start of the list in which case the
   return value for the list will be different from the input value.

   #FUNCTION blInsertResiduesInAAListAt()
   Inserts a set of residues after the specified position in the
   list. Residues are numbered from 1. If the position is > length of
   the list then the residue will be added at the end. If the position
   is zero, it will be at the start of the list in which case the
   return value for the list will be different from the input value.

   #FUNCTION blBuildAAList()
   Converts a sequence string into a linked list

   #FUNCTION blFindAAListOffsetByResnum()
   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   returns the position of that residue in the list (numbered from 1)

   #FUNCTION blFindAAListItemByResnum()
   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   returns a pointer to that item in the list.

   #FUNCTION blSetAAListFlagByResnum()
   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   sets the flag in that item in the linked list

   #FUNCTION blBuildFlagSeqFromAAList()
   Builds a sequence string with blanks except where the flag in the
   sequence structure is set. At these positions the character specified
   in ch is used instead.

   #FUNCTION blGetAAListLen()
   Returns the number of items in the linked list
*/
/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "macros.h"
#include "aalist.h"

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
/*>AA *blInsertNextResiduesInAAList(AA *a, char res, int nres)
   -----------------------------------------------------------
*//**

   \param[in]   *a    Sequence linked list
   \param[in]   res   Residue to insert
   \param[in]   nres  Number of residues to insert
   \return            Pointer to the residue that has just been
                      inserted

   Inserts a set of identical residues after the current position in
   the linked list.  The returned value is the last residue which has
   been inserted so this can be called again on the returned aa to
   insert another aa

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blInsertNextResiduesInAAList(AA *a, char res, int nres)
{
   int i;
   for(i=0; i<nres; i++)
   {
      a = blInsertNextResidueInAAList(a, res);
   }
   return(a);
}


/************************************************************************/
/*>AA *blInsertNextResidueInAAList(AA *a, char res)
   ------------------------------------------------
*//**

   \param[in]   *a    Sequence linked list
   \param[in]   res   Residue to insert
   \return            Pointer to the residue that has just been
                      inserted

   Inserts a residues after the current position in the linked list.
   The returned value is the residue which has been inserted so this
   can be called again on the returned aa to insert another aa

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blInsertNextResidueInAAList(AA *a, char res)
{
   AA *b = NULL;

   if(a!=NULL)
   {
      INITPREV(b, AA);
      if(b==NULL)
      {
         FREELIST(a, AA);
         return(NULL);
      }
      b->seqnum = (-1);
      b->res = res;

      b->next = a->next;
      b->prev = a;
      if(b->next != NULL)
         b->next->prev = b;
      a->next = b;
   }
   
   return(b);
}


/************************************************************************/
/*>char *blBuildSeqFromAAList(AA *aa)
   ----------------------------------
*//**

   \param[in]   *aa       Sequence linked list
   \return                Sequence as a string (malloc'd)

   Converts the linked list back into a string which is malloc'd

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blBuildSeqFromAAList(AA *aa)
{
   AA *a;
   char *seq=NULL;
   int count=0;
   
   count = blGetAAListLen(aa);
   if((seq=(char *)malloc((1+count)*sizeof(char)))!=NULL)
   {
      count = 0;
      for(a=aa; a!=NULL; NEXT(a))
      {
         seq[count++] = a->res;
      }
      seq[count] = '\0';
   }
   return(seq);
}


/************************************************************************/
/*>AA *blInsertResidueInAAListAt(AA *aa, char res, int pos)
   --------------------------------------------------------
*//**

   \param[in]   *aa         Sequence linked list
   \param[in]   res         Residue to insert
   \param[in]   pos         Position at which to insert (from 1...)
   \return                  Updated sequence linked list

   Inserts a residue after the specified position in the
   list. Residues are numbered from 1. If the position is > length of
   the list then the residue will be added at the end. If the position
   is zero, it will be at the start of the list in which case the
   return value for the list will be different from the input value.

-  21.08.06 Original   By: ACRM
-  18.06.08 Set inserted residue's flag to FALSE.  By: CTP
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blInsertResidueInAAListAt(AA *aa, char res, int pos)
{
   AA *a, *b;
   int count=0;
   
   for(a=aa, count=0; a!=NULL && count<pos-1; count++, NEXT(a));

   if(a!=NULL)
   {
      INITPREV(b, AA);
      if(b==NULL)
      {
         FREELIST(aa, AA);
         return(NULL);
      }
      b->seqnum = (-1);
      b->res = res;
      b->flag = FALSE;

      if(pos==0)
      {
         b->next = aa;
         b->prev = NULL;
         aa->prev = b;
         aa = b;
      }
      else
      {
         b->next = a->next;
         b->prev = a;
         if(b->next != NULL)
            b->next->prev = b;
         a->next = b;
      }
   }
   else     /* Append to end of sequence                                */
   {
      a = aa;
      LAST(a);
      ALLOCNEXTPREV(a, AA);
      if(a==NULL)
      {
         FREELIST(aa, AA);
         return(NULL);
      }
      a->seqnum = (-1);
      a->res = res;
      a->flag = FALSE;
   }
   
   return(aa);
}


/************************************************************************/
/*>AA *blInsertResiduesInAAListAt(AA *aa, char res, int nres, int pos)
   -------------------------------------------------------------------
*//**

   \param[in]   *aa         Sequence linked list
   \param[in]   res         Residue to insert
   \param[in]   nres        Number of residues to insert
   \param[in]   pos         Position at which to insert (from 1...)
   \return                  Updated sequence linked list

   Inserts a set of residues after the specified position in the
   list. Residues are numbered from 1. If the position is > length of
   the list then the residue will be added at the end. If the position
   is zero, it will be at the start of the list in which case the
   return value for the list will be different from the input value.

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blInsertResiduesInAAListAt(AA *aa, char res, int nres, int pos)
{
   int i;
   for(i=0; i<nres; i++)
   {
      aa = blInsertResidueInAAListAt(aa, res, pos);
   }
   return(aa);
}


/************************************************************************/
/*>AA *blBuildAAList(char *seq)
   ----------------------------
*//**

   \param[in]   *seq            The sequence as a string
   \return                      A linked list representation

   Converts a sequence string into a linked list

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blBuildAAList(char *seq)
{
   AA *aa = NULL, 
      *a  = NULL;
   int seqnum = 1;
   
   while(*seq)
   {
      if(aa == NULL)
      {
         INITPREV(aa, AA);
         a = aa;
      }
      else
      {
         ALLOCNEXTPREV(a, AA);
      }
      if(a==NULL)
      {
         FREELIST(aa, AA);
         return(NULL);
      }
      a->res    = *(seq++);
      a->seqnum = seqnum++;
      a->flag   = FALSE;
   }
   
   return(aa);
}

   
/************************************************************************/
/*>int blFindAAListOffsetByResnum(AA *aa, int resnum)
   --------------------------------------------------
*//**

   \param[in]   *aa           Sequence linked list
   \param[in]   resnum        Residue number
   \return                    Linked list offset

   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   returns the position of that residue in the list (numbered from 1)

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blFindAAListOffsetByResnum(AA *aa, int resnum)
{
   int count=1;
   AA  *a;
   
   for(a=aa; a!=NULL; NEXT(a))
   {
      if(a->seqnum == resnum)
         break;
      count++;
   }

   if(a==NULL)
      return(-1);

   return(count);
}


/************************************************************************/
/*>AA *blFindAAListItemByResnum(AA *aa, int resnum)
   ------------------------------------------------
*//**

   \param[in]   *aa           Sequence linked list
   \param[in]   resnum        Residue number
   \return                    Linked list item

   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   returns a pointer to that item in the list.

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
AA *blFindAAListItemByResnum(AA *aa, int resnum)
{
   AA  *a;
   
   for(a=aa; a!=NULL; NEXT(a))
   {
      if(a->seqnum == resnum)
         break;
   }

   return(a);
}


/************************************************************************/
/*>void blSetAAListFlagByResnum(AA *aa, int resnum)
   ------------------------------------------------
*//**

   \param[in]   *aa           Sequence linked list
   \param[in]   resnum        Residue number
   
   Searches the linked list of the specified resnum (i.e. the original
   residue number in the sequence before any insertions were made) and 
   sets the flag in that item in the linked list

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blSetAAListFlagByResnum(AA *aa, int resnum)
{
   AA *a;
   if((a = blFindAAListItemByResnum(aa, resnum))!=NULL)
      a->flag = TRUE;
}


/************************************************************************/
/*>char *blBuildFlagSeqFromAAList(AA *aa, char ch)
   -----------------------------------------------
*//**

   \param[in]   *aa        Sequence linked list
   \param[in]   ch         Character to use in the sequence
   \return                 Sequence string (malloc'd)

   Builds a sequence string with blanks except where the flag in the
   sequence structure is set. At these positions the character specified
   in ch is used instead.

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blBuildFlagSeqFromAAList(AA *aa, char ch)
{
   AA *a;
   char *seq=NULL;
   int count=0;

   count = blGetAAListLen(aa);
   if((seq=(char *)malloc((1+count)*sizeof(char)))!=NULL)
   {
      count = 0;
      for(a=aa; a!=NULL; NEXT(a))
      {
         seq[count++] = ((a->flag)?ch:' ');
      }
      seq[count] = '\0';
   }
   return(seq);
}


/************************************************************************/
/*>int blGetAAListLen(AA *aa)
   --------------------------
*//**

   \param[in]   *aa        Sequence linked list
   \return                 Length of sequence linked list

   Returns the number of items in the linked list

-  21.08.06 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blGetAAListLen(AA *aa)
{
   AA *a;
   int count;
   for(a=aa, count=0; a!=NULL; count++, NEXT(a));
   return(count);
}

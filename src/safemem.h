/************************************************************************/
/**

   \file       safemem.h
   
   \version    V1.4
   \date       14.08.14
   \brief      Safe malloc()/free() routines which check for array 
               overflow on free.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1995-2014
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
   N.B. This header file should be included *after* macros.h

**************************************************************************

   Revision History:
   =================
-  V1.0  23.06.95 Original
-  V1.1  03.07.06 Added safeleaks() prototype
-  V1.2  07.07.14 Use bl prefix for functions By: CTP
-  V1.3  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP
-  V1.4  14.08.14 Moved deprecated function prototypes to deprecated.h 
                  Use blSafefree instead of safefree for macros.
                  By: CTP

*************************************************************************/
#ifndef _SAFEMEM_H
#define _SAFEMEM_H

/* Includes                                                             */
#include "SysDefs.h"

/* Prototypes                                                           */
void *blSafemalloc(int nbytes);
BOOL blSafefree(void *ptr);
void blSafeleaks(void);


/************************************************************************/
/* Include deprecated functions                                         */
#define _SAFEMEM_H_DEPRECATED
# include "deprecated.h" 
/************************************************************************/


/************************************************************************/
/* Undefine memory macros defined by macros.h                           */
#ifdef _MACROS_H
#undef INIT
#undef INITPREV
#undef ALLOCNEXT
#undef ALLOCNEXTPREV
#undef FREELIST
#undef DELETE
#endif

/* Redefine macros to use safe versions of malloc()/free()              */
#define INIT(x,y) do { x=(y *)blSafemalloc(sizeof(y));                   \
                    if(x != NULL) x->next = NULL; } while(0)
#define INITPREV(x,y) do { x=(y *)blSafemalloc(sizeof(y));               \
                       if(x != NULL) {x->next = NULL; x->prev = NULL;} } \
                      while(0)
#define ALLOCNEXT(x,y) do { (x)->next=(y *)blSafemalloc(sizeof(y));      \
                         if((x)->next != NULL) { (x)->next->next=NULL; } \
                         NEXT(x); } while(0)
#define ALLOCNEXTPREV(x,y) do { (x)->next=(y *)blSafemalloc(sizeof(y));  \
                             if((x)->next != NULL)                       \
                             { (x)->next->prev = (x);                    \
                               (x)->next->next=NULL; }                   \
                               NEXT(x);} while(0)
/* FREELIST takes 2 parameters:
   y: name of list
   z: type of list
*/
#define FREELIST(y,z)   while((y)!=NULL)                                 \
                        {  z *_freelist_macro_q;                         \
                           _freelist_macro_q = (y)->next;                \
                           blSafefree((char *)(y));                      \
                           (y) = _freelist_macro_q;                      \
                        }
#define ORDFREELIST(y,z)   while((y)!=NULL)                              \
                        {  z *_freelist_macro_q;                         \
                           _freelist_macro_q = (y)->next;                \
                           free((char *)(y));                            \
                           (y) = _freelist_macro_q;                      \
                        }

/*>DELETE(start, item, type)
   -------------------------
   Deletes item from a linked list
   start may be modified if item is the first in the list
   item is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N time
   to delete N items. If start or item is NULL, does nothing

-  16.02.95 Original    By: ACRM
*/
#define DELETE(x, y, z)                                                  \
do {                                                                     \
   z *_delete_macro_p,                                                   \
     *_delete_macro_prev = NULL,                                         \
     *_delete_macro_temp,                                                \
     *_delete_macro_temp2;                                               \
   if((x)!=NULL && (y)!=NULL)                                            \
   {                                                                     \
      for(_delete_macro_p=(x);                                           \
          _delete_macro_p!=NULL;                                         \
          NEXT(_delete_macro_p))                                         \
      {                                                                  \
         if(_delete_macro_p == (y))                                      \
         {                                                               \
            _delete_macro_temp2 = (y)->next;                             \
            if(_delete_macro_prev == NULL)                               \
            {                                                            \
               _delete_macro_temp = (x)->next;                           \
               blSafefree(x);                                            \
               (x) = _delete_macro_temp;                                 \
               break;                                                    \
            }                                                            \
            else                                                         \
            {                                                            \
               _delete_macro_prev->next = _delete_macro_p->next;         \
               blSafefree(_delete_macro_p);                              \
            }                                                            \
         }                                                               \
         _delete_macro_prev = _delete_macro_p;                           \
      }                                                                  \
      (y) = _delete_macro_temp2;                                         \
   }                                                                     \
}  while(FALSE)

#define SAFEINIT(x,y) INIT(x,y)
#define SAFEINITPREV(x,y) INITPREV(x,y)
#define SAFEALLOCNEXT(x,y) ALLOCNEXT(x,y)
#define SAFEALLOCNEXTPREV(x,y) ALLOCNEXTPREV(x,y)
#define SAFEFREELIST(y,z) FREELIST(y,z)
#define SAFEDELETE(x, y, z) DELETE(x, y, z)

#endif

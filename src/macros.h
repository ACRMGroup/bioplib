/************************************************************************/
/**

   \file       macros.h
   
   \version    V2.25
   \date       04.11.15
   \brief      Useful macros
   
   \copyright  (c) Dr. Andrew C.R. Martin / UCL 1991-2015
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

   If not Amiga defines abs().
   Defines max(), min() and PI if not done.
   Defines list handling macros.
   Defines newline() and toggle() macros.
   
**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  06.02.91 Original
-  V1.1  15.02.91 Moved PI definition to non-Amiga's only
-  V1.2  21.03.91 Added RANGECHECK
-  V1.3  06.09.91 Added DIST, DISTSQ and Vec3f
-  V1.4  09.09.91 Fixed multi-command macros with {}
-  V1.5  24.01.92 Fixed for 32 bit addresses and added malloc checks.
-  V1.6  03.04.92 Small change to ALLOCNEXT and ALLOCNEXTPREV, so
                  will do a NEXT() even if malloc() fails.
-  V1.7  06.05.92 Added TERMINATE()
-  V1.8  06.07.92 Added MAX(), MIN() and ABS()
-  V1.9  22.07.92 Fixed ABS()
-  V1.10 28.09.92 Added TRUE & FALSE and UPPER()
-  V1.11 03.11.92 Changed TOGGLE and newline is now upper case.   
-  V1.12 16.11.92 Added KILLLEADSPACES()
-  V1.13 18.11.92 Fixed UPPER() for MicrosoftC which returns strlen()
                  as unsigned
-  V1.14 20.11.92 ABS() now uses 0 rather than 0.0, so we don't
                  try to use floats with ints...
-  V2.0  24.11.92 Removed all small letter macros
-  V2.1  12.07.93 Added double include check, moved math definitions
                  to MathType.h and added LOWER()
-  V2.2  07.10.93 UPPER() and LOWER() check case first for ESV
                  compatibility
-  V2.3  23.05.94 Added D(BUG)
-  V2.4  14.07.94 Added do{}while(0) bracketing of all multi-line macros
-  V2.5  21.11.94 ABS, MAX and MIN check that they're not already defined
-  V2.6  16.02.95 Added DELETE()
-  V2.7  21.02.95 Updated some internal variable names
-  V2.8  02.08.95 Added TESTINARRAY(), FINDINARRAY(), 
                  SET(), UNSET() and ISSET()
-  V2.9  20.11.95 Added TERMAT()
-  V2.10 06.02.96 Added KILLTRAILSPACES()
-  V2.11 14.06.96 Added PROMPT()
-  V2.12 23.07.96 Added PADMINTERM()
-  V2.13 19.09.96 Include ctype for UPPER() etc
-  V2.14 13.03.99 Added DELETEDOUBLE()
-  V2.15 01.03.01 Added DOTIFY() DEDOTIFY() PADCHARMINTERM() SUBSCHAR()
-  V2.16 25.01.06 Added FINDPREV()
-  V2.17 10.04.08 Fixed bug in DELETE() - the break was not properly
                  stopping prev from being changed
-  V2.18 29.04.14 Added DEPRECATED()   By: CTP
-  V2.19 07.05.14 Moved DEPRECATED() to deprecated.h  By: CTP
-  V2.20 07.07.14 Use bl prefix for functions - change padterm() to 
                  blPadterm() By: CTP
-  V2.21 24.07.14 Initialize list pointers for DELETE macro. By: CTP
-  V2.22 25.02.15 LAST() is now safe if the pointer is NULL
-  V2.23 26.06.15 Added STRNCPYNOSPACES(out, in, mx)
-  V2.24 28.08.15 Added FREE()
-  V2.25 04.11.15 Added FCLOSE()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP Handling linked lists

   #FUNCTION INIT(x,y)
   Macro: Initialise list of name x and type y. Set x->next to NULL
   #FUNCTION INITPREV(x,y)        
   Macro: Initialise list of name x and type y. 
          Set x->next and x->prev to NULL
   #FUNCTION NEXT(x)              
   Macro: Step on in linked list
   #FUNCTION PREV(x)              
   Macro: Step back in doubly linked list
   #FUNCTION ALLOCNEXT(x,y)       
   Macro: Allocate next item in list and step on
   #FUNCTION ALLOCNEXTPREV(x,y)   
   Macro: Allocate next item in doubly linked list and step on. 
   #FUNCTION LAST(x)              
   Macro: Move to end of list
   #FUNCTION FREELIST(y,z)        
   Macro: Free list y of type z
   #FUNCTION FREE(x)
   Macro: Free memory if non-NULL and set the variable to NULL
   #FUNCTION DELETE(lst,itm,type) 
   Macro: Deletes (itm) from linked list (lst) of type (type)
   #FUNCTION FINDPREV(p, start, q) 
   Set p to item in linked list start before q
   #FUNCTION DELETEDOUBLE(lst,itm,type) 
   Macro: Deletes (itm) from a doubly linked list (lst) of type (type)

   #SUBGROUP Miscellaneous
   #FUNCTION NEWLINE              
   Macro: Print a newline character to stdout
   #FUNCTION RANGECHECK(x,y,z)    
   Macro: Return x constrained to range y to z

   #SUBGROUP Flags and bitwise operations
   #FUNCTION TOGGLE(x)            
   Macro: Toggle a flag
   #FUNCTION SET(x,y)             
   Macro: Sets bit y (a hex value) in variable x
   #FUNCTION UNSET(x,y)           
   Macro: Clears bit y (a hex value) in variable x
   #FUNCTION ISSET(x,y)           
   Macro: Tests bit y (a hex value) in variable x

   #SUBGROUP String handling
   #FUNCTION TERMINATE(x)         
   Macro: Terminate a string at the first \n
   #FUNCTION UPPER(x)             
   Macro: Converts a string to upper case
   #FUNCTION KILLLEADSPACES(x,y)  
   Macro: Makes x a pointer into string y after any spaces or tabs.
   #FUNCTION TERMAT(x,y)          
   Macro: Terminates character string x at first character y
   #FUNCTION KILLTRAILSPACES(x)   
   Macro: Terminate string to remove any trailing white space
   #FUNCTION PADMINTERM(str,len)  
   Macro: Pads a string to len chars only if it is shorter
   #FUNCTION PADCHARMINTERM(str,char,len) 
   Macro: Pads a string to len chars with specified
          character only if it is shorter
   #FUNCTION DOTIFY(str)          
   Macro: Replace ' ' with '.' in string
   #FUNCTION DEDOTIFY(str)        
   Macro: Replace '.' with ' ' in string
   #FUNCTION STRNCPYNOSPACES(out, in, maxlen)
   Macro: Like strncpy() but skips spaces

   #SUBGROUP Maths
   #FUNCTION MAX(x,y)             
   Macro: max() as macro
   #FUNCTION MIN(x,y)             
   Macro: min() as macro
   #FUNCTION ABS(x,y)             
   Macro: abs() as macro

   #SUBGROUP Debugging
   #FUNCTION D(BUG)               
   Macro: Prints the BUG string if DEBUG is defined first

   #SUBGROUP Array handling
   #FUNCTION TESTINARRAY(x,l,y,r) 
   Macro: Tests whether value (y) is in array (x) if length (l) returning the result in (r)
   #FUNCTION FINDINARRAY(x,l,y,r) 
   Macro: Finds value (y) is in array (x) if length
      (l) returning the offset in (r). Offset is -1 if not found

   #SUBGROUP User interaction
   #FUNCTION PROMPT(fp,x)         
   Macro: Issue a prompt to stdout if fp is a terminal

   #SUBGROUP File handling
   #FUNCTION FCLOSE(fp)
   Macro: close a file pointer if it is non-NULL and not stdin/out/err.
   Set the file pointer to NULL afterwards

*/
/************************************************************************/
#ifndef _MACROS_H
#define _MACROS_H

/***************************** Includes *********************************/
#include <ctype.h>

/**************************** Definitions *******************************/
#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/***************************** Maths macros *****************************/
#define RANGECHECK(x,y,z) ((x)<(y)) ? (y) : ((x)>(z)) ? (z) : (x)
#define DISTSQ(a,b) (((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                     ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                     ((a)->z - (b)->z) * ((a)->z - (b)->z))
#define DIST(a,b) sqrt(((a)->x - (b)->x) * ((a)->x - (b)->x) + \
                       ((a)->y - (b)->y) * ((a)->y - (b)->y) + \
                       ((a)->z - (b)->z) * ((a)->z - (b)->z))
#ifndef ABS
#define ABS(x)   (((x)<0)   ? (-(x)) : (x))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b)) ? (a)    : (b))
#define MIN(a,b) (((a)<(b)) ? (a)    : (b))
#endif

/***************************** List macros ******************************/
#define INIT(x,y) do { x=(y *)malloc(sizeof(y)); \
                    if(x != NULL) x->next = NULL; } while(0)
#define INITPREV(x,y) do { x=(y *)malloc(sizeof(y));\
                        if(x != NULL) {x->next=NULL; x->prev=NULL;} } \
                      while(0)
#define NEXT(x) (x)=(x)->next
#define PREV(x) (x)=(x)->prev
#define ALLOCNEXT(x,y) do { (x)->next=(y *)malloc(sizeof(y));\
                         if((x)->next != NULL) { (x)->next->next=NULL; }\
                         NEXT(x); } while(0)
#define ALLOCNEXTPREV(x,y) do { (x)->next=(y *)malloc(sizeof(y));\
                             if((x)->next != NULL)\
                             { (x)->next->prev = (x); \
                               (x)->next->next=NULL; }\
                               NEXT(x);} while(0)
#define LAST(x)   while(((x)!=NULL) && ((x)->next != NULL)) NEXT(x)
/* FREELIST takes 2 parameters:
   y: name of list
   z: type of list
*/
#define FREELIST(y,z)   while((y)!=NULL) \
                        {  z *_freelist_macro_q; \
                           _freelist_macro_q = (y)->next; \
                           free((char *)(y)); \
                           (y) = _freelist_macro_q; \
                        }
#define FREE(x) if((x)!=NULL) free(x); (x) = NULL

/*>DELETE(start, item, type)
   -------------------------
   Deletes (item) from a linked list.
   (start) will be modified if (item) is the first in the list.
   (item) is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N times
   to delete N items. If (start) or (item) is NULL, does nothing

-  16.02.95 Original    By: ACRM
-  10.04.08 Fixed position of break. By: CTP
-  24.07.14 Initialize list pointers. By: CTP
*/
#define DELETE(x, y, z)                                                  \
do {                                                                     \
   z *_delete_macro_p = NULL,                                            \
     *_delete_macro_prev = NULL,                                         \
     *_delete_macro_temp = NULL,                                         \
     *_delete_macro_temp2 = NULL;                                        \
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
               free(x);                                                  \
               (x) = _delete_macro_temp;                                 \
            }                                                            \
            else                                                         \
            {                                                            \
               _delete_macro_prev->next = _delete_macro_p->next;         \
               free(_delete_macro_p);                                    \
            }                                                            \
            break;                                                       \
         }                                                               \
         _delete_macro_prev = _delete_macro_p;                           \
      }                                                                  \
      (y) = _delete_macro_temp2;                                         \
   }                                                                     \
}  while(FALSE)


/*>DELETEDOUBLE(start, item, type)
   -------------------------------
   Deletes (item) from a doubly linked list.
   (start) will be modified if (item) is the first in the list.
   (item) is returned as the pointer to the next item in the list (i.e.
   as item->next). One can therefore simply call the routine N times
   to delete N items. If (start) or (item) is NULL, does nothing

-  13.03.99 Original    By: ACRM
*/
#define DELETEDOUBLE(s, x, y)                                            \
        do { y *_deleteandnext_macro_temp;                               \
             if(((s)!=NULL) && ((x)!=NULL))                              \
             {  if((x)==(s)) (s) = (x)->next;                            \
                _deleteandnext_macro_temp = (x)->next;                   \
                if((x)->prev != NULL) (x)->prev->next = (x)->next;       \
                if((x)->next != NULL) (x)->next->prev = (x)->prev;       \
                free(x);                                                 \
                (x) = _deleteandnext_macro_temp;                         \
        }  } while(0)

/*>FINDPREV(ptr, start, item)
   --------------------------
   Searches a linked list beginning at (start) to find the item which
   preceeds (item). Its address is put into (ptr). If (item) is the
   same as (start) or (item) is not found, then the routine returns
   NULL in (ptr)
   This is used when wanting to look at the previous item in a singly
   linked list.

-  26.01.06 Original    By: ACRM
*/
#define FINDPREV(p, s, l)                                                \
        do { p = (s);                                                    \
             if((s)==(l))                                                \
             { p = NULL; } else                                          \
             {                                                           \
               while((p != NULL) && (p->next != (l)))                    \
               {  p = p->next;                                           \
           } } } while(0)


/***************************** Misc. macros *****************************/
#define NEWLINE printf("\n")

#define TOGGLE(x) (x) = (x) ? FALSE : TRUE

#define TERMINATE(x) do {  int _terminate_macro_j;                    \
                        for(_terminate_macro_j=0;                     \
                            (x)[_terminate_macro_j];                  \
                            _terminate_macro_j++)                     \
                        {  if((x)[_terminate_macro_j] == '\n')        \
                           {  (x)[_terminate_macro_j] = '\0';         \
                              break;                                  \
                     }  }  }  while(0)
#define TERMAT(x, y) do {  int _termat_macro_j;                       \
                        for(_termat_macro_j=0;                        \
                            (x)[_termat_macro_j];                     \
                            _termat_macro_j++)                        \
                        {  if((x)[_termat_macro_j] == (y))            \
                           {  (x)[_termat_macro_j] = '\0';            \
                              break;                                  \
                     }  }  }  while(0)
#define UPPER(x) do {  int _upper_macro_i;                            \
                    for(_upper_macro_i=0;                             \
                        _upper_macro_i<(int)strlen(x) &&              \
                           (x)[_upper_macro_i];                       \
                        _upper_macro_i++)                             \
                           if(islower((x)[_upper_macro_i]))           \
                              (x)[_upper_macro_i] =                   \
                                 (char)toupper((x)[_upper_macro_i]);  \
                    }  while(0)
#define LOWER(x) do {  int _lower_macro_i;                            \
                    for(_lower_macro_i=0;                             \
                        _lower_macro_i<(int)strlen(x) &&              \
                           (x)[_lower_macro_i];                       \
                        _lower_macro_i++)                             \
                           if(isupper((x)[_lower_macro_i]))           \
                              (x)[_lower_macro_i] =                   \
                                 (char)tolower((x)[_lower_macro_i]);  \
                    }  while(0)
#define KILLLEADSPACES(y,x)                                           \
                 do \
                 {  for((y)=(x); *(y) == ' ' || *(y) == '\t'; (y)++) ; } \
                 while(0)


#define KILLTRAILSPACES(x)                                              \
do {  int _kts_macro_i;                                                 \
      _kts_macro_i = strlen(x) - 1;                                     \
      while(((x)[_kts_macro_i] == ' ' ||                                \
             (x)[_kts_macro_i] == '\t') &&                              \
            _kts_macro_i>=0)                                            \
         (_kts_macro_i)--;                                              \
      (x)[++(_kts_macro_i)] = '\0';                                     \
   }  while(0)


/* Tests for the presence of (y) in array (x) of length (l). The result
   (TRUE or FALSE) is returned in (r)
-  02.08.95 Original
*/
#define TESTINARRAY(x, l, y, r)                                          \
do {                                                                     \
   int _inarray_macro_i;                                                 \
   (r) = FALSE;                                                          \
   if((x)==NULL) break;                                                  \
   for(_inarray_macro_i=0; _inarray_macro_i<(l); _inarray_macro_i++)     \
   {  if((x)[_inarray_macro_i] == (y))                                   \
      {  (r) = TRUE;                                                     \
         break;                                                          \
}  }  } while(FALSE)

/* Finds offset of item (y) in array (x) of length (l). The result
   is returned in (r) which is -1 if item not found
-  02.08.95 Original
*/
#define FINDINARRAY(x, l, y, r)                                          \
do {                                                                     \
   int _inarray_macro_i;                                                 \
   (r) = (-1);                                                           \
   if((x)==NULL) break;                                                  \
   for(_inarray_macro_i=0; _inarray_macro_i<(l); _inarray_macro_i++)     \
   {  if((x)[_inarray_macro_i] == (y))                                   \
      {  (r) = _inarray_macro_i;                                         \
         break;                                                          \
}  }  } while(FALSE)


/* Used just like padterm, but doesn't touch the string if it's already
   longer than len characters
*/
#define PADMINTERM(string, len)                                  \
        do {                                                     \
        if(strlen((string)) < (len)) blPadterm((string), (len)); \
        } while(0)

/************************************************************************/
/*>PADCHARMINTERM(string, char, length)
   ------------------------------------
*//**

   Pads a string to a specified length using char and terminates at that 
   point

-  13.03.99 Original   By: ACRM
*/
#define PADCHARMINTERM(s, c, l)                                          \
do {  int _padminterm_macro_i;                                           \
      if(strlen((s)) < (l))                                              \
      {  for(_padminterm_macro_i=strlen((s));                            \
             _padminterm_macro_i<(l);                                    \
             _padminterm_macro_i++)                                      \
            (s)[_padminterm_macro_i] = (c);                              \
         (s)[(l)] = '\0';                                                \
      }  }  while(0)


/************************************************************************/
/*>DOTIFY(char *str)
   -----------------
*//**

   Macro to replace ' ' in a string with '.'

-  21.04.99 Original   By: ACRM
*/
#define DOTIFY(str)                                                      \
do {                                                                     \
   char *_dotify_macro_chp;                                              \
   _dotify_macro_chp = str;                                              \
   while(*_dotify_macro_chp) {                                           \
      if(*_dotify_macro_chp==' ') *_dotify_macro_chp = '.';              \
      _dotify_macro_chp++;                                               \
}  } while(0)

/************************************************************************/
/*>DEDOTIFY(char *str)
   -------------------
*//**

   Macro to replace '.' in a string with ' '

-  21.04.99 Original   By: ACRM
*/
#define DEDOTIFY(str)                                                    \
do {                                                                     \
   char *_dedotify_macro_chp;                                            \
   _dedotify_macro_chp = str;                                            \
   while(*_dedotify_macro_chp) {                                         \
      if(*_dedotify_macro_chp=='.') *_dedotify_macro_chp = ' ';          \
      _dedotify_macro_chp++;                                             \
}  } while(0)


/************************************************************************/
/*>SUBSCHAR(s, x, y)
   -----------------
*//**

   Substitute character x by character y in string s

-  21.05.99 Original
*/
#define SUBSCHAR(s, x, y)                                                \
do {  char *_subschar_macro_ch = (s);                                    \
      while(*_subschar_macro_ch != '\0')                                 \
      {  if(*_subschar_macro_ch == (x)) *_subschar_macro_ch = (y);       \
         _subschar_macro_ch++;                                           \
   }  }  while(0)


/************************************************************************/
#define STRNCPYNOSPACES(out, in, mx)                                     \
do {  char *_chp;                                                      \
      int _pos = 0;                                                    \
      for(_chp=(in); *_chp != '\0'; _chp++)                            \
      {  if(*_chp != ' ')                                              \
         {   buffer[_pos++] = *_chp;                                   \
            if(_pos >= (mx)) break;                                    \
      }  }                                                             \
      if(_pos < (mx)) buffer[_pos] = '\0';                             \
   } while(0)


/* Bit-wise operators
-  02.08.95 Original
*/
#define SET(x, y)   (x) |= (y)
#define UNSET(x, y) (x) &= (~(y))
#define ISSET(x, y) ((BOOL)((x)&(y)))

/* A version of fclose() that checks the file pointer is non-NULL and
   not standard in/out/err. After closing the file, it sets the file
   pointer to NULL
   04.11.15 Original   By: ACRM
*/
#define FCLOSE(x)                                                        \
   do {                                                                  \
      if(((x) != NULL) && ((x) != stdin) &&                              \
         ((x) != stdout) && ((x) != stderr)) {                           \
         fclose((x));                                                    \
         (x) = NULL;                                                     \
      }                                                                  \
   } while (0)

#ifdef DEBUG
#define D(BUG) fprintf(stderr,"%s",BUG); fflush(stderr)
#else
#define D(BUG)
#endif


/************************** The PROMPT macro ****************************/
/* isatty() is not POSIX                                                */
#ifdef __unix
#  if defined(_POSIX_SOURCE) || !defined(_SVR4_SOURCE)
      #include <unistd.h>
#  endif
#endif

/* Default is just to print a string as a prompt                        */
#define PROMPT(in,x) printf("%s",(x))

/* More intelligent prompts for systems where we know the FILE structure*/
#ifdef __sgi
#  undef PROMPT
#  define PROMPT(in,x) do{if(isatty((in)->_file)) \
                      printf("%s",(x));}while(0)
#endif
#ifdef __linux__
#  undef PROMPT
#  define PROMPT(in,x) do{if(isatty((in)->_fileno)) \
                      printf("%s",(x));}while(0)
#endif

#endif /* _MACROS_H                                                     */


/************************************************************************/
/**

   \file       hash.c
   
   \version    V1.2
   \date       23.08.15
   \brief      Flexible hash functions
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2015
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

   Based very loosely on code from 
   http://www.sourcecodesworld.com/source/show.asp?ScriptID=1188

**************************************************************************

   Usage:
   ======
   Compile with -finline-functions

**************************************************************************

   Revision History:
   =================
-  V1.0   14.05.15  Original   By: ACRM
-  V1.1   03.06.15  Added wrappers to blSetHashValue()
-  V1.2   23.08.15  Cast NAN for int return

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP Hashes / Dictionaries

   #FUNCTION blInitializeHash()
   Initializes a hash structure

   #FUNCTION blFreeHashKeyList()
   Frees a list of hash keys from blGetHashKeyList()

   #FUNCTION blGetHashKeyList()
   Allocates an array containing the hash keys

   #FUNCTION blFreeHash()
   Frees all memory associated with a hash

   #FUNCTION blSetHashValue()
   Sets a hash key:value pair

   #FUNCTION blSetHashValueString()
   Sets a hash key:value pair for a string value

   #FUNCTION blSetHashValueInt()
   Sets a hash key:value pair for an int value

   #FUNCTION blSetHashValueDouble()
   Sets a hash key:value pair for a double value

   #FUNCTION blSetHashValuePointer()
   Sets a hash key:value pair for a general pointer (BPTR) value

   #FUNCTION blSetHashValueChar()
   Sets a hash key:value pair for a char value

   #FUNCTION blGetHashValue()
   Gets a hash value for a specified key

   #FUNCTION blGetHashValueInt()
   Wrapper to blGetHashValue() specifically for ints

   #FUNCTION blGetHashValueDouble()
   Wrapper to blGetHashValue() specifically for doubles

   #FUNCTION blGetHashValueChar()
   Wrapper to blGetHashValue() specifically for chars

   #FUNCTION blGetHashValueString()
   Wrapper to blGetHashValue() specifically for strings

   #FUNCTION blGetHashValuePointer()
   Wrapper to blGetHashValue() specifically for pointers

   #FUNCTION blDumpHash()
   Utility function to dump a hash to a file

   #FUNCTION blHashKeyDefined()
   Tests if a key is defined in the hash

   #FUNCTION blDeleteHashKey()
   Deletes an entry from the hash

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "hash.h"
#include "macros.h"
#include "MathUtil.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static ULONG      hash(char *s, ULONG size);
static _HASHENTRY *lookup(HASHTABLE *hashtable, char *n);
static ULONG      CalculateHashSize(ULONG hashsize);
#ifdef DEBUG
static void       dumpHashTable(FILE *out, HASHTABLE *hashtable);
#endif

/************************************************************************/
/*>HASHTABLE *blInitializeHash(ULONG hashsize)
   -------------------------------------------
*//**
   \param[in]   hashsize  An estimate of the size of the hash (or 0)
   \return                A pointer to a hashtable structure

   The size estimate is used to allocate a hash table. If too large you
   waste memory; if too small lookups will take longer. You can use 0
   for the default hash size of 1001

   A future version will implement dynamic hashing.

-  12.05.15  Original   By: ACRM
*/
HASHTABLE *blInitializeHash(ULONG hashsize)
{
   HASHTABLE *hashtable = NULL;

   /* Allocate the main structure                                       */
   if((hashtable = (HASHTABLE *)malloc(sizeof(HASHTABLE)))==NULL)
      return(NULL);

   /* Allocate the table                                                */
   hashtable->size = CalculateHashSize(hashsize);
   if((hashtable->table = 
       (_HASHENTRY **)malloc(hashtable->size*sizeof(_HASHENTRY)))!=NULL)
   {
      int i;
      for(i=0; i<hashtable->size; i++)
         hashtable->table[i] = NULL;
   }
   else
   {
      free(hashtable);
      hashtable = NULL;
   }

   return(hashtable);
}


/************************************************************************/
/*>void blFreeHashKeyList(char **keylist)
   --------------------------------------
*//**
   \param[in]    keylist    List of keys returned by blGetHashKeyList()

   Frees memory allocated for a list of keys in the hash

-  12.05.15  Original   By: ACRM
*/
void blFreeHashKeyList(char **keylist)
{
   int i;

   if(keylist == NULL)
      return;
   
   for(i=0; keylist[i] != NULL; i++)
      free(keylist[i]);

   free(keylist);
}


/************************************************************************/
/*>char **blGetHashKeyList(HASHTABLE *hashtable)
   ---------------------------------------------
*//**
   \param[in]   hashtable    The hashtable
   \return                   Array of strings representing hash keys

   Gets a list of the keys in the hash. Allocates memory which must be
   freed using blFreeHashKeyList()

-  12.05.15  Original   By: ACRM
*/
char **blGetHashKeyList(HASHTABLE *hashtable)
{
   int        i, j,
              nKeys;
   _HASHENTRY **hashtab;
   char       **keys = NULL;

   if(hashtable==NULL)
      return(NULL);

   hashtab = hashtable->table;

   /* First walk the hash table to see how many keys there are          */
   nKeys = 0;
   for(i=0; i<hashtable->size; i++)
   {
      if(hashtab[i]!=NULL)
      {
         _HASHENTRY *hashEntry;

         for(hashEntry=hashtab[i]; hashEntry!=NULL; NEXT(hashEntry))
            nKeys++;
      }
   }

   /* Allocate an array of this size (+1)                               */
   if((keys = (char **)malloc((nKeys+1)*sizeof(char *)))==NULL)
      return(NULL);

   /* Set the last one to NULL to indicate the end of the list          */
   keys[nKeys] = NULL;
   
   /* Walk the hash table again, copying in the keys                    */
   for(i=0, j=0; i<hashtable->size && j<nKeys; i++)
   {
      if(hashtab[i]!=NULL)
      {
         _HASHENTRY *hashEntry;

         for(hashEntry=hashtab[i]; hashEntry!=NULL; NEXT(hashEntry))
         {
            keys[j] = blStrdup(hashEntry->key);
            /* Free key list if allocation error and return NULL        */
            if(keys[j] == NULL)
            {
               blFreeHashKeyList(keys);
               return(NULL);
            }
            j++;
         }
      }
   }

   return(keys);
}


/************************************************************************/
/*>void blFreeHash(HASHTABLE *hashtable)
   -------------------------------------
*//**
   \param[in]   hashtable    The hash table

   Frees all memory used by the hash

-  12.05.15  Original   By: ACRM
*/
void blFreeHash(HASHTABLE *hashtable)
{
   int        i;
   _HASHENTRY *nextItem,
              **hashtab;

   if(hashtable == NULL)
      return;
   
   hashtab = hashtable->table;

   /* Walk the hash table                                               */
   for(i=0; i<hashtable->size; i++)
   {
      /* If this slot has been used                                     */
      if(hashtab[i]!=NULL)
      {
         _HASHENTRY *hashEntry;

         hashEntry=hashtab[i];

         /* Walk the linked list freeing entries                        */
         while(hashEntry!=NULL)
         {
            nextItem=hashEntry->next;
            if(hashEntry->key   != NULL) free(hashEntry->key);
            if(hashEntry->string != NULL) free(hashEntry->string);
            free(hashEntry);
            hashEntry=nextItem;
         }
      }
   }
   /* Finally free the table and the main hash structure                */
   free(hashtab);
   free(hashtable);
}


/************************************************************************/
/*>BOOL blSetHashValue(HASHTABLE *hashtable, char *key, int type, ...)
   -------------------------------------------------------------------
*//**
   \param[in]  hashtable   The hash table
   \param[in]  key         The hash key
   \param[in]  type        The value type
   \param[in]  ...         The value
   \return                 Success

   Sets a value in the hash. 'type' is one of
-   HASHTYPE_INT
-   HASHTYPE_STRING
-   HASHTYPE_DOUBLE
-   HASHTYPE_CHAR
-   HASHTYPE_POINTER

-  12.05.15  Original   By: ACRM
*/
BOOL blSetHashValue(HASHTABLE *hashtable, char *key, int type, ...)
{
   va_list    ap;
   ULONG      hashIndex;
   _HASHENTRY *hashEntry;
   char       *stringPtr;
   int        intValue,
              charValue;
   double     doubleValue;
   BPTR       ptrValue;
   
   if(hashtable == NULL)
      return(FALSE);
   
   /* If this key is not already in the hash                            */
   if((hashEntry=lookup(hashtable, key))==NULL)
   {
      _HASHENTRY **hashtab = hashtable->table;

      /* Calculate a hash index for this key                            */
      hashIndex = hash(key, hashtable->size);

      /* And allocate a structure to store the key and value            */
      if((hashEntry = (_HASHENTRY *)malloc(sizeof(_HASHENTRY)))==NULL)
         return(FALSE);

      /* Copy in the key                                                */
      if((hashEntry->key = blStrdup(key))==NULL)
         return(FALSE);

      /* Insert this entry at the start of the linked list for this
         hashtable index
      */
      hashEntry->next=hashtab[hashIndex];
      hashtab[hashIndex]=hashEntry;
   }
   else
   {
      /* The key is already in the hash. If it was a string, delete 
         its value             
      */
      if((hashEntry->type == HASHTYPE_STRING) && 
         (hashEntry->string != NULL))
         free(hashEntry->string);
   }

   /* Initialize values in the hash entry                               */
   hashEntry->string      = NULL;
   hashEntry->intValue    = 0;
   hashEntry->doubleValue = 0.0;
   hashEntry->charValue   = '\0';
   hashEntry->ptrValue    = NULL;
   
   
   /* Copy in the new value                                             */
   va_start(ap, type);
   switch(type)
   {
   case HASHTYPE_STRING:
      stringPtr = va_arg(ap, char *);
      if((hashEntry->string=blStrdup(stringPtr))==NULL)
         return(FALSE);
      break;
   case HASHTYPE_INT:
      intValue  = va_arg(ap, int);
      hashEntry->intValue = intValue;
      break;
   case HASHTYPE_CHAR:
      charValue = va_arg(ap, int);
      hashEntry->charValue = (char)charValue;
      break;
   case HASHTYPE_DOUBLE:
      doubleValue = va_arg(ap, double);
      hashEntry->doubleValue = doubleValue;
      break;
   case HASHTYPE_POINTER:
      ptrValue = va_arg(ap, BPTR);
      hashEntry->ptrValue = ptrValue;
      break;
   }
   va_end(ap);

   hashEntry->type = type;

   return(TRUE);
}


/************************************************************************/
/*>BOOL blSetHashValueString(HASHTABLE *hashtable, char *key, char *value)
   -----------------------------------------------------------------------
*//**
   \param[in]    *hashtable    The hash table
   \param[in]    *key          The key for the hash entry
   \param[in]    *value        The value for the hash entry
   \return                     Success

   Set a string value in a hash
   Wrapper to blSetHashValue()

-  03.06.15   Original   By: ACRM
*/
BOOL blSetHashValueString(HASHTABLE *hashtable, char *key, char *value)
{
   return(blSetHashValue(hashtable, key, HASHTYPE_STRING, value));
}


/************************************************************************/
/*>BOOL blSetHashValueInt(HASHTABLE *hashtable, char *key, int value)
   ------------------------------------------------------------------
*//**
   \param[in]    *hashtable    The hash table
   \param[in]    *key          The key for the hash entry
   \param[in]    *value        The value for the hash entry
   \return                     Success

   Set an int value in a hash 
   Wrapper to blSetHashValue()

-  03.06.15   Original   By: ACRM
*/
BOOL blSetHashValueInt(HASHTABLE *hashtable, char *key, int value)
{
   return(blSetHashValue(hashtable, key, HASHTYPE_INT, value));
}


/************************************************************************/
/*>BOOL blSetHashValueDouble(HASHTABLE *hashtable, char *key,
   double value)
   ----------------------------------------------------------   
*//**
   \param[in]    *hashtable    The hash table
   \param[in]    *key          The key for the hash entry
   \param[in]    *value        The value for the hash entry
   \return                     Success

   Set a double value in a hash
   Wrapper to blSetHashValue()

-  03.06.15   Original   By: ACRM
*/
BOOL blSetHashValueDouble(HASHTABLE *hashtable, char *key, double value)
{
   return(blSetHashValue(hashtable, key, HASHTYPE_DOUBLE, value));
}


/************************************************************************/
/*>BOOL blSetHashValuePointer(HASHTABLE *hashtable, char *key, BPTR ptr)
   ---------------------------------------------------------------------
*//**
   \param[in]    *hashtable    The hash table
   \param[in]    *key          The key for the hash entry
   \param[in]    *value        The value for the hash entry
   \return                     Success

   Set a general pointer (BPTR) value in a hash
   Wrapper to blSetHashValue()

-  03.06.15   Original   By: ACRM
*/
BOOL blSetHashValuePointer(HASHTABLE *hashtable, char *key, BPTR ptr)
{
   return(blSetHashValue(hashtable, key, HASHTYPE_POINTER, ptr));
}


/************************************************************************/
/*>BOOL blSetHashValueChar(HASHTABLE *hashtable, char *key, char value)
   --------------------------------------------------------------------
*//**
   \param[in]    *hashtable    The hash table
   \param[in]    *key          The key for the hash entry
   \param[in]    *value        The value for the hash entry
   \return                     Success

   Wrapper to blSetHashValue()

-  03.06.15   Original   By: ACRM
*/
BOOL blSetHashValueChar(HASHTABLE *hashtable, char *key, char value)
{
   return(blSetHashValue(hashtable, key, HASHTYPE_CHAR, value));
}


/************************************************************************/
/*>int blGetHashValueInt(HASHTABLE *hashtable, char *key)
   ------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \return                    The value

   Simple wrapper to blGetHashValue() for extracting integers

-  12.05.15  Original   By: ACRM
-  23.08.15  Added cast of NAN
*/
int blGetHashValueInt(HASHTABLE *hashtable, char *key)
{
   BPTR value;
   if((value = blGetHashValue(hashtable, key, NULL))!=NULL)
      return(*(int *)value);
#ifdef NAN   
   return((int)NAN);
#endif
   return(0);
}


/************************************************************************/
/*>double blGetHashValueDouble(HASHTABLE *hashtable, char *key)
   ------------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \return                    The value

   Simple wrapper to blGetHashValue() for extracting doubles

-  12.05.15  Original   By: ACRM
*/
double blGetHashValueDouble(HASHTABLE *hashtable, char *key)
{
   BPTR value;
   if((value = blGetHashValue(hashtable, key, NULL))!=NULL)
      return(*(double *)value);
#ifdef NAN   
   return(NAN);
#endif
   return(0.0);
}


/************************************************************************/
/*>char blGetHashValueChar(HASHTABLE *hashtable, char *key)
   --------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \return                    The value

   Simple wrapper to blGetHashValue() for extracting characters

-  12.05.15  Original   By: ACRM
*/
char blGetHashValueChar(HASHTABLE *hashtable, char *key)
{
   BPTR value;
   if((value = blGetHashValue(hashtable, key, NULL))!=NULL)
      return(*(char *)value);
   return('\0');
}


/************************************************************************/
/*>char *blGetHashValueString(HASHTABLE *hashtable, char *key)
   -----------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \return                    The value

   Simple wrapper to blGetHashValue() for extracting strings

-  12.05.15  Original   By: ACRM
*/
char *blGetHashValueString(HASHTABLE *hashtable, char *key)
{
   return((char *)blGetHashValue(hashtable, key, NULL));
}


/************************************************************************/
/*>BPTR blGetHashValuePointer(HASHTABLE *hashtable, char *key)
   -----------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \return                    The value

   Simple wrapper to blGetHashValue() for extracting pointers

-  12.05.15  Original   By: ACRM
*/
BPTR blGetHashValuePointer(HASHTABLE *hashtable, char *key)
{
   return((BPTR)blGetHashValue(hashtable, key, NULL));
}


/************************************************************************/
/*>BPTR blGetHashValue(HASHTABLE *hashtable, char *key, int *type)
   ---------------------------------------------------------------
*//**
   \param[in]    hashtable    The hash table
   \param[in]    key          The hash key
   \param[out]   type         The type for that item (or NULL)
   \return                    Pointer to the value

   Obtain an entry from the hash. Type 'type' parameter can be NULL
   if you are sure you know what the type is and don't want to check.
   The value is returned as a pointer. If the value is itself a pointer
   or a string then this is what you want, but you will need to cast
   it appropriately. If it is an int, double or char then you need to 
   cast and dereference this:

   string:  (char *)value
   pointer: (my-type *)value
   int:     *(int *)value
   double:  *(double *)value
   char:    *(char *)value

-  12.05.15  Original   By: ACRM
*/
BPTR blGetHashValue(HASHTABLE *hashtable, char *key, int *type)
{
   _HASHENTRY *hashEntry;

   if(hashtable == NULL)
      return(FALSE);
   
   if((hashEntry = lookup(hashtable, key))==NULL)
   {
      return(NULL);
   }
   else
   {
      if(type != NULL)
         *type = hashEntry->type;

      switch(hashEntry->type)
      {
      case HASHTYPE_INT:
         return((BPTR)(&(hashEntry->intValue)));
         break;
      case HASHTYPE_DOUBLE:
         return((BPTR)(&(hashEntry->doubleValue)));
         break;
      case HASHTYPE_CHAR:
         return((BPTR)(&(hashEntry->charValue)));
         break;
      case HASHTYPE_STRING:
         return((BPTR)(hashEntry->string));
         break;
      case HASHTYPE_POINTER:
         return((BPTR)(hashEntry->ptrValue));
         break;
      }
   }

   return(NULL);
}


/************************************************************************/
/*>BOOL blDumpHash(FILE *out, HASHTABLE *hashtable)
   ------------------------------------------------
*//**
   \param[in]    out         File pointer
   \param[in]    hashtable   The hash table
   \return                   Success

   Dumps the contents of the hash to the specified file pointer

-  12.05.15  Original   By: ACRM
*/
BOOL blDumpHash(FILE *out, HASHTABLE *hashtable)
{
   char      **keyList = NULL;

   if((keyList = blGetHashKeyList(hashtable))!=NULL)
   {
      int i;
      for(i=0; keyList[i]!=NULL; i++)
      {
         BPTR value;
         int  type;
         
         fprintf(out, "%s => ", keyList[i]);
         value = blGetHashValue(hashtable, keyList[i], &type);
         switch(type)
         {
         case HASHTYPE_INT:
            fprintf(out, "%d\n", *(int *)value);
            break;
         case HASHTYPE_CHAR:
            fprintf(out, "%c\n", *(char *)value);
            break;
         case HASHTYPE_DOUBLE:
            fprintf(out, "%f\n", *(double *)value);
            break;
         case HASHTYPE_STRING:
            fprintf(out, "%s\n", (char *)value);
            break;
         case HASHTYPE_POINTER:
         default:
            fprintf(out, "%lu\n", (unsigned long)value);
            break;
         }
      }
      blFreeHashKeyList(keyList);
      keyList = NULL;
      return(TRUE);
   }

   return(FALSE);
}


/************************************************************************/
/*>BOOL blHashKeyDefined(HASHTABLE *hashtable, char *key)
   ------------------------------------------------------
*//**
   \param[in]    hashtable   The hash table
   \param[in]    key         Key for which we are searching
   \return                   Key found

   Checks whether a hash key has been defined

-  14.05.15  Original   By: ACRM
*/
BOOL blHashKeyDefined(HASHTABLE *hashtable, char *key)
{
   
   /* If this key is not already in the hash                            */
   if(lookup(hashtable, key)==NULL)
   {
      return(FALSE);
   }
   return(TRUE);
}


/************************************************************************/
/*>void blDeleteHashKey(HASHTABLE *hashtable, char *key)
   -----------------------------------------------------
*//**
   \param[in]  hashtable   The hash table
   \param[in]  key         The hash key

   Deletes a hash key

-  14.05.15  Original   By: ACRM
*/
void blDeleteHashKey(HASHTABLE *hashtable, char *key)
{
   ULONG      hashIndex;
   _HASHENTRY *hashEntry, 
              *prev, 
              *h,
              **hashtab;

   if(hashtable == NULL)
      return;

   hashtab = hashtable->table;
   
   /* If this key is in the hash then we need to delete it.             */
   if((hashEntry=lookup(hashtable, key))!=NULL)
   {
      /* If it's a string, free the memory for it                       */
      if((hashEntry->type == HASHTYPE_STRING) && 
         (hashEntry->string != NULL))
         free(hashEntry->string);

      hashIndex = hash(key, hashtable->size);
      prev=NULL;
      for(h=hashtab[hashIndex]; h!=NULL; NEXT(h))
      {
         if(h==hashEntry)
         {
            if(prev==NULL)
            {
               hashtab[hashIndex] = h->next;
            }
            else
            {
               prev->next = h->next;
            }
            break;
         }
         prev=h;
      }

      free(hashEntry->key);
      free(hashEntry);
   }
}


/************************************************************************/
/**                                                                    **/
/** Static functions below here                                        **/
/**                                                                    **/
/************************************************************************/

#ifdef DEBUG
/************************************************************************/
/*>static void dumpHashTable(FILE *out, HASHTABLE *hashtable)
   ----------------------------------------------------------
*//**
   \param[in]   out        A file pointer
   \param[in]   hashtable  The hash table
   \return

   A simple debugging function to display the hash table as
   (key => value) pairs. Prints the table itself with all entries in
   each slot

-  12.05.15  Original   By: ACRM
*/
static void dumpHashTable(FILE *out, HASHTABLE *hashtable)
{
   int        i;
   _HASHENTRY **hashtab;

   if(hashtable == NULL)
      return;
   
   hashtab = hashtable->table;

   for(i=0; i<hashtable->size; i++)
   {
      if(hashtab[i]==NULL)
      {
         fprintf(out, "()\n");
      }
      else
      {
         _HASHENTRY *hashEntry;
         
         fprintf(out, "(");
         for(hashEntry=hashtab[i]; hashEntry!=NULL; NEXT(hashEntry))
         {
            switch(hashEntry->type)
            {
            case HASHTYPE_INT:
               fprintf(out, "(%s => %d) ",
                       hashEntry->key, hashEntry->intValue);
               break;
            case HASHTYPE_STRING:
               fprintf(out, "(%s => %s) ",
                       hashEntry->key, hashEntry->string);
               break;
            case HASHTYPE_DOUBLE:
               fprintf(out, "(%s => %f) ",
                       hashEntry->key, hashEntry->doubleValue);
               break;
            case HASHTYPE_CHAR:
               fprintf(out, "(%s => %c) ",
                       hashEntry->key, hashEntry->charValue);
               break;
            case HASHTYPE_POINTER:
               fprintf(out, "(%s => %lu) ",
                       hashEntry->key, (unsigned long)hashEntry->ptrValue);
               break;
            default:
               break;
            }
         }
         
         fprintf(out, ")\n");
      }
   }
}
#endif


/************************************************************************/
/*>static ULONG hash(char *string, ULONG size)
   -------------------------------------------
*//**
   \param[in]    string    The string to be hashed
   \param[in]    size      The size of the hash table
   \return                 Index into the hash table

   Simple hashing function for strings

-  12.05.15  Original   By: ACRM
*/
static ULONG hash(char *string, ULONG size)
{
   ULONG h;

   for(h=0; *string; string++)
      h=(*string)+(h*31);

   return(h % size);
}


/************************************************************************/
/*>static _HASHENTRY *lookup(HASHTABLE *hashtable, char *key)
   --------------------------------------------------------
*//**
   \param[in]   hashtable     The hash table
   \param[in]   key           A key in the table
   \return                    A pointer to an entry in the hash

   Basic lookup function to access a _HASHENTRY structure in the hash.
   First looks up what slot we are in using the hash() function then
   searches that slot for the specified key.

-  12.05.15  Original   By: ACRM
*/
static _HASHENTRY *lookup(HASHTABLE *hashtable, char *key)
{
   ULONG       hashIndex;
   _HASHENTRY  *hashEntry,
               **hashtab;

   if(hashtable == NULL)
      return(NULL);

   hashtab   = hashtable->table;
   hashIndex = hash(key, hashtable->size);

   for(hashEntry=hashtab[hashIndex]; hashEntry!=NULL; NEXT(hashEntry))
   {
      if(!strcmp(hashEntry->key,key))
         return(hashEntry);
   }

   return(NULL);
}


/************************************************************************/
/*>static ULONG CalculateHashSize(ULONG hashsize)
   ----------------------------------------------
*//**
   \param[in]     hashsize    Basic size of hash
   \return                    Better size

   Calculates a good size for the hash - this is a prime number
   greater than the specified size.

   Currently our prime numbers only go up to 60013, so in cases
   larger than this we just return an odd number >= the specified
   size.

-  12.05.15  Original   By: ACRM
-  15.05.15  Now returns the next largest prime
*/
static ULONG CalculateHashSize(ULONG hashsize)
{
   ULONG newSize;

   if(hashsize == 0)
      return(DEF_HASHSIZE);

   if((newSize = blFindNextPrime(hashsize, FALSE)) == 0)
      return(1+(2*((int)(hashsize/2))));    /* Odd number >= hashsize   */

   return(newSize);
}

#ifdef TEST
/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Test code

-  12.05.15  Original   By: ACRM
*/
int main(int argc, char **argv)
{
   int       i;
   char      *keys[]={"name","address","phone","room101",NULL};
   char      *values[]={"Andrew","London","02071234567","Contents",NULL};
   HASHTABLE *hashtable = NULL;

   if((hashtable = blInitializeHash(0))==NULL)
   {
      fprintf(stderr,"No memory for hash table\n");
      return(1);
   }
  
   for(i=0; keys[i]!=NULL; i++)
      blSetHashValue(hashtable, keys[i], HASHTYPE_STRING, values[i]);
   
   printf("\n\nInitial hash...\n");
   blDumpHash(stdout, hashtable);

   blSetHashValue(hashtable, "phone",   HASHTYPE_STRING, "020898765432");
   blSetHashValue(hashtable, "an int",  HASHTYPE_INT, 10);
   blSetHashValue(hashtable, "double",  HASHTYPE_DOUBLE, 123.456);
   blSetHashValue(hashtable, "pointer", HASHTYPE_POINTER, &i);
   blSetHashValue(hashtable, "char",    HASHTYPE_CHAR, 'a');
   /* Change this one again                                             */
   blSetHashValue(hashtable, "an int",  HASHTYPE_INT, 99);

   blSetHashValueString(hashtable, "phone2",    "020898765432");
   blSetHashValueInt(hashtable,     "int2",     10);
   blSetHashValueDouble(hashtable,  "double2",  123.456);
   blSetHashValuePointer(hashtable, "pointer2", (BPTR)&i);
   blSetHashValueChar(hashtable,    "char2",    'a');

   printf("\n\nAfter updating...\n");
   blDumpHash(stdout, hashtable);
   
#ifdef DEBUG
   dumpHashTable(stdout, hashtable);
#endif

   printf("\n\nSome individual values:\n");
   printf("an int  : %d\n",  
          blGetHashValueInt(hashtable,     "an int"));
   printf("double  : %f\n",  
          blGetHashValueDouble(hashtable,  "double"));
   printf("pointer : %lu\n", 
          (unsigned long)blGetHashValuePointer(hashtable, "pointer"));
   printf("char    : '%c'\n",  
          blGetHashValueChar(hashtable,    "char"));
   printf("phone   : \"%s\"\n",  
          blGetHashValueString(hashtable,  "phone"));


   printf("\nChecking key existence:\n");
   printf("Key fubar: %s exist\n", 
          (blHashKeyDefined(hashtable, "fubar")?"DOES":"does NOT"));
   printf("Key phone: %s exist\n", 
          (blHashKeyDefined(hashtable, "phone")?"DOES":"does NOT"));


   printf("\nDeleting entries: \"phone\" and \"double\"\n");
   blDeleteHashKey(hashtable, "phone");
   blDeleteHashKey(hashtable, "double");
   blDumpHash(stdout, hashtable);

   blFreeHash(hashtable); 

   return(0);
}
#endif



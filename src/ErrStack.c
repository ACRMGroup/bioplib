/************************************************************************/
/**

   \file       ErrStack.c
   
   \version    V1.1
   \date       07.07.14
   \brief      Build and print an error stack for program failure.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2014
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


   This set of routines allows a stack of errors to be created. 
   When a program has a fatal error, the StoreError() routine is called
   to place the error on the stack. As the program un-winds, each
   routine which fails stores it's error. Finally, before the program
   actually exits, it calls ShowErrors() to display the error stack.

**************************************************************************

   Usage:
   ======

   blStoreError(char *routine, char *error)

   The routine is called with the name of the routine at fault and the
   description of the fault.

   blShowErrors(void *PrintRoutine, BOOL Trace)

   The routine is called with a pointer to the routine which is to
   do the actual error display and a flag to indicate whether the
   faulty routine names should be displyed. This is only of use if
   the user has access to the source code so should be used for
   debugging purposes only.

   If PrintRoutine is supplied as NULL, the simple PrintAnError()
   routine will be used which displays the error on stderr. More
   complex routines could, for example, show the error in a requester
   or output to a window.

   If a print routine is specified, the routine is called with:

   blShowErrors(void *)MyRoutine, TRUE);

**************************************************************************

   Revision History:
   =================
-  V1.0  31.08.94 Original    By: ACRM
-  V1.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    General Programming
   #SUBGROUP Handling errors
   #ROUTINE  blStoreError()
   Stores an error on the error stack.

   #ROUTINE  blShowErrors()
   Display the error stack using the supplied print routine or the
   simple default one if NULL is given.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SysDefs.h"
#include "macros.h"

#include "ErrStack.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct _errorstack
{
   struct _errorstack *next;
   char               *error,
                      *routine;
}  ERRORSTACK;

/************************************************************************/
/* Globals
*/
static ERRORSTACK *sErrorStack = NULL;

/************************************************************************/
/* Prototypes
*/
static void PrintAnError(char *error);

/************************************************************************/
/*>void blStoreError(char *routine, char *error)
   ---------------------------------------------
*//**

   \param[in]     *routine        Name of the routine generating the error
   \param[in]     *error          Description of the error

   Stores an error on the error stack.

-  31.08.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blStoreError(char *routine, char *error)
{
   static ERRORSTACK *p = NULL;
   
   if(sErrorStack == NULL)
   {
      INIT(sErrorStack, ERRORSTACK);
      p = sErrorStack;
   }
   else
   {
      ALLOCNEXT(p, ERRORSTACK);
   }
   
   if(p != NULL)
   {
      if((p->error = (char *)malloc((strlen(error)+1) * sizeof(char)))
         != NULL)
      {
         strcpy(p->error, error);
      }
      if((p->routine = (char *)malloc((strlen(routine)+1) * sizeof(char)))
         != NULL)
      {
         strcpy(p->routine, routine);
      }
   }
}

/************************************************************************/
/*>void blShowErrors(void *PrintRoutine(char *), BOOL Trace)
   ---------------------------------------------------------
*//**

   \param[in]     *PrintRoutine           The print routine or NULL
   \param[in]     Trace                   Flag to print routine names

   Display the error stack using the supplied print routine or the
   simple default one if NULL is given.

-  31.08.94 Original    By: ACRM
-  06.09.94 No longer tries to set PrintRoutine if was NULL (strict ANSI
            compliance)
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blShowErrors(void *PrintRoutine(char *), BOOL Trace)
{
   ERRORSTACK *p;
   char       buffer[160];

   for(p=sErrorStack; p!=NULL; NEXT(p))
   {
      if(Trace)
         sprintf(buffer,"%s :  %s",p->routine,p->error);
      else
         strcpy(buffer,p->error);
      
      if(PrintRoutine == NULL)
         PrintAnError(buffer);
      else
         (*PrintRoutine)(buffer);
   }
}

/************************************************************************/
/*>static void PrintAnError(char *string)
   ---------------------------------------
*//**

   \param[in]     *string        A string to be printed

   A simple error printing routine used if NULL given as a parameter
   to ShowErrors()

-  31.08.94 Original    By: ACRM
*/
static void PrintAnError(char *string)
{
   fprintf(stderr,"%s\n",string);
}


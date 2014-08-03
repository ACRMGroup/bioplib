/************************************************************************/
/**

   \file       deprecated.h
   
   \version    v1.2
   \date       31.07.14
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

   Contains the DEPRECATED() macro which prints a warning message when a 
   call is made to a deprecated function.


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   
-  V1.0  07.05.14 Original By: CTP
-  V1.1  07.07.14 Rename functions with 'bl' prefix. By: CTP
-  V1.2  31.07.14 Rewrite of deprecation system. Moved deprecated function
                  code to deprecated.c  By: CTP 

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
*//**

   The DEPRECATED macro gives a warning message if a function is 
   deprecated and indicates the replacement function.
   
   The default option is to give a warning message unless an environment 
   variable, BIOPLIB_DEPRECATED_QUIET, is set.
   
   Alternatively the compile options: -D BIOPLIB_DEPRECATED_CHECK or 
   -D BIOPLIB_DEPRECATED_QUIET will set the DEPRECATED macro to ignore the
   BIOPLIB_DEPRECATED_QUIET environment variable. 
   
   -D BIOPLIB_DEPRECATED_CHECK will display the warning message. 
   -D BIOPLIB_DEPRECATED_QUIET will always silence the warning message.
   
-  29.04.14 Original    By: CTP
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
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/

#endif

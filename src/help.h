/************************************************************************/
/**

   \file       help.h
   
   \version    V1.2
   \date       31.07.14
   \brief      Include file for help functions
   
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
-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP

*************************************************************************/
#ifndef _HELP_H
#define _HELP_H

/*#include "deprecated.h"*/

/* Prototypes */
void blHelp(char *string, char *HelpFile);
void blDoHelp(char *string, char *HelpFile);

/************************************************************************/
/* Deprecated functions: help.h                                         */
/** \cond deprecated                                                    */

void Help(char *string, char *HelpFile);
void DoHelp(char *string, char *HelpFile);

/* \endcond                                                             */
/************************************************************************/

#endif

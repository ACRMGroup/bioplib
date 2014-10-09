/************************************************************************/
/**

   \file       parseresspec.h
   
   \version    V1.0
   \date       05.08.14
   \brief      Include file for ParseResSpec test suite.
   
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

   Test suite for blParseResSpec().

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.08.14 Original By: CTP

*************************************************************************/
#ifndef _PARSERESSPEC_H
#define _PARSERESSPEC_H

/* Includes for tests */
#include <stdlib.h>
#include <check.h>

/* Includes from source file */
#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "../../macros.h"
#include "../../SysDefs.h"
#include "../../pdb.h"

/* Prototypes */
Suite *parseresspec_suite(void);

#endif

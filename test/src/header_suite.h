/************************************************************************/
/**

   \file       header_suite.h
   
   \version    V1.0
   \date       05.05.15
   \brief      Include file for CONECT test suite.
   
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

   Test suite for reading and writing Header data for pdbml files.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.05.15 Original By: CTP

*************************************************************************/

#ifndef _HEADER_H
#define _HEADER_H

/* Includes for tests */
#include <stdlib.h>
#include <check.h>
#include <unistd.h>

/* Includes from source file */
#include "port.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "macros.h"
#include "general.h"
#include "pdb.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <time.h>

/* Prototypes */
Suite *header_suite(void);

#endif

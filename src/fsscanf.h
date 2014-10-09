/************************************************************************/
/**

   \file       fsscanf.h
   
   \version    V1.1R
   \date       01.03.94
   \brief      Include file for fsscanf()
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994
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

*************************************************************************/
/* Includes
*/

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

#ifndef _FSSCANF_H
#define _FSSCANF_H

int fsscanf(char *buffer, char *format, ...);

#endif

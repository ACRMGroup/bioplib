/************************************************************************/
/**

   \file       WindIO.h
   
   \version    V1.3R
   \date       18.10.95
   \brief      Header for window/normal interface routines

   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-5
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
#ifndef _WINDIO_H
#define _WINDIO_H

#include "SysDefs.h"

void screen(char *string);
void prompt(char *string);
void RePrompt(void);
void GetKybdString(char *string, int maxlen);
void PagingOn(void);
void PagingOff(void);
void WindowMode(BOOL mode);
void WindowInteractive(BOOL mode);
int YorN(char deflt);

#endif

/************************************************************************/
/**

   \file       SysDefs.h
   
   \version    V1.3
   \date       14.05.15
   \brief      System-type variable type definitions
   
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


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  01.03.94 Original    By: ACRM
-  V1.1  02.08.95 Added UCHAR
-  V1.2  01.02.96 Added UBYTE
-  V1.3  14.05.15 Added BPTR

*************************************************************************/
#ifndef _SYSDEFS_H
#define _SYSDEFS_H

#ifndef EXEC_TYPES_H    /* Commodore Amiga; defines in <exec/types.h>   */
typedef void            *APTR;
typedef char            *BPTR;

#ifndef SYS_TYPES_H     /* Unix: <sys/types.h>, MS-DOS: <sys\types.h>   */
#ifndef _TYPES_         /* Ditto                                        */
typedef short           BOOL;
typedef long            LONG;
typedef unsigned long   ULONG;
typedef short           SHORT;
typedef unsigned short  USHORT;
typedef unsigned char   UCHAR;
typedef unsigned char   UBYTE;
#endif
#endif
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifdef _ESV_
typedef long            time_t;   /* Required on E&S System V           */
typedef long            clock_t;  /* Ditto                              */
#define CLOCKS_PER_SEC  1000000   /* Ditto                              */
#endif

#endif

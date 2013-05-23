/*************************************************************************

   Program:    
   File:       port.h
   
   Version:    V1.0R
   Date:       27.02.98
   Function:   Port-specific defines to allow us to use things like 
               popen() in a clean compile
   
   Copyright:  (c) SciTech Software 1988-1998
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======
   This must be included before the system includes (i.e. stdio.h etc)

**************************************************************************

   Revision History:
   =================
   V1.0  27.02.98  Original

*************************************************************************/
/***
 *** The following are necessary for getting popen() to work cleanly      
 ***/

/* Silicon graphics IRIX                                                */
#ifdef __sgi     
/*  These are for Irix6                                                 */
#   ifndef __EXTENSIONS__
#      define __EXTENSIONS__
#   endif

#   ifdef _POSIX_C_SOURCE
#      undef _POSIX_C_SOURCE
#   endif
#   define _POSIX_C_SOURCE 2

/*  These are for Irix5                                                 */
#   ifdef _POSIX_SOURCE
#      undef _POSIX_SOURCE
#   endif
#endif

/* DEC OSF/1 (Alpha)                                                    */
#ifdef __osf__
#   ifndef _XOPEN_SOURCE
#      define _XOPEN_SOURCE
#   endif
#   ifndef __STDC__
#      define __STDC__
#   endif
#endif

/* SunOS - doesn't need anythging for SunOS 4.1.2                       */
#if defined(sun) && defined(sparc) && !defined(__svr4)
#endif

/* Linux                                                                */
#ifdef linux
#   define _POSIX_SOURCE
#endif


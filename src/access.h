/************************************************************************/
/**

   \file       access.h
   
   \version    V1.2
   \date       17.06.16
   \brief      Accessibility calculation code
   
   \copyright  (c) UCL, Dr. Andrew C.R. Martin, 1999-2015
   \author     Dr. Andrew C.R. Martin
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
-  V1.0  21.04.99 Original   By: ACRM
-  V1.1  17.07.14 Extracted from XMAS code
-  V1.2  17.06.15 Added scAccess and scRelAccess to RESACCESS structure

*************************************************************************/
#ifndef _ACCESS_H_
#define _ACCESS_H_ 1

#define ACCESS_MAX_ATOMS_PER_RESIDUE 50
#define ACCESS_DEF_INTACC            0.05

#ifndef VERY_SMALL
#define VERY_SMALL            (REAL)1e-6
#endif

/* Used to store the atom radii and standard accessibility data         */
typedef struct _resrad
{
   struct _resrad *next;
   REAL  stdaccess,
         radius[ACCESS_MAX_ATOMS_PER_RESIDUE];
   int   natoms;
   char  resnam[8],
         atnam[ACCESS_MAX_ATOMS_PER_RESIDUE][8];
}  RESRAD;

/* Used to store residue accessibility values                           */
typedef struct _resaccess
{
   struct _resaccess *next;
   REAL resAccess,
        relAccess,
        scAccess,
        scRelAccess;
   int  resnum;
   char chain[8],
        resnam[8],
        insert[8];
}  RESACCESS;

/* Prototypes                                                           */
RESRAD *blSetAtomRadii(PDB *pdb, FILE *fpRad);
BOOL blCalcAccess(PDB *pdb, int natoms, 
                  REAL integrationAccuracy, REAL probeRadius,
                  BOOL doAccessibility);
RESACCESS *blCalcResAccess(PDB *pdb, RESRAD *resrad);

#endif



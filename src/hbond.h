/************************************************************************/
/**

   \file       hbond.h
   
   \version    V1.0
   \date       26.01.96
   \brief      Header file for hbond determining code
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1996
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
-  V1.0  26.01.96 Original    By: ACRM

*************************************************************************/
#ifndef _hbond_h
#define _hbond_h

/************************************************************************/
/* Defines and macros
*/
#define HBOND_BACK1 0x01
#define HBOND_BACK2 0x02
#define HBOND_SIDE1 0x04
#define HBOND_SIDE2 0x08

#define HBOND_BB        (HBOND_BACK1 | HBOND_BACK2)
#define HBOND_BS        (HBOND_BACK1 | HBOND_SIDE2)
#define HBOND_SS        (HBOND_SIDE1 | HBOND_SIDE2)
#define HBOND_SB        (HBOND_SIDE1 | HBOND_BACK2)
#define HBOND_SIDECHAIN (HBOND_SIDE1 | HBOND_SIDE2 | HBOND_BACK2)
#define HBOND_BACKBONE  (HBOND_BACK1 | HBOND_SIDE2 | HBOND_BACK2)

#define HBOND_ANY (HBOND_BACK1 | HBOND_BACK2 | HBOND_SIDE1 | HBOND_SIDE2)

/************************************************************************/
/* Prototypes
*/
int  IsHBonded(PDB *res1, PDB *res2, int type);
BOOL ValidHBond(PDB *AtomH, PDB *AtomD, PDB *AtomA, PDB *AtomP);
int IsMCDonorHBonded(PDB *res1, PDB *res2, int type);
int IsMCAcceptorHBonded(PDB *res1, PDB *res2, int type);


#endif

/************************************************************************/
/**

   \file       cssr.h
   
   \version    V1.2
   \date       31.07.14
   \brief      Defines for CSSR handling
   
   \copyright  SciTech Software 1991-2014
   \author     Andrew C. R. Martin
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

   Structure definitions for ReadCSSR()

****************************************************************************

   Usage:
   ======

****************************************************************************

   Revision History:
   =================
-  V1.1  07.07.14 Use bl prefix for functions By: CTP
-  V1.2  31.07.14 Updated deprecation: Removed deprecated.h and added 
                  prototypes for renamed functions. By: CTP

***************************************************************************/
#ifndef _CSSR_H
#define _CSSR_H

#include "MathType.h"
#include "pdb.h"

struct cssr_entry
{
   REAL  charge;
   REAL  x,y,z;
   struct cssr_entry *next;
   int   atnum;
   int   group;
   int   link[8];
   char  atnam[8];
} ;

typedef struct cssr_entry CSSR;

#define CLEAR_CSSR(p) p->atnum=0; \
                      strcpy(p->atnam,"    "); \
                      p->group=0; \
                      p->x = 0.0; p->y = 0.0; p->z = 0.0; \
                      p->charge = 0.0; \
                      p->link[0] = p->link[1] = p->link[2] = p->link[3] = 0; \
                      p->link[4] = p->link[5] = p->link[6] = p->link[7] = 0; \
                      p->next = NULL

/************************************************************************/
/* Prototypes
*/
CSSR *blReadCSSR(FILE *fp, int *natom, char *name, char *title);
PDB *blReadCSSRasPDB(FILE *fp, int *natom);
void blNormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta,
                     REAL gamma);
void blNormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta,
                    REAL gamma);
void blOrtho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
           REAL amatrx[3][3], int isw, int ncode);
/* void blPadterm(char *string, int len);*/ /* defined in general.h */
void blWriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title);

/************************************************************************/
/* Deprecated functions: cssr.h                                         */
/** \cond deprecated                                                    */

CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title);
PDB *ReadCSSRasPDB(FILE *fp, int *natom);
void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta,
                   REAL gamma);
void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta,
                  REAL gamma);
void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
           REAL amatrx[3][3], int isw, int ncode);
void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title);

/* \endcond                                                             */
/************************************************************************/

#endif

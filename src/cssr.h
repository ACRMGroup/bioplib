/************************************************************************/
/**

   \file       cssr.h
   
   \version    V1.0R
   \date       09.09.91
   \brief      Defines for CSSR handling
   
   \copyright  SciTech Software 1991
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
CSSR *ReadCSSR(FILE *fp, int *natom, char *name, char *title);
PDB *ReadCSSRasPDB(FILE *fp, int *natom);
void NormaliseCSSR(CSSR *cssr, REAL cell[3], REAL alpha, REAL beta,
                   REAL gamma);
void NormalisePDB(PDB *pdb, REAL cell[3], REAL alpha, REAL beta,
                  REAL gamma);
void ortho(REAL cell[3], REAL alpha, REAL beta, REAL gamma,
           REAL amatrx[3][3], int isw, int ncode);
void padterm(char *string, int len);
void WriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title);

#endif

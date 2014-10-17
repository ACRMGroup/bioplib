/************************************************************************/
/**

   \file       WrtCSSR.c
   
   \version    V1.4
   \date       07.07.14
   \brief      Write a CSSR file
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1991-2014
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

   blWriteCSSR(fp,cssr,name,title)

   This subroutine will write a CSSR file from a linked list of structures 
   of type cssr_entry. 

   The strucure is set up by including the file "cssr.h". For details of 
   the structure, see this file.

**************************************************************************

   Usage:
   ======

\code
   WriteCSSR(fp,cssr,name,title)
\endcode

   \param[in]     *fp      A pointer to type FILE in which the
                              CSSR file is stored.
   \param[in]     *cssr    A pointer to the first allocated item of
                              the CSSR linked list
   \param[in]     *name    The molecule's name.
   \param[in]     *title   Title on the molecule.

**************************************************************************

   Revision History:
   =================
-  V1.0  22.09.91 Original
-  V1.1  01.06.92 Autodoc'd. Added FPU check.
-  V1.2  10.06.93 void return; float->REAL
-  V1.3  27.07.93 %f -> %lf
-  V1.3  01.03.94 %lf -> %f  (!)
-  V1.4  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling CSSR Data
   #SUBGROUP File IO
   #FUNCTION  blWriteCSSR()
   Write a CSSR file from a CSSR linked list.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "MathType.h"

#include "cssr.h"             /* Note cssr includes pdb                 */
#include "macros.h"

/************************************************************************/
/* Defines and macros
*/
#define CR 13
#define LF 10

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void blWriteCSSR(FILE *fp, CSSR *cssr, char *name, char *title)
   ---------------------------------------------------------------
*//**

   \param[in]     *fp      A pointer to type FILE in which the
                           CSSR file is stored.
   \param[in]     *cssr    A pointer to the first allocated item of
                           the CSSR linked list
   \param[in]     *name    The molecule's name.
   \param[in]     *title   Title on the molecule.

   Write a CSSR file from a CSSR linked list.

-  22.09.91 Original
-  01.06.92 Autodoc'd
-  10.06.93 void return; float->REAL
-  27.07.93 %f -> %lf
-  01.03.94 %lf -> %f (!)
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blWriteCSSR(FILE  *fp,
                 CSSR  *cssr,
                 char  *name,
                 char  *title)
{
   REAL  cell[3],
         alpha, beta, gamma;
   int   natom;
   CSSR  *p;
   
   /* Initilialise cell parameters                                      */
   cell[0] = cell[1] = cell[2] = 1.0;
   alpha = beta = gamma = 90.0;
   
   /*** Record 1, assume no REFCODE                                   ***/
   fprintf(fp,"                                      %8.3f%8.3f%8.3f\n",
           cell[0],cell[1],cell[2]);
   
   /*** Record 2                                                      ***/
   fprintf(fp,"                     %8.3f%8.3f%8.3f\n",
           alpha,beta,gamma);

   /*** Record 3                                                      ***/
   /* Assume orthogonal.
      Count the atoms   
   */
   for(p=cssr, natom=0; p; NEXT(p)) natom++;
   fprintf(fp,"%4d   1 %-60s\n",natom,name);
   
   /*** Record 4                                                      ***/
   /* We'll assume the charges are valid (even if 0.0)                  */
   fprintf(fp,"     2 %-60s\n",title);
   
   /*** Remaining records represent the atoms...                      ***/
   for(p=cssr; p; NEXT(p))
   {
      fprintf(fp,"%4d %-4s  %9.5f %9.5f %9.5f \
%4d%4d%4d%4d%4d%4d%4d%4d %8.4f %2d\n",
              p->atnum,
              p->atnam,
              p->x,
              p->y,
              p->z,
              p->link[0],
              p->link[1],
              p->link[2],
              p->link[3],
              p->link[4],
              p->link[5],
              p->link[6],
              p->link[7],
              p->charge,
              p->group);
   }
}


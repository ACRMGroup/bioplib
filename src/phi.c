/************************************************************************/
/**

   \file       phi.c
   
   \version    V1.8
   \date       17.07.14
   \brief      Calculate a torsion angle
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2014
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

-  V1.7  07.07.14 Use bl prefix for functions By: CTP
-  V1.8  17.07.14 Removed unused varables  By: ACRM

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Maths
   #SUBGROUP Geometry
   #FUNCTION  blPhi()
   Calculates the torsion angle described by 4 sets of coordinates.
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"

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
/*>REAL blPhi(REAL xi,REAL yi,REAL zi,REAL xj,REAL yj,REAL zj,
              REAL xk,REAL yk,REAL zk,REAL xl,REAL yl,REAL zl)
   ---------------------------------------------------------
*//**

   \param[in]     xi          Input coordinate
   \param[in]     yi          Input coordinate
   \param[in]     zi          Input coordinate
   \param[in]     xj          Input coordinate
   \param[in]     yj          Input coordinate
   \param[in]     zj          Input coordinate
   \param[in]     xk          Input coordinate
   \param[in]     yk          Input coordinate
   \param[in]     zk          Input coordinate
   \param[in]     xl          Input coordinate
   \param[in]     yl          Input coordinate
   \param[in]     zl          Input coordinate
   \return                    The torsion angle between the 4 atoms

   Calculates the torsion angle described by 4 sets of coordinates.

-  04.03.91 Original    By: ACRM
-  16.06.93 Changed float to REAL
-  07.07.14 Use bl prefix for functions By: CTP
-  17.07.14 Removed unused variables  By: ACRM
*/
REAL blPhi(REAL xi,
           REAL yi,
           REAL zi,
           REAL xj,
           REAL yj,
           REAL zj,
           REAL xk,
           REAL yk,
           REAL zk,
           REAL xl,
           REAL yl,
           REAL zl)
{
   REAL xij,yij,zij,
        xkj,ykj,zkj,
        xkl,ykl,zkl,
        dxi,dyi,dzi,
        gxi,gyi,gzi,
        bi,bk,ct,
        z1,z2,ap,s;
   
   /* Calculate the vectors C,B,C                                       */
   xij = xi - xj;
   yij = yi - yj;
   zij = zi - zj;
   xkj = xk - xj;
   ykj = yk - yj;
   zkj = zk - zj;
   xkl = xk - xl;
   ykl = yk - yl;
   zkl = zk - zl;

   /* Calculate the normals to the two planes n1 and n2
      this is given as the cross products:
       AB x BC
      --------- = n1
      |AB x BC|

       BC x CD
      --------- = n2
      |BC x CD|
   */
   dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
   dyi = zij * xkj - xij * zkj;
   dzi = xij * ykj - yij * xkj;
   gxi = zkj * ykl - ykj * zkl;     /* Normal to plane 2                */
   gyi = xkj * zkl - zkj * xkl;
   gzi = ykj * xkl - xkj * ykl;

   /* Calculate the length of the two normals                           */
   bi = dxi * dxi + dyi * dyi + dzi * dzi;
   bk = gxi * gxi + gyi * gyi + gzi * gzi;
   ct = dxi * gxi + dyi * gyi + dzi * gzi;

   bi   = (REAL)sqrt((double)bi);
   bk   = (REAL)sqrt((double)bk);

   z1   = 1./bi;
   z2   = 1./bk;

   ct   = ct * z1 * z2;
   if (ct >  1.0)   ct = 1.0;
   if (ct < (-1.0)) ct = -1.0;
   ap   = acos(ct);

   s = xkj * (dzi * gyi - dyi * gzi)
     + ykj * (dxi * gzi - dzi * gxi)
     + zkj * (dyi * gxi - dxi * gyi);

   if (s < 0.0) ap = -ap;

   ap = (ap > 0.0) ? PI-ap : -(PI+ap);

   return(ap);
}


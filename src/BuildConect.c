/************************************************************************/
/**

   \file       BuildConect.c
   
   \version    V1.4
   \date       22.07.15
   \brief      Build connectivity information in PDB linked list
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2002-2015
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
-  V1.0  19.02.15 Original
-  V1.1  16.03.15 Added blDeleteAConect() and blDeleteAConectByNum()
-  V1.2  12.05.15 Fixed to add conects between C and N if either is
                  a HETATM
-  V1.3  23.06.15 Added blAreResiduesBonded() and 
                  blAreResiduePointersBonded()
-  V1.4  22.07.15 Added blIsConnected()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Miscellaneous functions

   #FUNCTION blBuildConectData()
   Rebuild all CONECT data using covalent radii of the atoms

   #FUNCTION blAddConect()
   Adds a CONECT in both directions between two specified atoms

   #FUNCTION blAddOneDirectionConect()
   Adds a CONECT in one direction between two specified atoms

   #FUNCTION blDeleteAConect()
   Deletes a CONECT between two specified atoms

   #FUNCTION blDeleteAConectByNum()
   Deletes a specified CONECT from an atom

   #FUNCTION blDeleteAtomConects()
   Deletes all CONECT information associated with a PDB pointer. Also
   deletes the relevant CONECTs (back to this atom) from the partner 
   atoms

   #FUNCTION blCopyConects()
   Updates the CONECT information in the out PDB linked list based on
   those in the in linked list. This is used when we have copied a PDB
   linked list (e.g. in the bl...AsCopy() routines) to make sure the
   CONECT data points to records in the new linked list instead of the
   old one.

   #FUNCTION blAreResiduesBonded()
   Tests whether residues specified by chain, number and insert are
   bonded

   #FUNCTION blAreResiduePointersBonded()
   Tests whether residues specified by pointers to the start of the
   residues

   #FUNCTION blIsBonded()
   Test whether two atoms are bonded

   #FUNCTION blIsConected()
   Tests whether two atoms are listed as connected in the CONECT records
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include "macros.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
*/
struct _covalentradii
{
   char element[8];
   REAL radius;
};
   
/************************************************************************/
/* Globals
*/
static struct _covalentradii covalentRadii[] = 
{
   {"C",  0.76}, {"O",  0.66}, {"N",  0.71}, {"S",  1.05}, {"P",  1.11},
   {"H",  0.31}, {"CL", 1.02}, {"K",  2.03}, {"NA", 1.66}, {"MN", 1.61},
   {"CO", 1.50}, {"NI", 1.24}, {"CU", 1.32}, {"ZN", 1.22}, {"SE", 1.20},
   {"BR", 1.20}, {"MG", 1.41}, {"CA", 1.76}, {"FE", 1.52}, {"LI", 1.33},
   {"BE", 1.02}, {"B",  0.85}, {"F",  0.64}, {"NE", 0.67}, {"AL", 1.26},
   {"SI", 1.16}, {"AR", 1.07}, {"SC", 1.70}, {"TI", 1.60}, {"V",  1.53},
   {"CR", 1.39}, {"GA", 1.24}, {"GE", 1.21}, {"AS", 1.21}, {"KR", 1.21},
   {"RB", 2.20}, {"SR", 1.95}, {"Y",  1.90}, {"ZR", 1.75}, {"NB", 1.64},
   {"MO", 1.54}, {"TC", 1.47}, {"RU", 1.46}, {"RH", 1.42}, {"PD", 1.39},
   {"AG", 1.45}, {"CD", 1.44}, {"IN", 1.46}, {"SN", 1.40}, {"SB", 1.40},
   {"TE", 1.38}, {"I",  1.39}, {"XE", 1.40}, {"CS", 2.44}, {"BA", 2.15},
   {"LU", 1.75}, {"HF", 1.70}, {"TA", 1.62}, {"W",  1.51}, {"RE", 1.44},
   {"OS", 1.41}, {"IR", 1.36}, {"PT", 1.36}, {"AU", 1.32}, {"HG", 1.45},
   {"TL", 1.46}, {"PB", 1.48}, {"BI", 1.51}, {"PO", 1.50}, {"AT", 1.50},
   {"RN", 1.50}, {"FR", 2.60}, {"RA", 2.21}, {"LR", 1.61}, {"RF", 1.57},
   {"DB", 1.49}, {"SG", 1.43}, {"BH", 1.41}, {"HS", 1.34}, {"MT", 1.29},
   {"DS", 1.28}, {"RG", 1.21}, {"CN", 1.37}, {"FL", 1.43}, {"LV", 1.75},
   {"LA", 2.07}, {"CE", 2.04}, {"PR", 2.03}, {"ND", 2.01}, {"PM", 1.99},
   {"SM", 1.98}, {"EU", 1.98}, {"GD", 1.96}, {"TB", 1.94}, {"DY", 1.92},
   {"HO", 1.92}, {"ER", 1.89}, {"TM", 1.90}, {"YB", 1.87}, {"AC", 2.15},
   {"TH", 2.06}, {"PA", 2.00}, {"U",  1.96}, {"NP", 1.90}, {"PU", 1.87},
   {"AM", 1.80}, {"CM", 1.69}, {"BK", 1.68}, {"CF", 1.68}, {"ES", 1.65},
   {"FM", 1.67}, {"MD", 1.73}, {"NO", 1.76}, {"HE", 0.46}, {"\0", 0.00}
};


/************************************************************************/
/* Prototypes
*/
static REAL findCovalentRadius(char *element);


/************************************************************************/
/*>BOOL blAddConect(PDB *p, PDB *q)
   --------------------------------
*//**

   \param[in,out]   *p    First PDB item
   \param[in,out]   *q    Second PDB item
   \return                Success?

   Adds a conect bteween p and q (i.e. in both directions)
   Fails if there are too many CONECTs 

-  19.02.15  Original   By: ACRM
*/
BOOL blAddConect(PDB *p, PDB *q)
{
   BOOL retval = TRUE;

   if(!blAddOneDirectionConect(p,q))
      retval=FALSE;
   if(!blAddOneDirectionConect(q,p))
      retval=FALSE;

   return(retval);
}


/************************************************************************/
/*>BOOL blAddOneDirectionConect(PDB *p, PDB *q)
   --------------------------------------------
*//**

   \param[in,out]   *p    First PDB item
   \param[in,out]   *q    Second PDB item
   \return                Success?

   Adds a conect from p to q (i.e. one direction only)     
   Fails if there are too many CONECTs 

-  19.02.15  Original   By: ACRM
*/
BOOL blAddOneDirectionConect(PDB *p, PDB *q)
{
   int i;
   BOOL gotConect;
   
   /* See if we have the conect already                                 */
   gotConect=FALSE;
   for(i=0; i<p->nConect; i++)
   {
      if(p->conect[i] == q)
      {
         gotConect = TRUE;
         break;
      }
   }

   /* If we haven't got it already store it                             */
   if(!gotConect)
   {
      if(p->nConect < MAXCONECT)
      {
         p->conect[p->nConect] = q;
         (p->nConect)++;
      }
      else
      {
         return(FALSE);
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>BOOL blBuildConectData(PDB *pdb, REAL tol)
   ------------------------------------------
*//**
   \param[in,out]   *pdb   PDB linked list
   \param[in]       tol    Tolerence for distance between atoms
   \return                 Were all CONECTs added OK

   Deletes all current connectivity data and rebuilds it using covalent
   radii data. A return of FALSE indicates that there were too many
   connections for an atom. If this happens, MAXCONECT needs to be 
   increased in pdb.h

-  19.02.15  Original   By: ACRM
-  26.02.15  Added tol paramater
-  12.05.15  Conects are built involving backbone C and N if either atom
             is a HETATM
*/
BOOL blBuildConectData(PDB *pdb, REAL tol)
{
   PDB  *p, 
        *q,
        *res,
        *nextRes;
   BOOL retval=TRUE;

   /* Clear all current connect data                                    */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->nConect = 0;
   }

   for(res=pdb; res!=NULL; res=nextRes)
   {
      nextRes = blFindNextResidue(res);

      /* Check for any HETATM connections within this residue           */
      for(p=res; p!=nextRes; NEXT(p))
      {
         for(q=p->next; q!=nextRes; NEXT(q))
         {
            if(!strncmp(p->record_type, "HETATM", 6) ||
               !strncmp(q->record_type, "HETATM", 6))
            {
               if(blIsBonded(p, q, tol))
               {
                  if(!blAddConect(p,q))
                     retval=FALSE;
               }
            }
         }
      }

      /* Check for connections between residues which don't involve backbone
         C or N, or do involve HETATMS
      */
      for(p=res; p!=nextRes; NEXT(p))
      {
         for(q=nextRes; q!=NULL; NEXT(q))
         {
            if(strncmp(p->atnam, "C   ", 4) ||
               strncmp(q->atnam, "N   ", 4) ||
               !strncmp(p->record_type, "HETATM", 6) ||
               !strncmp(q->record_type, "HETATM", 6))
            {
               if(blIsBonded(p, q, tol))
               {
                  if(!blAddConect(p,q))
                     retval=FALSE;
               }
            }
         }
      }
   }

   return(retval);
}


/************************************************************************/
/*>BOOL blIsBonded(PDB *p, PDB *q, REAL tol)
   -----------------------------------------
*//**
   \param[in]    *p   First PDB atom
   \param[in]    *q   Second PDB atom
   \param[in]    tol  Telerance for separation between atoms
   \return            Bonded?

   Test whether two atoms are bonded

-  19.02.15  Original   By: ACRM
-  26.02.15  Added tol parameter. Changed to used squared distances
*/
BOOL blIsBonded(PDB *p, PDB *q, REAL tol)
{
   REAL r1, r2, bondDist;
   r1 = findCovalentRadius(p->element);
   r2 = findCovalentRadius(q->element);
   bondDist = (r1+r2+tol);

   if(DISTSQ(p,q) <= bondDist*bondDist)
      return(TRUE);
   return(FALSE);
}


/************************************************************************/
/*>BOOL blAreResiduesBonded(PDB *pdb, 
                            char *chain1, int resnum1, char *insert1,
                            char *chain2, int resnum2, char *insert2,
                            REAL tol)
   ------------------------------------------------------------------
*//**
   \param[in] *pdb       PDB linked list
   \param[in] *chain1    First chain label
   \param[in] resnum1    First residue number
   \param[in] *insert1   First insert code
   \param[in] *chain2    Second chain label
   \param[in] resnum2    Second residue number
   \param[in] *insert2   Second insert code
   \param[in] tol        Tolerance for distances
   \return               Are they bonded

   Tests whether two residue are bonded

-  23.06.15 Original   By: ACRM
*/
BOOL blAreResiduesBonded(PDB *pdb, 
                         char *chain1, int resnum1, char *insert1,
                         char *chain2, int resnum2, char *insert2,
                         REAL tol)
{
   PDB *res1, 
       *res2;

   /* Find the residues                                                 */
   if((res1 = blFindResidue(pdb, chain1, resnum1, insert1))!=NULL)
   {
      if((res2 = blFindResidue(pdb, chain2, resnum2, insert2))!=NULL)
      {
         return(blAreResiduePointersBonded(res1, res2, tol));
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>BOOL blAreResiduePointersBonded(PDB *res1, PDB *res2, REAL tol)
   ---------------------------------------------------------------
*//**
   \param[in] *res1      Start of first residue
   \param[in] *res2      Start of second residue
   \param[in] tol        Tolerance for distances
   \return               Are they bonded

   Tests whether two residue are bonded using pointers to start of
   residues

-  23.06.15 Original   By: ACRM
*/
BOOL blAreResiduePointersBonded(PDB *res1, PDB *res2, REAL tol)
{
   PDB *res1next,
       *res2next,
       *p,
       *q;

   if((res1 != NULL) && (res2 != NULL))
   {
      /* Find the following residues                                    */
      res1next = blFindNextResidue(res1);
      res2next = blFindNextResidue(res2);

      /* Step through the atoms in each residue and see if they are
         bonded
      */
      for(p=res1; p!=res1next; NEXT(p))
      {
         for(q=res2; q!=res2next; NEXT(q))
         {
            if(blIsBonded(p, q, tol))
               return(TRUE);
         }
      }
   }
   return(FALSE);
}

/************************************************************************/
/*>static REAL findCovalentRadius(char *element)
   ---------------------------------------------
*//**
   \param[in]   *element   The element type
   \return                 The covalent bonding radius of the atom

-  19.02.15  Original   By: ACRM
*/
static REAL findCovalentRadius(char *element)
{
   int i;
   
   for(i=0; covalentRadii[i].element[0] != '\0'; i++)
   {
      if(!strcmp(element, covalentRadii[i].element))
      {
         return(covalentRadii[i].radius);
      }
   }
   return((REAL)1.0);
}


/************************************************************************/
/*>BOOL blDeleteAConect(PDB *p, PDB *q)
   -----------------------------------
*//**
   \param[in]    *p     First PDB pointer
   \param[in]    *q     Second PDB pointer
   \return              Success

   Deletes the CONECT information between the two specified atoms

-  16.03.15  Original   By: ACRM
*/
BOOL blDeleteAConect(PDB *p, PDB *q)
{
   int  cNum;
   BOOL retval = TRUE;
   
   for(cNum=0; cNum < p->nConect; cNum++)
   {
      if(p->conect[cNum] == q)
      {
         if(!blDeleteAConectByNum(p, cNum))
            retval = FALSE;
         break;
      }
   }

   for(cNum=0; cNum < q->nConect; cNum++)
   {
      if(q->conect[cNum] == p)
      {
         if(!blDeleteAConectByNum(q, cNum))
            retval = FALSE;
         break;
      }
   }

   return(retval);
}

/************************************************************************/
/*>BOOL blDeleteAConectByNum(PDB *pdb, int cNum)
   ------------------------------------------
*//**
   \param[in]    *pdb   PDB pointer
   \param[in]    cNum   Index into the ->conect[] array for the CONECT
                        to be deleted
   \return              Success

   Deletes the link for the specified CONECT. Other CONECTs are shuffled
   down.

-  16.03.15  Original   By: ACRM
*/
BOOL blDeleteAConectByNum(PDB *pdb, int cNum)
{
   int i;
   
   /* Check that the CONECT exists                                      */
   if((cNum >= pdb->nConect) || (pdb->nConect == 0))
      return(FALSE);
   
   /* Shuffle the CONECTs down                                          */
   for(i=cNum; i<pdb->nConect; i++)
   {
      PDB *next = NULL;
      if((i+1) < pdb->nConect)
         next = pdb->conect[i+1];

      pdb->conect[i] = next;
   }

   /* Decrement the number of CONECTs and return                        */
   pdb->nConect--;
   return(TRUE);
}


/************************************************************************/
/*>void blDeleteAtomConects(PDB *pdb)
   ----------------------------------
*//**
   \param[in]     *pdb     PDB pointer

   Deletes all CONECT information associated with a PDB pointer. Also
   deletes the relevant CONECTs (back to this atom) from the partner 
   atoms

-  17.03.15  Original   By: ACRM
*/
void blDeleteAtomConects(PDB *pdb)
{
   int i;

   if(pdb!=NULL)
   {
      /* For each CONECT (if there are any)                             */
      for(i=0; i<pdb->nConect; i++)
      {
         /* Find the partner                                               */
         PDB *conect = pdb->conect[i];
         
         if(conect!=NULL)
         {
            blDeleteAConect(pdb, conect);
         }
         pdb->conect[i] = NULL;
      }
      pdb->nConect = 0;
   }
}

/************************************************************************/
/*>BOOL blCopyConects(PDB *out, PDB *in)
   -------------------------------------
*//**
   \param[out]   *out   PDB linked list containing new CONECT records
   \param[in]    *in    PDB linked list containing original CONECT records
   \return              Success

   Updates the CONECT information in the out PDB linked list based on
   those in the in linked list. This is used when we have copied a PDB
   linked list (e.g. in the bl...AsCopy() routines) to make sure the
   CONECT data points to records in the new linked list instead of the
   old one.

-  17.04.15  Original   By: ACRM
*/
BOOL blCopyConects(PDB *out, PDB *in)
{
   PDB **idxOut = NULL,
       *p;
   int indexSize = 0;
   
   /* Create an indices by atom number                                  */
   if((idxOut = blIndexAtomNumbersPDB(out, &indexSize))==NULL)
      return(FALSE);
   
   /* Step through the output linked list                               */
   for(p=out; p!=NULL; NEXT(p))
   {
      int i;
      
      /* Step through the CONECT data for this atom - which will be
         from the OLD linked list
      */
      for(i=0; i<p->nConect; i++)
      {
         if(p->conect[i] != NULL)
         {
            int atomNum;
            PDB *newPtr = NULL;
         
            /* Find the atom number from old linked list                */
            atomNum = (p->conect[i])->atnum;
            
            /* Look up the pointer for this in the new linked list index*/
            if(atomNum < indexSize)
               newPtr = idxOut[atomNum];

            /* Reset pointer                                            */
            p->conect[i] = newPtr;
         }
      }
   }

   free(idxOut);
   return(TRUE);
}

/************************************************************************/
/*>BOOL blIsConected(PDB *p, PDB *q)
   ---------------------------------
*//**
   \param[in]       *p           First atom
   \param[in]       *q           Second atom
   \return                       Connected?

   Tests whether there is a link between the specified atoms in the
   CONECT list

-  07.06.99 Original   By: ACRM
*/
BOOL blIsConected(PDB *p, PDB *q)
{  int i;
   for(i=0; i<p->nConect; i++)
   {
      if(p->conect[i] == q)
         return(TRUE);
   }

   for(i=0; i<q->nConect; i++)
   {
      if(q->conect[i] == p)
         return(TRUE);
   }
   
   return(FALSE);
}


/************************************************************************/
#ifdef TEST
#include <stdio.h>

int main(int argc, char **argv)
{
   char     *filename = "test.pdb";
   WHOLEPDB *wpdb;
   FILE     *fp;
   PDB      *pdb1, *pdb2, *p;
   int      natoms;

   if((fp=fopen(filename, "r"))==NULL)
   {
      fprintf(stderr,"Can't open %s\n", filename);
      return(1);
   }
   
   if(((wpdb=blReadWholePDB(fp))==NULL)||(wpdb->pdb == NULL))
   {
      fprintf(stderr,"Can't read atoms from %s\n", filename);
      return(1);
   }
   fclose(fp);
   
   pdb1 = wpdb->pdb;
   if((pdb2 = blStripHPDBAsCopy(pdb1, &natoms))==NULL)
   {
      fprintf(stderr,"Stripping hydrogens failed\n");
      return(1);
   }
   
   if(!blCopyConects(pdb2, pdb1))
   {
      fprintf(stderr,"Failed to copy CONECT data\n");
      return(1);
   }

   printf(">>> Original CONECT data\n");
   wpdb->pdb = pdb1;
   blWriteWholePDBTrailer(stdout, wpdb, 1);

   /* Reset atom numbers in original list to prove we aren't picking up
      the old atoms
   */
   for(p=pdb1; p!=NULL; NEXT(p))
      p->atnum = 0;
   
   printf(">>> New CONECT data\n");
   wpdb->pdb = pdb2;
   blWriteWholePDBTrailer(stdout, wpdb, 1);
   
   return(0);
}
#endif


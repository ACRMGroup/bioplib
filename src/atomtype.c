/************************************************************************/
/**

   \file       atomtype.c
   
   \version    V2.0
   \date       22.07.15
   \brief      Set atom types in a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1999-2015
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
   This code sets atomtypes for PDB atoms. See pdb.h/ATOMTYPE_xxxx
   for the types

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0   23.03.99 Original  By: ACRM
-  V1.1   27.04.99 Fixes in SetPDBAtomTypesNSResidues() with O3* and P
-  V1.2   18.06.99 Fixes in SetPDBAtomTypesNSResidues() BOOL return and
                   force parameter
-  V1.3   21.06.99 Fixes in SetPDBAtomTypesNSResidues() Force now for
                   individual issues
-  V1.4   06.09.99 Fixes in SetPDBAtomTypesNSResidues() Further simple 
                   check for nucleotides.
-  V2.0   22.07.15 All force code etc removed - now returns a list
                   of warnings

*************************************************************************/
/* Includes
*/
#include "general.h"
#include "pdb.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 240

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static STRINGLIST *SetPDBAtomTypesNSResidues(PDB *pdb);
static void InitializePDBAtomTypes(PDB *pdb);
static void SetPDBAtomTypesModifiers(PDB *pdb);
static void SetPDBAtomTypesWaterAndNucleotides(PDB *pdb);
static void SetPDBAtomTypesMetals(PDB *pdb);

/************************************************************************/
/*>STRINGLIST *blSetPDBAtomTypes(PDB *pdb)
   ---------------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list
   \return                  STRINGLIST of any warning messages
                            NULL if all OK

-  21.07.15  Original   By: ACRM
*/
STRINGLIST *blSetPDBAtomTypes(PDB *pdb)
{
   InitializePDBAtomTypes(pdb);
   SetPDBAtomTypesMetals(pdb);
   SetPDBAtomTypesWaterAndNucleotides(pdb);
   SetPDBAtomTypesModifiers(pdb);
   return(SetPDBAtomTypesNSResidues(pdb));
}


/************************************************************************/
/*>static void SetPDBAtomTypesMetals(PDB *pdb)
   -------------------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list

   Identifies and sets metal atoms

-  21.07.15  Original   By: ACRM
*/
static void SetPDBAtomTypesMetals(PDB *pdb)
{
   PDB *p;
   /* Update HETATMs to metals and waters */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atomtype == ATOMTYPE_HETATM)
      {
         /* This is a list of non-metals in something like the order of
            likelihood of occurrence. We don't include noble gases since
            if these are found (unlikely!) they will be unbound and can
            be thought of as metals.
         */
         if(strcmp(p->element,"C")  &&
            strcmp(p->element,"N")  &&
            strcmp(p->element,"O")  &&
            strcmp(p->element,"H")  &&
            strcmp(p->element,"S")  &&
            strcmp(p->element,"P")  &&
            strcmp(p->element,"CL") &&
            strcmp(p->element,"BR") &&
            strcmp(p->element,"I")  &&
            strcmp(p->element,"F")  &&
            strcmp(p->element,"B")  &&
            strcmp(p->element,"SI") &&
            strcmp(p->element,"AS") &&
            strcmp(p->element,"SE") &&
            strcmp(p->element,"TE") &&
            strcmp(p->element,"AT"))
         {
            p->atomtype = ATOMTYPE_METAL;
         }
      }
   }
}


/************************************************************************/
/*>static void SetPDBAtomTypesWaterAndNucleotides(PDB *pdb)
   --------------------------------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list

   Sets atom types for waters and nucleotides - i.e. modifies the 
   ATOM and HETATM settings

-  21.07.15  Original   By: ACRM
*/
static void SetPDBAtomTypesWaterAndNucleotides(PDB *pdb)
{
   PDB *p;

   /* Update HETATMs to metals and waters */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(ISWATER(p))
      {
         p->atomtype = ATOMTYPE_WATER;
      }
      else if(p->atomtype == ATOMTYPE_ATOM)
      {
         if(!strncmp(p->resnam,"A  ",3) ||
            !strncmp(p->resnam,"C  ",3) ||
            !strncmp(p->resnam,"G  ",3) ||
            !strncmp(p->resnam,"I  ",3) ||
            !strncmp(p->resnam,"T  ",3) ||
            !strncmp(p->resnam,"Y  ",3) ||
            !strncmp(p->resnam,"U  ",3) ||
            !strncmp(p->resnam,"DA ",3) ||
            !strncmp(p->resnam,"DC ",3) ||
            !strncmp(p->resnam,"DT ",3) ||
            !strncmp(p->resnam,"DG ",3) ||
            !strncmp(p->resnam,"+A ",3) ||
            !strncmp(p->resnam,"+C ",3) ||
            !strncmp(p->resnam,"+G ",3) ||
            !strncmp(p->resnam,"+I ",3) ||
            !strncmp(p->resnam,"+T ",3) ||
            !strncmp(p->resnam,"+Y ",3) ||
            !strncmp(p->resnam,"+U ",3))
         {
            p->atomtype = ATOMTYPE_NUC;
         }
      }
   }
   
}


/************************************************************************/
/*>static void SetPDBAtomTypesModifiers(PDB *pdb)
   ----------------------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list

   Updates HETATMs to indicate if they are bound to ATOMs - i.e. they
   are residue modifiers

-  21.07.15  Original   By: ACRM
*/
static void SetPDBAtomTypesModifiers(PDB *pdb)
{
   PDB *p, *q; 
   int i;
   BOOL doneMod;

   /* Now look for connections between HETATMs and ATOMs */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atomtype == ATOMTYPE_HETATM)
      {
         for(i=0; i<p->nConect; i++)
         {
            q = p->conect[i];
            if(p->atomtype == ATOMTYPE_ATOM)
            {
               q->atomtype = ATOMTYPE_MODPROT;
            }
            else if(p->atomtype == ATOMTYPE_NUC)
            {
               q->atomtype = ATOMTYPE_MODNUC;
            }
         }
      }
   }
   
   /* Now update all the HETATMs that are connected to MODPROT or 
      MODNUC 
   */
   do
   {
      doneMod = FALSE;
      for(p=pdb; p!=NULL; NEXT(p))
      {
         for(i=0; i<p->nConect; i++)
         {
            q = p->conect[i];
            
            if(p->atomtype == ATOMTYPE_MODPROT &&
               q->atomtype == ATOMTYPE_HETATM)
            {
               q->atomtype = ATOMTYPE_MODPROT;
               doneMod = TRUE;
            }
            else if(p->atomtype == ATOMTYPE_MODNUC &&
                    q->atomtype == ATOMTYPE_HETATM)
            {
               q->atomtype = ATOMTYPE_MODNUC;
               doneMod = TRUE;
            }
            else if(p->atomtype == ATOMTYPE_BOUNDHET &&
                    q->atomtype == ATOMTYPE_HETATM)
            {
               q->atomtype = ATOMTYPE_BOUNDHET;
               doneMod = TRUE;
            }
         }
      }
   }  while(doneMod);
}


/************************************************************************/
/*>static STRINGLIST *SetPDBAtomTypesNSResidues(PDB *pdb)
   ------------------------------------------------------
*//**
   \param[in,out]     pdb    PDB linked list
   \return                   A STRINGLIST containing any warning 
                             messages

   Goes through BOUNDHET atoms and changes them to Non-standard residue
   atoms if they are linked via the backbone.

-  23.03.99 Original  By: ACRM
-  27.04.99 .O3* and .P.. were the wrong way round so NSNUC not being
            identified!
-  18.06.99 Added BOOL return and check on chain being correct (to catch
            1ubs) and force parameter
-  21.06.99 Force now works for individual issues rather than on/off
            for everything
-  06.09.99 Further simple check for nucleotides. If a non-standard
            residue (not recognised as a nucleotide) contains a 
            phosphorus and has a nucleotide on either side in the same
            chain, then set it to a nonstandard nucleotide. Fixes 1gsg
            (which is backbone only), 1ser (where the distance is too
            long to get flagged as boundhet).  
*/
static STRINGLIST *SetPDBAtomTypesNSResidues(PDB *pdb)
{
   STRINGLIST *warnings = NULL;
   PDB        *p, *q,
              *res1     = NULL,
              *res2     = NULL,
              *res3     = NULL,
              *res0     = NULL;
   int        nsVal     = 0;
   
   for(res1=pdb; res1!=NULL; res1=res2)
   {
      /* Replacement non-standard value. 0 indicates not found to be
         a non-standard residue
      */
      nsVal=0;
         
      if(res1!=NULL) res2 = blFindNextResidue(res1);
      if(res2!=NULL) res3 = blFindNextResidue(res2);

      /* If it's a bound het:
         If bound to following N or preceeding C, set type to NONSTDAA
         If bound to following P or preceeding O3*, set type to NONSTDNUC
      */
      if(res1->atomtype == ATOMTYPE_BOUNDHET)
      {
         for(p=res1; p!=res2; NEXT(p))
         {
            /* Search next residue                                      */
            for(q=res2; q!=res3; NEXT(q))
            {
               if(!strncmp(q->atnam,"N   ",4))
               {
                  if(blIsConected(p, q))
                  {
                     /* 18.06.99 Check they are in the same chain       */
                     if(!PDBCHAINMATCH(p, q))
                     {
                        char buffer[MAXBUFF];
                        sprintf(buffer,"Warning: Apparent \
non-standard amino acid has different chain\n\
          label from amino acid Nitrogen to which it is connected\n\
          Residue %s %s%d%s\n",
                                res1->resnam,
                                res1->chain,
                                res1->resnum,
                                res1->insert);
                        warnings = blStoreString(warnings, buffer);
                     }
                     
                     nsVal = ATOMTYPE_NONSTDAA;
                     p=NULL;
                     break;
                  }
               }
               if(!strncmp(q->atnam,"P   ",4))
               {
                  if(blIsConected(p, q))
                  {
                     /* 18.06.99 Check they are in the same chain       */
                     if(!PDBCHAINMATCH(p, q))
                     {
                        char buffer[MAXBUFF];
                        sprintf(buffer,"Warning: Apparent \
non-standard nucleotide has different chain\n\
          label from nucleotide phosphorus to which it is connected\n\
          Residue %s %s%d%s\n",
                                res1->resnam,
                                res1->chain,
                                res1->resnum,
                                res1->insert);
                        warnings = blStoreString(warnings, buffer);
                     }

                     nsVal = ATOMTYPE_NONSTDNUC;
                     p=NULL;
                     break;
                  }
               }
            }
            if(nsVal==0)
            {
               /* Search previous residue                               */
               for(q=res0; q!=NULL && q!=res1; NEXT(q))
               {
                  if(!strncmp(q->atnam,"C   ",4))
                  {
                     if(blIsConected(p, q))
                     {
                        /* 18.06.99 Check they are in the same chain    */
                        if(!PDBCHAINMATCH(p, q))
                        {
                           char buffer[MAXBUFF];
                           sprintf(buffer,"Warning: Apparent \
non-standard amino acid has different chain\n\
          label from amino acid Carbon to which it is connected\n\
          Residue %s %s%d%s\n",
                                   res1->resnam,
                                   res1->chain,
                                   res1->resnum,
                                   res1->insert);
                           warnings = blStoreString(warnings, buffer);
                        }

                        nsVal = ATOMTYPE_NONSTDAA;
                        p=NULL;
                        break;
                     }
                  }
                  if(!strncmp(q->atnam,"O3* ",4)) 
                  {
                     if(blIsConected(p, q))
                     {
                        /* 18.06.99 Check they are in the same chain    */
                        if(!PDBCHAINMATCH(p, q))
                        {
                           char buffer[MAXBUFF];
                           sprintf(buffer,"Warning: Apparent \
non-standard nucleotide has different chain\n\
          label from nucleotide O3* to which it is connected\n\
          Residue %s %s%d%s\n",
                                   res1->resnam,
                                   res1->chain,
                                   res1->resnum,
                                   res1->insert);
                           warnings = blStoreString(warnings, buffer);
                        }
                        nsVal = ATOMTYPE_NONSTDNUC;
                        p=NULL;
                        break;
                     }
                  }
               }
            }
            if(p==NULL)
               break;
         }
      }

      /* 06.09.99 Further simple check for nucleotides. If res1 is of
         type ATOM with res0 or res2 in the same chain and of type
         NUC, then we check that res1 contains a Phosphorus and, if
         so, we assume res1 is a nonstandard nucleotide. Fixes 1gsg
         (which is backbone only), 1ser (where the distance is too
         long to get flagged as boundhet).  
      */
      if(!nsVal &&                           /* Not done already       */
         (res1 != NULL) && (res1->atomtype == ATOMTYPE_ATOM) &&
         (((res0 != NULL) &&                 /* residue before         */
           ((res0->atomtype == ATOMTYPE_NUC) || 
            (res0->atomtype == ATOMTYPE_NONSTDNUC)) &&
           PDBCHAINMATCH(res0, res1)) ||
          ((res2 != NULL) &&                 /* residue after          */
           ((res2->atomtype == ATOMTYPE_NUC) || 
            (res2->atomtype == ATOMTYPE_NONSTDNUC)) &&
           PDBCHAINMATCH(res2, res1))))
      {
         for(p=res1; p!=res2; NEXT(p))
         {
            if(!strcmp(p->element,"P"))
            {
               nsVal = ATOMTYPE_NONSTDNUC;
               break;
            }
         }
      }

      /* If we found one of these links then modify all atoms in
         this residue
      */
      if(nsVal)
      {
         for(p=res1; p!=res2; NEXT(p))
         {
            p->atomtype = nsVal;
         }
      }

      res0=res1;
   }

   return(warnings);
}


/************************************************************************/
/*>static void InitializePDBAtomTypes(PDB *pdb)
   --------------------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list

   Initializes atom types based on ATOM or HETATM information

-  21.07.15  Original   By: ACRM
*/
static void InitializePDBAtomTypes(PDB *pdb)
{
   PDB *p;
   
   /* Initialize atom types based on record types   */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->record_type, "ATOM  ", 6))
      {
         p->atomtype = ATOMTYPE_ATOM;
      }
      else if(!strncmp(p->record_type, "HETATM", 6))
      {
         p->atomtype = ATOMTYPE_HETATM;
      }
      else
      {
         p->atomtype = ATOMTYPE_UNDEF;
      }
   }
}

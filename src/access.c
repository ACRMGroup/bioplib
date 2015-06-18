/************************************************************************/
/**

   \file       access.c
   
   \version    V1.2
   \date       17.06.15
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

   Calculation of solvent accessibility by the method of Lee and Richards.
   Based loosely on PMCL code by Peter McLaughlin

**************************************************************************

   Usage:
   ======
\code   
   RESRAD *blSetAtomRadii(PDB *pdb, FILE *fpRad)
\endcode
      Read the radius file and set the radii in the PDB linked list.
      Also returns the data from the radius file for use by
      blCalcResAccess()
      
\code
   BOOL blCalcAccess(PDB *pdb, int natoms, 
                     REAL integrationAccuracy, REAL probeRadius,
                     BOOL doAccessibility)
\endcode
      Does the accessibilty calculations. integrationAccuracy can be set
      to zero to use the default value

\code
   RESACCESS *blCalcResAccess(PDB *pdb, RESRAD *resrad)
\endcode
      Calculates residue accessibility and relative accessibility

**************************************************************************

   Revision History:
   =================
-  V1.0  21.04.99 Original   By: ACRM
-  V1.1  17.07.14 Extracted from XMAS code
-  V1.2  17.06.15 Added sidechain residues access

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Analyzing structures
   #FUNCTION blSetAtomRadii()
   Set atom radii from the radius file in the PDB linked list
   Returns the radius lookup information since it also contains the 
   standard accessibilities

   #FUNCTION  blCalcAccess()
   Allocates arrays and calls routines to populate them, do the access
   calculations and populate into the PDB linked list

   #FUNCTION  blCalcResAccess()
   Calculates and populates the residue totals and relative values
   using standards stored in resrad
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "macros.h"
#include "SysDefs.h"
#include "pdb.h"
#include "access.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF            160 /* I/O buffer                            */
#define MAX_INTERSECT    10000 /* Initial max no. of intersections of 
                                  neighbouring spheres - expands as
                                  required                              */
#define MAX_INTERSECT_EXP  100 /* Amount to expand intersect array by   */
#define MAX_ATOM_IN_CUBE   100 /* Initial max no. of atoms in a cube -
                                  expands as required                   */

#define FREE_ACCESS_STORAGE                                              \
do {                                                                     \
   if(cube)         free(cube);                                          \
   if(radii)        free(radii);                                         \
   if(radiiSquared) free(radiiSquared);                                  \
   if(neighbours)   free(neighbours);                                    \
   if(flag)         free(flag);                                          \
   if(arci)         free(arci);                                          \
   if(arcf)         free(arcf);                                          \
   if(deltaX)       free(deltaX);                                        \
   if(deltaY)       free(deltaY);                                        \
   if(dist)         free(dist);                                          \
   if(distSquared)  free(distSquared);                                   \
   if(atomTable)    free(atomTable);                                     \
   if(atomsInCube)  {                                                    \
      for(i=0;i<=maxAtomInCube;i++)                                      \
      {  if(atomsInCube[i]) free(atomsInCube[i]);                        \
         else        break;                                              \
      }                                                                  \
      free(atomsInCube);                                                 \
   }                                                                     \
}  while(0);


#define EXPAND_INTERSECT_ARRAYS                                          \
do {  int newMaxIntersect = (maxIntersect+MAX_INTERSECT_EXP),            \
          iexpand;                                                       \
   neighbours = (int  *)realloc(neighbours,                              \
                                (newMaxIntersect+1)*sizeof(int));        \
   flag = (int  *)realloc(flag, (newMaxIntersect+1)*sizeof(int));        \
   arci = (REAL *)realloc(arci, (newMaxIntersect+1)*sizeof(REAL));       \
   arcf = (REAL *)realloc(arcf, (newMaxIntersect+1)*sizeof(REAL));       \
   deltaX = (REAL *)realloc(deltaX, (newMaxIntersect+1)*sizeof(REAL));   \
   deltaY = (REAL *)realloc(deltaY, (newMaxIntersect+1)*sizeof(REAL));   \
   dist   = (REAL *)realloc(dist,   (newMaxIntersect+1)*sizeof(REAL));   \
   distSquared = (REAL *)realloc(distSquared,                            \
                                 (newMaxIntersect+1)*sizeof(REAL));      \
   if(neighbours==NULL || flag==NULL   || arci==NULL || arcf==NULL ||    \
      deltaX==NULL     || deltaY==NULL || dist==NULL ||                  \
      distSquared==NULL)                                                 \
   {  FREE_ACCESS_STORAGE;                                               \
      return(FALSE);                                                     \
   }                                                                     \
   for(iexpand=maxIntersect+1; iexpand<=newMaxIntersect; iexpand++)      \
   {  neighbours[iexpand]  = 0;                                          \
      flag[iexpand]        = 0;                                          \
      arci[iexpand]        = 0.0;                                        \
      arcf[iexpand]        = 0.0;                                        \
      deltaX[iexpand]      = 0.0;                                        \
      deltaY[iexpand]      = 0.0;                                        \
      dist[iexpand]        = 0.0;                                        \
      distSquared[iexpand] = 0.0;                                        \
   }                                                                     \
   maxIntersect = newMaxIntersect;                                       \
} while(0)
   
/************************************************************************/
/* Prototypes
*/
static void SortArcEndpoints(REAL *a, int n, int *flag);
static BOOL doCalcAccess(int numAtoms, REAL integrationAccuracy, 
                         REAL probeRadius, BOOL access,
                         REAL *AtomRadius,
                         REAL *x, REAL *y, REAL *z,
                         REAL *accessResults);
static void FillArrays(PDB *pdb, REAL *x, REAL *y, REAL *z, REAL *r);
static RESRAD *GetResidueRadii(RESRAD *resrad, char *resnam);
static RESRAD *ReadRadiusFile(FILE *fpRad);
static REAL GetStandardAccess(char *resnam, RESRAD *resrad);
static REAL DefaultRadius(char *element);
static void SetPDBAccess(PDB *pdb, REAL *accessArray);
static char *blGetElement(PDB *p);
static REAL GetStandardAccessSC(char *resnam, RESRAD *resrad);


/************************************************************************/
/*>static char *blGetElement(PDB *p)
   ---------------------------------
*//**
   \param[in]   *p   PDB structure pointer
   \return           Element assignment

   This will be replaced by Craig's routine in Bioplib
   Note: Atomic symbol is not stored in PDB data structure.
   Value set is based on columns 13-14 of pdb-formated text 
   file. 

-  17.07.14  Original   By: ACRM
 */
static char *blGetElement(PDB *p)
{
   char buffer[8];
   char *buffer_ptr;
   char *rawAtomName = p->atnam_raw;
   
   strncpy(buffer,rawAtomName,2);
   buffer[2] = '\0';
   
   KILLLEADSPACES(buffer_ptr,buffer);
      
   /* remove digits                                                     */
   if(strlen(buffer_ptr) == 2 && isdigit(buffer[1]))
   {
      buffer[1] = '\0';
   }
   else if(strlen(buffer_ptr) == 2 && !isalpha(buffer[0]))
   {
      buffer_ptr += 1;
   }
   
   /* fix hydrogens and carbons                                         */
   if(strlen(buffer_ptr) == 2 && rawAtomName[3] != ' ' &&
      (rawAtomName[0] == 'H' || rawAtomName[0] == 'C' ||
       rawAtomName[0] == 'N' || rawAtomName[0] == 'O' ||
       rawAtomName[0] == 'P'))
   {
      if(!isalpha(rawAtomName[2]) || !isalpha(rawAtomName[3]))
      {
         buffer[1] = '\0';
      }
   }
   
   return(buffer_ptr);
}

/************************************************************************/
/*>RESRAD *blSetAtomRadii(PDB *pdb, FILE *fpRad)
   ---------------------------------------------
*//**
   \param[in,out] *pdb       PDB linked list
   \param[in]     *fpRad     Radius file pointer
   \return                   Linked list of radius information

   Set atom radii from the radius file in the PDB linked list
   Returns the radius lookup information since it also contains the 
   standard accessibilities

-  22.04.99 Original   By: ACRM
-  16.06.99 Initialise radii to NULL
-  22.06.99 Changed to call DefaultRadius() with element type rather
            than atom name
-  16.07.14 Rewritten to work outside XMAS format and now takes the
            file pointer to the radii file rather than the filename
*/
RESRAD *blSetAtomRadii(PDB *pdb, FILE *fpRad)
{
   RESRAD *resrad,
          *radii = NULL;
   int    i;
   char   currentResnam[8];
   PDB    *p;


   resrad=ReadRadiusFile(fpRad);

   strcpy(currentResnam,"    ");
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If residue name has changed, find pointer to this residue in
         the radius list
      */
      if(strncmp(p->resnam, currentResnam, 3))
      {
         strcpy(currentResnam, p->resnam);
         radii = GetResidueRadii(resrad,p->resnam);
      }

      if(radii == NULL)
      {
         /* Residue not found, use default for element                  */
         p->radius = DefaultRadius(blGetElement(p));
      }
      else
      {
         for(i=0; i<radii->natoms; i++)
         {
            if(!strncmp(p->atnam_raw, radii->atnam[i], 4))
            {
               p->radius = radii->radius[i];
               break;
            }
         }
         
         if(i==radii->natoms)
         {
            /* Didn't find this atom in the residue - use default atom 
               radius
            */
            p->radius = DefaultRadius(blGetElement(p));
         }
      }
   }

   return(resrad);
}


/************************************************************************/
/*>BOOL blCalcAccess(PDB *pdb, int natoms, 
                     REAL integrationAccuracy, 
                     REAL probeRadius,
                     BOOL doAccessibility)
   --------------------------------------------------------------
*//**
   \param[in,out]    *pdb                  PDB linked list
   \param[in]        natoms                Number of atoms
   \param[in]        integrationAccuracy   Integration accuracy
   \param[in]        probeRadius           Probe radius
   \param[in]        doAccessibility       Accessibility or contact area
   \return                                 Success

   Allocates arrays and calls routines to populate them, do the access
   calculations and populate into the PDB linked list

-  22.04.99 Original   By: ACRM
*/
BOOL blCalcAccess(PDB *pdb, int natoms, 
                  REAL integrationAccuracy, REAL probeRadius,
                  BOOL doAccessibility)
{
   REAL *x = NULL, 
        *y = NULL, 
        *z = NULL, 
        *radii = NULL, 
        *accessArray = NULL;
   BOOL retval = FALSE;

   if(integrationAccuracy < VERY_SMALL)
      integrationAccuracy = ACCESS_DEF_INTACC;
   
   /* Allocate arrays                                                   */
   if((x=(REAL *)malloc(natoms * sizeof(REAL)))!=NULL)
   {
      if((y=(REAL *)malloc(natoms * sizeof(REAL)))!=NULL)
      {
         if((z=(REAL *)malloc(natoms * sizeof(REAL)))!=NULL)
         {
            if((radii=(REAL *)malloc(natoms * sizeof(REAL)))!=NULL)
            {
               if((accessArray=(REAL *)malloc(natoms * sizeof(REAL)))
                  !=NULL)
               {
                  retval = TRUE;
                  
                  /* Populate arrays from PDB structure, do the 
                     accessibility run and put the results back into the 
                     PDB structure
                  */
                  FillArrays(pdb, x, y, z, radii);
                  doCalcAccess(natoms, integrationAccuracy, probeRadius, 
                               doAccessibility,
                               radii, x, y, z,
                               accessArray);
                  SetPDBAccess(pdb, accessArray);
               }
            }
         }
      }
   }
   
   /* Free the allocated memory                                         */
   if(x!=NULL)           free(x);
   if(y!=NULL)           free(y);
   if(z!=NULL)           free(z);
   if(radii!=NULL)       free(radii);
   if(accessArray!=NULL) free(accessArray);
   
   return(retval);
}



/************************************************************************/
/*>static void FillArrays(PDB *pdb, REAL *x, REAL *y, REAL *z, REAL *r)
   --------------------------------------------------------------------
*//**
   \param[in]   *pdb      PDB linked list
   \param[out]  *x        X coordinate array
   \param[out]  *y        Y coordinate array
   \param[out]  *z        Z coordinate array
   \param[out]  *r        radius array

   Fills the coordinate and atom arrays from the PDB linked list

-  22.04.99 Original   By: ACRM
*/
static void FillArrays(PDB *pdb, REAL *x, REAL *y, REAL *z, REAL *r)
{
   PDB *p;
   int  count = 0;
   

   for(p=pdb; p!=NULL; NEXT(p))
   {
      x[count] = p->x;
      y[count] = p->y;
      z[count] = p->z;
      r[count] = p->radius;
      
      count++;
   }
}


/************************************************************************/
/*>static RESRAD *GetResidueRadii(RESRAD *resrad, char *resnam)
   --------------------------------------------------------------
*//**
   \param[in]   *resrad    Linked list of atom radii information
   \param[in]   *resnam    Residue name we are looking for
   \return                 Pointer to information on this residue
                           (or NULL)

   Gets the radius information for the specified residue.

-  22.04.99 Original   By: ACRM
*/
static RESRAD *GetResidueRadii(RESRAD *resrad, char *resnam)
{
   RESRAD *r;

   /* Search through the residue types to find this residue             */
   for(r=resrad; r!=NULL; NEXT(r))
   {
      if(!strncmp(r->resnam, resnam, 3))
      {
         return(r);
         break;
      }
   }

   return(NULL);
}

/************************************************************************/
/*>static RESRAD *ReadRadiusFile(FILE *fpRad)
   -----------------------------------------
*//**
   \param[in]   *fpRad    Atom radius file pointer
   \return                Linked list of residue and atom radius
                          information

   Reads the specified file into the linked list of atom radii for each
   residue type. The file also contains standard residue accessibility.

-  22.04.99 Original   By: ACRM
-  16.06.99 Initialise atomIndex to 0 and r to NULL
-  16.07.14 Changed to passing in file pointer
-  17.06.15 Reads stdAccessSc
*/
static RESRAD *ReadRadiusFile(FILE *fpRad)
{
   char   buffer[MAXBUFF],
          *chp,
          junk[8];
   RESRAD *resrad = NULL, 
          *r = NULL;
   int    atomCount = 0,
          atomIndex = 0;

   while(fgets(buffer,MAXBUFF,fpRad))
   {
      if((chp=strchr(buffer,'#'))!=NULL)  /* Strip comments             */
         *chp = '\0';

      KILLLEADSPACES(chp, buffer);
      if(!strlen(chp))                    /* Skip blank lines           */
         continue;

      /* Beginning of a new residue                                     */
      if(!atomCount)
      {
         if(resrad==NULL)
         {
            INIT(resrad, RESRAD);
            r = resrad;
         }
         else
         {
            ALLOCNEXT(r, RESRAD);
         }
         if(r==NULL)
         {
            FREELIST(resrad, RESRAD);
            return(NULL);
         }
         
         sscanf(buffer,"%s %d %lf %lf", 
                r->resnam, &(r->natoms), 
                &(r->stdAccess), &(r->stdAccessSC));
         atomCount = r->natoms;
         atomIndex = 0;
      }
      else
      {
         strncpy(r->atnam[atomIndex],buffer,5);
         r->atnam[atomIndex][5] = '\0';
#ifdef NO_LEAD_SPACE
         if(r->atnam[atomIndex][0] == ' ')
         {
            for(i=0; i<4; i++)
               r->atnam[atomIndex][i] = r->atnam[atomIndex][i+1];
         }
#endif
         /* Replace dots with spaces                                    */
         DEDOTIFY(r->atnam[atomIndex]);
         
         sscanf(buffer,"%s %lf",junk,&(r->radius[atomIndex]));
         atomIndex++;
         atomCount--;
      }
   }
   
   return(resrad);
}


/************************************************************************/
/*>static REAL GetStandardAccess(char *resnam, RESRAD *resrad)
   ------------------------------------------------------------
*//**
   \param[in]   *resnam  Residue name
   \param[in]   *resrad  Residue/atom radii and standard 
                         accessibilities
   \return               Standard accessibility for this residue

   Gets the standard accessibility for the specified residue type

-  22.04.99 Original   By: ACRM
*/
static REAL GetStandardAccess(char *resnam, RESRAD *resrad)
{
   RESRAD *r;
   
   /* Search through the residue types to find this residue             */
   for(r=resrad; r!=NULL; NEXT(r))
   {
      if(!strncmp(r->resnam, resnam, 3))
      {
         return(r->stdAccess);
      }
   }

   return((REAL)0.0);
}

/************************************************************************/
/*>static REAL GetStandardAccessSC(char *resnam, RESRAD *resrad)
   --------------------------------------------------------------
*//**
   \param[in]   *resnam  Residue name
   \param[in]   *resrad  Residue/atom radii and standard 
                         accessibilities
   \return               Standard s/c accessibility for this residue

   Gets the standard sidechain accessibility for the specified residue
   type

-  17.06.15 Original   By: ACRM
*/
static REAL GetStandardAccessSC(char *resnam, RESRAD *resrad)
{
   RESRAD *r;
   
   /* Search through the residue types to find this residue             */
   for(r=resrad; r!=NULL; NEXT(r))
   {
      if(!strncmp(r->resnam, resnam, 3))
      {
         return(r->stdAccessSC);
      }
   }

   return((REAL)0.0);
}

/************************************************************************/
/*>static REAL DefaultRadius(char *element)
   ----------------------------------------
*//**
   \param[in]   *element   Element
   \return                 Default radius

   Returns the default radius for a specified atom type

-  22.04.99 Original   By: ACRM
-  22.06.99 Modified to work with element rather than atom name
*/
static REAL DefaultRadius(char *element)
{
   if(!strcmp(element,"C"))
      return((REAL)1.80);
   else if(!strcmp(element,"N"))
      return((REAL)1.60);
   else if(!strcmp(element,"S"))
      return((REAL)1.85);
   else if(!strcmp(element,"O"))
      return((REAL)1.40);
   else if(!strcmp(element,"P"))
      return((REAL)1.90);
   else if(!strcmp(element,"CA"))
      return((REAL)2.07);
   else if(!strcmp(element,"FE"))
      return((REAL)1.47);
   else if(!strcmp(element,"CU"))
      return((REAL)1.78);
   else if(!strcmp(element,"ZN"))
      return((REAL)1.39);
   else if(!strcmp(element,"MG"))
      return((REAL)1.73);
   else 
      return((REAL)1.80);
}


/************************************************************************/
/*>static void SetPDBAccess(PDB *pdb, REAL *accessArray)
   -----------------------------------------------------
*//**
   \param       *pdb           PDB linked list
   \param[in]   *accessArray   Accessibility information

   Puts accessibility information from the array into the PDB linked list.

-  16.07.14 Original   By: ACRM
*/
static void SetPDBAccess(PDB *pdb, REAL *accessArray)
{
   PDB *p;
   int i = 0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->access = accessArray[i++];
   }
}




/************************************************************************/
/*>RESACCESS *blCalcResAccess(PDB *pdb, RESRAD *resrad)
   ----------------------------------------------
*//**
   \param[in,out] *pdb      PDB linked list
   \param[in]     *resrad   Linked list of atom radius information
   \return                  Linked list of residue accessibilities

   Calculates and populates the residue totals and relative values
   using standards stored in resrad

-  22.04.99 Original   By: ACRM
-  04.08.99 Changed check on return from GetStandardAccess() to check
            for < VERY_SMALL rather than ==0.0
            Set relative access to -1 if the standard accessibility is
            unknown rather than to 0.0
-  17.06.15 Added calculation of sidechain accessibility
*/
RESACCESS *blCalcResAccess(PDB *pdb, RESRAD *resrad)
{
   PDB *start, *stop, *p;
   REAL resAccess,
        relAccess,
        scAccess,
        scRelAccess,
        stdAccess,
        stdAccessSC;

   RESACCESS *residues = NULL, 
             *r = NULL;

   for(start=pdb; start!=NULL; start=stop)
   {
      stop = blFindNextResidue(start);

      /* Add up accessibility for this residue                          */
      resAccess = (REAL)0.0;
      scAccess  = (REAL)0.0;
      for(p=start; p!=stop; NEXT(p))
      {
         resAccess += p->access;
         if(strncmp(p->atnam, "N   ", 4) &&
            strncmp(p->atnam, "CA  ", 4) &&
            strncmp(p->atnam, "C   ", 4) &&
            strncmp(p->atnam, "O   ", 4) &&
            strncmp(p->atnam, "OXT ", 4))
         {
            scAccess += p->access;
         }
      }

      /* Get the standard accessibility for this amino acid and calculate
         relative accessibility
      */
      stdAccess   = GetStandardAccess(start->resnam, resrad);
      stdAccessSC = GetStandardAccessSC(start->resnam, resrad);

      if(stdAccess<VERY_SMALL)
         relAccess   = -1.0;
      else
         relAccess   = 100.0 * resAccess / stdAccess;

      if(stdAccessSC < VERY_SMALL)
         scRelAccess = -1.0;
      else
         scRelAccess = 100.0 * scAccess  / stdAccessSC;
      
      /* Create space to store the values                               */
      if(residues == NULL)
      {
         INIT(residues, RESACCESS);
         r = residues;
      }
      else
      {
         ALLOCNEXT(r, RESACCESS);
      }
      if(r==NULL)
      {
         FREELIST(residues, RESACCESS);
         return(NULL);
      }

      /* and store them                                                 */
      strcpy(r->resnam, start->resnam);
      strcpy(r->insert, start->insert);
      strcpy(r->chain,  start->chain);
      r->resnum      = start->resnum;
      r->resAccess   = resAccess;
      r->relAccess   = relAccess;
      r->scAccess    = scAccess;
      r->scRelAccess = scRelAccess;
   }

   return(residues);
}


/************************************************************************/
/*>static void SortArcEndpoints(REAL *endpoints, int nEndpoints, 
                                int *flag)
   -------------------------------------------------------------
*//**
   \param[in,out]    *endpoints     Arc endpoints for sorting
   \param[in]        nEndpoints     Number of endpoints
   \param[out]       *flag          Flagged endpoints

   Sorts the arc endpoints

-  21.04.99 Original   By: ACRM
*/
static void SortArcEndpoints(REAL *endpoints, int nEndpoints, int *flag)
{
   int   tmpFlagVal, midpoint,
         il[16+1],
         iu[16+1],
         i, j, k, l, m;
   REAL  tmpEndpointVal, tmpEndpointVal2;
   
   for(i=1; i<=nEndpoints; i++)
      flag[i]=i;

   m=1;
   i=1;
   j=nEndpoints;
   
jump5:
   if(i < j)
   {

jump10:
      k=i;
      midpoint=(j+i)/2;
      tmpEndpointVal=endpoints[midpoint];

      if(endpoints[i] > tmpEndpointVal)
      {
         endpoints[midpoint]= endpoints[i];
         endpoints[i]=tmpEndpointVal;
         tmpEndpointVal=endpoints[midpoint];
         tmpFlagVal=flag[midpoint];
         flag[midpoint]=flag[i];
         flag[i]=tmpFlagVal;
      }

      l=j;

      if(endpoints[j] < tmpEndpointVal)
      {
         endpoints[midpoint]=endpoints[j];
         endpoints[j]=tmpEndpointVal;
         tmpEndpointVal=endpoints[midpoint];
         tmpFlagVal=flag[midpoint];
         flag[midpoint]=flag[j];
         flag[j]=tmpFlagVal;
         
         if(endpoints[i] > tmpEndpointVal)
         {
            endpoints[midpoint]=endpoints[i];
            endpoints[i]=tmpEndpointVal;
            tmpEndpointVal=endpoints[midpoint];
            tmpFlagVal=flag[midpoint];
            flag[midpoint]=flag[i];
            flag[i]=tmpFlagVal;
            
            goto exit40;
            
jump30:
            endpoints[l]=endpoints[k];
            endpoints[k]=tmpEndpointVal2;
            tmpFlagVal=flag[l];
            flag[l]=flag[k];
            flag[k]=tmpFlagVal;
exit40:
            ;
         }
      }
      
      do
      {
         l--;
      }  while((l>0) && (endpoints[l] > tmpEndpointVal));
      tmpEndpointVal2=endpoints[l];
      
      do
      {
         k++;
      }  while((k<=nEndpoints) && (endpoints[k] < tmpEndpointVal));

      if(k<=l) goto jump30;

      if((l-i) > (j-k))
      {
         il[m]=i;
         iu[m]=l;
         i=k;
         m++;
         
         goto exit80;
      }

      il[m]=k;
      iu[m]=j;
      j=l;
      m++;
      
      goto exit80;
   }

   do
   {
      m--;
      if(m==0) return;
      i=il[m];
      j=iu[m];

exit80:
      if(j-i>=1) goto jump10;
      if(i==1) goto jump5;

      i--;

jump90:
      i++;
   }  while(i==j);
   
   tmpEndpointVal=endpoints[i+1];
   if(endpoints[i]<=tmpEndpointVal) goto jump90;
   
   tmpFlagVal=flag[i+1];
   k=i;
   
   do
   {
      endpoints[k+1]=endpoints[k];
      flag[k+1]=flag[k];
      k--;
   }  while(tmpEndpointVal < endpoints[k]);
 
   endpoints[k+1]=tmpEndpointVal;
   flag[k+1]=tmpFlagVal;
   
   goto jump90;
}
      

/************************************************************************/
/*>static BOOL doCalcAccess(int numAtoms, REAL integrationAccuracy,
                            REAL probeRadius, 
                            BOOL access, REAL *atomRadii,
                            REAL *x, REAL *y, REAL *z,
                            REAL *accessResults)
   --------------------------------------------------------------------
*//**
   \param[in]   numAtoms             Number of atoms
   \param[in]   integrationAccuracy  Integration accuracy
   \param[in]   probeRadius          Radius of probe atom
   \param[in]   access               Do solvent accessibility
                                     rather than contact surface
   \param[in]   *atomRadii           Array of atom radii
   \param[in]   *x                   Array of x coordinates
   \param[in]   *y                   Array of y coordinates
   \param[in]   *z                   Array of z coordinates
   \param[out]  *accessResults       Array of accessibility results
   \return                           Success?

   Does the real work of calculating accessibility

-  21.04.99 Original   By: ACRM
-  08.06.99 Fixed allocation of second dimension of atomsInCube[][] 
            to njidim rather than numAtoms
*/
static BOOL doCalcAccess(int numAtoms, REAL integrationAccuracy,
                         REAL probeRadius, 
                         BOOL access, REAL *atomRadii,
                         REAL *x, REAL *y, REAL *z,
                         REAL *accessResults)
{
   int   *cube   = NULL,
         *neighbours   = NULL,
         *flag   = NULL,
         *atomTable   = NULL,
         **atomsInCube  = NULL;
   REAL  *radii  = NULL, *radiiSquared=NULL,
         *arci   = NULL, *arcf=NULL,
         *deltaX = NULL, *deltaY=NULL,
         *dist   = NULL, *distSquared=NULL;
         
   int   i, j, k, l, m, n,
         jj, kk,
         cubeAtom, io, keyAtom,
         cubeIndex,
         nzp, karc, idim, jidim, kjidim,
         mkji, nm,
         maxIntersect  = MAX_INTERSECT,
         maxAtomInCube = MAX_ATOM_IN_CUBE;
   REAL  xmin  =  999999.0,    
         ymin  =  999999.0,
         zmin  =  999999.0,
         xmax  = -999999.0,
         ymax  = -999999.0,
         zmax  = -999999.0,
         pi    = acos(-1.0),
         twoPi = 2.0*acos(-1.0),
         maxRadius, totalArea, tmpArea,
         intersect,
         xr, yr, zr, 
         radius, radiusX2, radiusSquared, 
         zres, zgrid,
         t, tf, ti, tt, 
         partialArea,
         rsec2r, rsecr, 
         rsec2n, rsecn,
         alpha, beta,
         arcsum;
   BOOL  SkipAccess = FALSE;

#ifdef DEBUG
   int   maxAtomsSeenInCube = 0,
         maxIntersectsSeen  = 0;
#endif

   /* Reset arrays to count from 1 instead of 0                         */
   atomRadii--;
   x--; y--; z--;
   accessResults--;

   /* Allocate memory for arrays based on number of atoms               */
   cube         = (int  *)malloc((numAtoms+1)*sizeof(int));
   radii        = (REAL *)malloc((numAtoms+1)*sizeof(REAL));
   radiiSquared = (REAL *)malloc((numAtoms+1)*sizeof(REAL));
   
   /* Allocate memory for arrays based on number of intersects          */
   neighbours        = (int  *)malloc((MAX_INTERSECT+1)*sizeof(int));
   flag        = (int  *)malloc((MAX_INTERSECT+1)*sizeof(int));
   arci        = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));
   arcf        = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));
   deltaX      = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));
   deltaY      = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));
   dist        = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));
   distSquared = (REAL *)malloc((MAX_INTERSECT+1)*sizeof(REAL));

   /* Check allocations                                                 */
   if(cube         == NULL ||
      radii        == NULL ||
      radiiSquared == NULL ||
      neighbours   == NULL ||
      flag         == NULL ||
      arci         == NULL ||
      arcf         == NULL ||
      deltaX       == NULL ||
      deltaY       == NULL ||
      dist         == NULL ||
      distSquared  == NULL)
   {
      FREE_ACCESS_STORAGE;
      return(FALSE);
   }

   /* Find the limits of the surrounding box                            */
   maxRadius = 0.0;
   for(i=1; i<=numAtoms; i++)
   {
      radii[i]         = atomRadii[i] + probeRadius;
      radiiSquared[i]  = radii[i] * radii[i];
      accessResults[i] = 0.0;

      if(radii[i] > maxRadius) maxRadius=radii[i];
      if(x[i] < xmin)          xmin=x[i];
      if(y[i] < ymin)          ymin=y[i];
      if(z[i] < zmin)          zmin=z[i];
      if(x[i] > xmax)          xmax=x[i];
      if(y[i] > ymax)          ymax=y[i];
      if(z[i] > zmax)          zmax=z[i];
   }
   maxRadius *= 2.0;

   /* Set up cubes containing the atoms. The dimension of a cube edge is 
      equal to the radius of the largest atom sphere.
   */
   idim = (xmax-xmin)/maxRadius + 1.0;
   if(idim < 3) idim = 3;

   jidim = (ymax-ymin)/maxRadius + 1.0;
   if(jidim < 3) jidim = 3;
   jidim *= idim;

   kjidim = (zmax-zmin)/maxRadius + 1.0;
   if(kjidim < 3) kjidim = 3;
   kjidim *= jidim;

#ifdef DEBUG
   fprintf(stderr,"Number of cubes: %d\n", kjidim);
#endif

   /* Allocate memory for atomsInCube 2D array
      08.06.99 Corrected to inner dimension being kjidim rather than
               numAtom
   */
   if((atomsInCube = (int **)malloc((MAX_ATOM_IN_CUBE+1)*sizeof(int *)))
      ==NULL)
   {
      FREE_ACCESS_STORAGE;
      return(FALSE);
   }

   for(i=0; i<=MAX_ATOM_IN_CUBE; i++)
   {
      if((atomsInCube[i] = (int *)malloc((kjidim+1)*sizeof(int)))==NULL)
      {
         FREE_ACCESS_STORAGE;
         return(FALSE);
      }
   }

   /* Prepare the cubes
      -----------------
      Allocate memory for the cubes. Each cube may contain upto 
      MAX_ATOM_IN_CUBE atoms. We count through the cubes with cubeIndex 
      and store the atom indices in the atomTable[] array.
   */
   if((atomTable = (int *)malloc((kjidim+1)*sizeof(int)))==NULL)
   {
      FREE_ACCESS_STORAGE;
      return(FALSE);
   }
   for(l=1; l<=kjidim; l++)
      atomTable[l]=0;

   for(l=1; l<=numAtoms; l++)
   {
      i = (x[l]-xmin)/maxRadius + 1.0;
      j = (y[l]-ymin)/maxRadius;
      k = (z[l]-zmin)/maxRadius;

      cubeIndex = k*jidim + j*idim + i;
      n         = atomTable[cubeIndex] + 1;

#ifdef DEBUG
      if(n > maxAtomsSeenInCube)
         maxAtomsSeenInCube = n;
#endif

      /* If we have too many atoms in the cube, expand the atomsInCube 
         array   
      */
      if(n > maxAtomInCube)
      {
         int newMaxAtomInCube = maxAtomInCube + MAX_ATOM_IN_CUBE,
             iexpand;
         
         if((atomsInCube = (int **)realloc(atomsInCube, 
                                    (newMaxAtomInCube+1)*sizeof(int *)))
            ==NULL)
         {
            FREE_ACCESS_STORAGE;
            return(FALSE);
         }

         for(iexpand=maxAtomInCube+1; 
             iexpand<=newMaxAtomInCube; 
             iexpand++)
         {
            /* 08.06.99 Corrected numAtoms to kjidim                    */
            if((atomsInCube[iexpand] = 
                (int *)malloc((kjidim+1) * sizeof(int)))
               ==NULL)
            {
               FREE_ACCESS_STORAGE;
               return(FALSE);
            }
         }
         
         maxAtomInCube = newMaxAtomInCube;
      }

      atomTable[cubeIndex]      = n;
      atomsInCube[n][cubeIndex] = l;
      cube[l]                   = cubeIndex;
   }

#ifdef DEBUG
   fprintf(stderr,"Max number of atoms in a cube: %d\n", 
           maxAtomsSeenInCube);
#endif


   /* Perform the actual accessibility calculations 
      ---------------------------------------------
      We cycle through each atom in turn
   */
   for(keyAtom=1; keyAtom<=numAtoms; keyAtom++)
   {
      cubeIndex     = cube[keyAtom];
      io            = 0;
      totalArea     = 0.0;
      xr            = x[keyAtom];
      yr            = y[keyAtom];
      zr            = z[keyAtom];
      radius        = radii[keyAtom];
      radiusX2      = radius*2.0;
      radiusSquared = radiiSquared[keyAtom];
      
      /* Find the 'mkji' cubes neighboring the cubeIndex cube           */
      for(kk=1; kk<=3; kk++)
      {
         k=kk-2;
         
         for(jj=1; jj<=3; jj++)
         {
            j=jj-2;
            
            for(i=1; i<=3; i++)
            {
               mkji = cubeIndex + k*jidim + j*idim + i - 2;
               if(mkji >= 1)
               {
                  if(mkji > kjidim)
                  {
                     jj=5;    /* Force exit from for(jj) loop           */
                     kk=5;    /* Force exit from for(kk) loop           */
                     break;   /* out of for(i) loop                     */
                  }
                  
                  nm = atomTable[mkji];
                  
                  if(nm >= 1)
                  {
                     /* Create neighbours[] which is a list of the atoms 
                        which are neighbours of atom keyAtom
                     */
                     for(m=1; m<=nm; m++)
                     {
                        cubeAtom = atomsInCube[m][mkji];
                        if(cubeAtom != keyAtom)
                        {
                           io++;
#ifdef DEBUG
                           if(io > maxIntersectsSeen)
                              maxIntersectsSeen = io;
#endif
                           if(io > maxIntersect)
                           {
                              EXPAND_INTERSECT_ARRAYS;
                           }
                           
                           deltaX[io] = xr - x[cubeAtom];
                           deltaY[io] = yr - y[cubeAtom];
                           
                           distSquared[io] = deltaX[io]*deltaX[io] +
                                             deltaY[io]*deltaY[io];
                           dist[io]        = sqrt(distSquared[io]);
                           neighbours[io]        = cubeAtom;
                        }
                     }
                  }
               }
            }
         }
      }
      
      if(io == 0)
      {
         totalArea = twoPi * radiusX2;
      }
      else
      {
         /* Calculate the z resolution                                  */
         nzp   = 1.0/integrationAccuracy + 0.5;
         zres  = radiusX2 / nzp;
         zgrid = z[keyAtom] - radius - zres/2.0;
         
         /* Take a section of atom spheres which is perpendicular to the 
            z axis 
         */
         for(i=1; i<=nzp; i++)
         {
            zgrid += zres;
            
            /* Calculate the radius of the circle of intersection of the 
               keyAtom sphere on the current z-plane 
            */
            rsec2r = radiusSquared - (zgrid-zr)*(zgrid-zr);
            rsecr  = sqrt(rsec2r);
            
            for(k=1; k<=maxIntersect; k++)
               arci[k] = 0.0;
            
            karc=0;
            
            for(j=1; j<=io; j++)
            {
               cubeAtom = neighbours[j];
               
               /* Find the radius of the circle locus                   */
               rsec2n = radiiSquared[cubeAtom] - 
                  (zgrid-z[cubeAtom])*(zgrid-z[cubeAtom]);
               
               if(rsec2n > 0.0)
               {
                  rsecn = sqrt(rsec2n);
                  
                  /* Find the intersections of the n circles with the 
                     keyAtom circles in this section 
                  */
                  if(dist[j] < (rsecr+rsecn))
                  {
                     /* Test whether the the circles intersect, or
                        whether one circle is completely inside the
                        other in which case we don't calculate the
                        accessibility!  
                     */
                     intersect = rsecr - rsecn;
                     
                     SkipAccess = FALSE;
                     if(dist[j] <= fabs(intersect)) 
                     {
                        if(intersect <= 0.0)
                        {
                           SkipAccess = TRUE;
                           break;               /* out of the j loop    */
                        }
                        continue;               /* the j loop           */
                     }
                     
                     /* Expand the intersect arrays if we have too many */
                     if(++karc >= maxIntersect)
                     {
                        EXPAND_INTERSECT_ARRAYS;
                     }
#ifdef DEBUG
                     if(karc > maxIntersectsSeen)
                        maxIntersectsSeen = karc;
#endif
                     
                     /* If the circles do intersect, then we find the 
                        points of intersection.
                        
                        The initial and final arc endpoints are found
                        for the keyAtom circle intersected by a 
                        neighboring circle contained in the same plane. 
                        The initial endpoint of the enclosed arc is stored
                        in arci, and the final arc in arcf. This uses
                        the cosine law.
                        
                        Calculate alpha which is the angle between a
                        line containing a point of intersection, the
                        reference circle center and the line
                        containing both circle centers.  
                     */
                     alpha = acos((distSquared[j] + rsec2r - rsec2n) /
                                  (2.0 * dist[j] * rsecr));
                     
                     /* Calculate beta which is the angle between the
                        line containing both circle centers and the
                        x-axis 
                     */
                     beta = atan2(deltaY[j], deltaX[j]) + pi;
                     
                     ti = beta - alpha;
                     tf = beta + alpha;
                     if(ti < 0.0)
                        ti += twoPi;
                     
                     if(tf > twoPi)
                        tf -= twoPi;
                     
                     arci[karc] = ti;
                     
                     /* If the arc crosses zero, then it is broken
                        into two segments. The first ends at twoPi and
                        the second begins at zero 
                     */
                     if(tf < ti) 
                     {
                        arcf[karc] = twoPi;
                        karc++;
                     }
                     
                     arcf[karc] = tf;
                  }
               }
            }
            
            if(SkipAccess)
            {
               /* Accessibility skipped because circle was within another
                  one 
               */
               SkipAccess = FALSE;
            }
            else
            {
               /* Find the accessible contact surface area for the
                  sphere keyAtom on this section 
               */
               if(karc == 0)
               {
                  arcsum = twoPi;
               }
               else
               {
                  /* Sort the arc endpoints on the value of the
                     initial arc endpoint 
                  */
                  SortArcEndpoints(arci, karc, flag);
                  
                  /* Calculate the length of the accessible arc         */
                  arcsum = arci[1];
                  t      = arcf[flag[1]];
                  
                  if(karc != 1)
                  {
                     for(k=2; k<=karc; k++)
                     {
                        if(t < arci[k])
                        {
                           arcsum += (arci[k]-t);
                        }
                        
                        tt = arcf[flag[k]];
                        if(tt > t)
                        {
                           t = tt;
                        }
                     }
                  }
                  
                  arcsum += (twoPi-t);
               }
               
               /* Calculate the partial accessible area for this atom
                  on this section. The area/radius is equal to the
                  accessible arc length x the section thickness.  
               */
               partialArea = arcsum * zres;
               
               /* ...and add this to the total area for this atom       */
               totalArea += partialArea;
            }
         }
      }
      
      /* Scale the area to Van der Waals shell                          */
      tmpArea = totalArea * (radius-probeRadius) * (radius-probeRadius) /
         radius;
      
      /* Convert from the contact area to the accessible surface area
         if required 
      */
      if(access)
      {
         tmpArea *= (radius*radius) / 
            ((radius-probeRadius) * (radius-probeRadius));
      }
      
      accessResults[keyAtom] = tmpArea;
   }

#ifdef DEBUG
   fprintf(stderr,"Maximum intersects: %d\n",maxIntersectsSeen);
#endif

   FREE_ACCESS_STORAGE;
   
   return(TRUE);
}


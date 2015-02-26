/************************************************************************/
/**

   \file       RdSSPDB.c
   
   \version    V1.2
   \date       26.02.15
   \brief      Read disulphide information from header records of a PDB
               file
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1993-2015
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

   Read the header records of a PDB file to find disulphide information.
   Builds a linked list of type DISULPHIDE.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   14.10.93 Original   By: ACRM
   V1.1   07.07.14 Use bl prefix for functions By: CTP
   V1.2   26.02.15 Added blReadDisulphidesWholePDB()

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP File IO
   #FUNCTION  blReadDisulphidesPDB()
   Searches a PDB file for SSBOND records and constructs a linked list
   of information from these records.

   #FUNCTION blReadDisulphidesWholePDB()
   Builds a linked list of SSBOND record data from the header stored in
   a WHOLEPDB structure

*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>

#include "SysDefs.h"

#include "pdb.h"
#include "macros.h"
#include "fsscanf.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>DISULPHIDE *blReadDisulphidesPDB(FILE *fp, BOOL *error)
   -------------------------------------------------------
*//**

   \param[in]     *fp     PDB file pointer
   \param[out]    *error  Success
   \return                Linked list of disulphide information.
                          NULL if none found or error (Check flag)

   Searches a PDB file for SSBOND records and constructs a linked list
   of information from these records.
   Returns NULL if no disulphide information found. If memory allocation
   fails, the DISULPHIDE linked list formed thus far is returned and 
   the error flag is set to TRUE

-  14.10.93 Original   By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
DISULPHIDE *blReadDisulphidesPDB(FILE *fp, BOOL *error)
{
   DISULPHIDE *dis = NULL,
              *p   = NULL;
   char       buffer[MAXBUFF];

   *error = FALSE;

   while(fgets(buffer,MAXBUFF,fp))
   {
      /* Exit as soon as we reach an ATOM record                        */
      if(!strncmp(buffer,"ATOM  ",6)) break;
      if(!strncmp(buffer,"SSBOND",6))
      {
         /* Allocate memory in linked list                              */
         if(dis==NULL)
         {
            INIT(dis,DISULPHIDE);
            p=dis;
         }
         else
         {
            ALLOCNEXT(p,DISULPHIDE);
         }

         if(p==NULL)
         {
            *error = TRUE;
            return(dis);
         }

         /* Read data out of SSBOND record                              */
         fsscanf(buffer,"%15x%1s%5d%1s%7x%1s%5d%1s",
                 p->chain1,&p->res1,p->insert1,
                 p->chain2,&p->res2,p->insert2);
      }
   }

   return(dis);
}


/************************************************************************/
/*>DISULPHIDE *blReadDisulphidesWholePDB(WHOLEPDB *wpdb, BOOL *error)
   ------------------------------------------------------------------
*//**

   \param[in]     *wpdb   Whole PDB structure
   \param[out]    *error  Success
   \return                Linked list of disulphide information.
                          NULL if none found or error (Check flag)

   Searches the headers in a WHOLEPDB structure for SSBOND records and 
   constructs a linked list of information from these records.
   Returns NULL if no disulphide information found. If memory allocation
   fails, the DISULPHIDE linked list formed thus far is returned and 
   the error flag is set to TRUE

-  26.02.15 Original based on blReadDisulphidesPDB()   By: ACRM
*/
DISULPHIDE *blReadDisulphidesWholePDB(WHOLEPDB *wpdb, BOOL *error)
{
   DISULPHIDE *dis = NULL,
              *p   = NULL;
   STRINGLIST *s   = NULL;

   *error = FALSE;

   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      /* Exit as soon as we reach an ATOM record                        */
      if(!strncmp(s->string,"SSBOND",6))
      {
         /* Allocate memory in linked list                              */
         if(dis==NULL)
         {
            INIT(dis,DISULPHIDE);
            p=dis;
         }
         else
         {
            ALLOCNEXT(p,DISULPHIDE);
         }

         if(p==NULL)
         {
            *error = TRUE;
            return(dis);
         }

         /* Read data out of SSBOND record                              */
         fsscanf(s->string,"%15x%1s%5d%1s%7x%1s%5d%1s",
                 p->chain1,&p->res1,p->insert1,
                 p->chain2,&p->res2,p->insert2);
      }
   }

   return(dis);
}


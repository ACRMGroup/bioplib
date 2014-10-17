/************************************************************************/
/**

   \file       ReadSimplePIR.c
   
   \version    V2.9
   \date       07.07.14
   \brief      
   
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


**************************************************************************

   Usage:
   ======

\code
   int ReadSimplePIR(FILE *fp, int maxres, char **seqs)
\endcode

   This version (previously called ReadPIR()) is maintained only for
   compatibility with the old version. It will only read minimal 
   specification PIR files.

**************************************************************************

   Revision History:
   =================
-  V1.0  01.06.92 Original
-  V2.0  08.03.94 Changed name of ReadPIR() to ReadSimplePIR()
                  Added new ReadPIR().
-  V2.1  18.03.94 getc() -> fgetc()
-  V2.2  11.05.94 Changes to ReadPIR() for better compatibility with
                  PIR V38.0 and V39.0
-  V2.3  28.02.95 Added ReadRawPIR()
-  V2.4  13.03.95 Fixed bug in reading text lines in ReadRawPIR()
-  V2.5  26.07.95 Removed unused variables
-  V2.6  30.10.95 Cosmetic
-  V2.7  06.02.96 Removes trailing spaces from comment line
-  V2.8  18.06.02 Added string.h
-  V2.9  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP File IO
   #ROUTINE  blReadSimplePIR()
   Read a PIR file containing multiple chains of up to maxres amino acids.
   Doesn't handle special PIR characters
*/
/************************************************************************/
/* Includes
*/
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
/*>int blReadSimplePIR(FILE *fp, int maxres, char **seqs)
   ------------------------------------------------------
*//**

   \param[in]     *fp       File pointer
   \param[in]     maxres    Max number of residues in chain.
   \param[out]    **seqs    Array of pointers to sequences
   \return                  Number of chains. 0 if error

   Read a PIR file containing multiple chains of up to maxres amino acids.
   Each chain is returned in seqs[].
   The number of chains is returned by the routine.
   0 is returned if a memory allocation failed
   
-  01.06.91 Original
-  03.03.94 Added check on case before toupper(). Changed name.
-  18.03.94 Changed getc() to fgetc()
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blReadSimplePIR(FILE *fp,
                    int  maxres,
                    char **seqs)
{
   char *buffer;
   int  rescount = 0,
        chain    = 0;
   
   /* Allocate space for the sequence                                   */
   buffer = (char *)malloc((maxres+1) * sizeof(char));
   if(!buffer) return(0);
   
   /* Read header lines from the file                                   */
   fgets(buffer,maxres-1,fp);
   fgets(buffer,maxres-1,fp);

   /* Now loop through to get the sequence                              */
   while(rescount<maxres && !feof(fp))
   {
      int ch;
      
      /* Get a character                                                */
      ch = fgetc(fp);
      if(ch==EOF) break;
      buffer[rescount] = ch;
      
      if(isalpha(buffer[rescount]))
      {
         /* If it's an alpha character, then toupper() it and 
            increment the counter.
         */
         buffer[rescount] = (isupper(buffer[rescount]) ? 
                             buffer[rescount]          :
                             toupper(buffer[rescount]));
         rescount++;
      }
      else if(buffer[rescount] == '*')
      {
         /* If it's a star, then it's the end of a chain, 
            so copy the chain 
         */
         buffer[rescount] = '\0';
         seqs[chain] = (char *)malloc((rescount+2)*sizeof(char));
         if(!seqs[chain]) return(-1);
         strcpy(seqs[chain],buffer);
         chain++;
         rescount=0;
      }
   }
   
   /* Check to see if the last chain ended without a *                  */
   if(rescount)
   {
      buffer[rescount] = '\0';
      seqs[chain] = (char *)malloc((rescount+2)*sizeof(char));
      if(!seqs[chain])
      {
         int i;
         for(i=0; i<=chain; i++)
         {
            if(seqs[i]) free(seqs[i]);
         }
         return(0);
      }
      strcpy(seqs[chain],buffer);
      chain++;
   }
   free(buffer);

   return(chain);
}



/************************************************************************/
/**

   \file       DNAtoAA.c
   
   \version    V1.1
   \date       07.07.14
   \brief      Convert DNA codons to amino acid 1-letter code
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1994-2014
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


      #include "bioplib/seq.h" to define prototype

**************************************************************************

   Revision History:
   =================
   
-  V1.1  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling Sequence Data
   #SUBGROUP Conversions
   #ROUTINE  blDNAtoAA()
   Converts a nucleic acid codon to the 1-letter amino acid equivalent.
   Termination codons are returned as X. No special action is taken
   for initiation codons.
*/
/************************************************************************/
/* Includes
*/
#include "macros.h"
#include <string.h>
#include <ctype.h>

/************************************************************************/
/* Defines and macros
*/
#define NUCINDEX(x) ((x)=='T' ? (0) : \
                    ((x)=='U' ? (0) : \
                    ((x)=='C' ? (1) : \
                    ((x)=='A' ? (2) : \
                     (3)))))

/************************************************************************/
/* Globals
*/
static char *sAACode[4][4] =
{  { "FFLL", "SSSS", "YYXX", "CCXW" },   /* TTX, TCX, TAX, TGX          */
   { "LLLL", "PPPP", "HHQQ", "RRRR" },   /* CTX, CCX, CAX, CGX          */
   { "IIIM", "TTTT", "NNKK", "SSRR" },   /* ATX, ACX, AAX, AGX          */
   { "VVVV", "AAAA", "DDEE", "GGGG" }    /* GTX, GCX, GAX, GGX          */
} ;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>char blDNAtoAA(char *dna)
   -------------------------
*//**

   \param[in]     *dna        DNA/RNA codon
   \return                    1-letter amino acid code (X=termination)

   Converts a nucleic acid codon to the 1-letter amino acid equivalent.
   Termination codons are returned as X. No special action is taken
   for initiation codons.

-  18.04.94 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
char blDNAtoAA(char *dna)
{
   char buffer[8], *p;
   int idx1, idx2, idx3;
   
   KILLLEADSPACES(p,dna);
   
   strncpy(buffer,p,8);
   buffer[7] = '\0';
   UPPER(buffer);

   idx1 = NUCINDEX(buffer[0]);
   idx2 = NUCINDEX(buffer[1]);
   idx3 = NUCINDEX(buffer[2]);
   
   return(sAACode[idx1][idx2][idx3]);
}


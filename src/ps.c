/************************************************************************/
/**

   \file       ps.c
   
   \version    V1.5
   \date       07.07.14
   \brief      PostScript plotting routines
   
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

   These routines supply rudimentry PostScript support with
   simple commands from C.

   This version does not support EPSF.

**************************************************************************

   Usage:
   ======
   You should #include "ps.h" in your program and assign
   values to the global variables PSxpicsize,PSypicsize,PSxoffset
   and PSyoffset. All these values are in inches. Your plot should
   then run between 0.0 and 1.0 (you need to look after scaling
   to fit within these boundaries yourself).

   Start with a call to PSInit which will set up the scaling and other
   routines within the PostScript program.

**************************************************************************

   Revision History:
   =================
-  V1.0  22.11.90 Original
-  V1.1  06.07.93 Modified for book
-  V1.2  27.07.93 Changed I/O precision to double
-  V1.4  22.06.94 The file pointer is now global rather than static
-  V1.5  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define _PS_MAIN

/************************************************************************/
/* Doxygen
   -------
   #GROUP    Graphics
   #SUBGROUP Postscript
   #FUNCTION blPSInit()
   Initialises the file writing the Prologue. The filename and creator
   are written into the Prologue and EPSFxoff and EPSFyoff are used to
   calculate the bounding box size for EPSF plots.

   #FUNCTION blPSThick()
   Set the line thickness

   #FUNCTION blPSMove()
   Move to X,Y

   #FUNCTION blPSDraw()
   Draw to X,Y

   #FUNCTION blPSSetDash()
   Set a line dash pattern which must be supplied as a string

   #FUNCTION blPSClearDash()
   Clear the dash pattern to a full line

   #FUNCTION blPSStroke()
   Actually draw what you've just done onto the paper

   #FUNCTION blPSFont()
   Set the font and size

   #FUNCTION blPSLText()
   Left justify text

   #FUNCTION blPSCBText()
   Centers a piece of text with X,Y being the Coords of the BOTTOM centre 
   point

   #FUNCTION blPSROffText()
   Right justify text with offset in device coordinates (points).

   #FUNCTION blPSLCText()
   Left justify text, centred on Y

   #FUNCTION blPSCTText()
   Centers a piece of text with X,Y being the Coords of the TOP centre 
   point

   #FUNCTION blPSVText()
   Write vertical text centred on x,y offset back along x by the size of
   label and by xoff in pts. Used, for example, to title the y-axis of
   a graph. The `label' specification is used to calculate an amount by
   which to move the text back. Typically this would be the longest data
   label on the graph's Y-axis.

   #FUNCTION blPSShowText()
   Displays text, processing it first if any control codes are found. Used
   by the various text positioning routines.

   #FUNCTION blPSEnd()
   End of page

   #FUNCTION blPSCorrectCase()
   Goes through a fontname and fixes case to match the required standard.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "MathType.h"
#include "SysDefs.h"

#include "plotting.h"
#include "general.h"
#include "ps.h"

/************************************************************************/
/* Globals
*/
static char sPSBuff[200];
static REAL sTextHeight = 10.0;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL blPSInit(char *FName, char *creator, char *AltFont)
   --------------------------------------------------------
*//**

   \param[in]     *FName     PostScript filename
   \param[in]     *creator   Creator string
   \param[in]     *AltFont   Alternate font (normall greek style)
   \return                   Success?

   Initialises the file writing the Prologue. The filename and creator
   are written into the Prologue and EPSFxoff and EPSFyoff are used to
   calculate the bounding box size for EPSF plots.

-  08.05.92 Added definitions of raise, lower and greek. New parameter
            to specify name of alternate font.
-  23.06.92 Padded BoundingBox with spaces and setting of PSXMin, etc.
-  26.07.92 Correctly divide EPSFxoff and EPSFyoff by 72 when writing
            offsets in prologue.
-  15.09.92 Added support for Amiga reencoding. Starting dimensions
            for BoundingBox set to dimensions set by Paper.
-  06.07.93 Opens the file
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blPSInit(char    *FName,
              char    *creator,
              char    *AltFont)
{
   struct tm   *time_struc;
   time_t      time_value;
   
   if((gPSFile = fopen(FName,"w")) == NULL) return(FALSE);
   
   /* Header information                                                */
   fputs("%!PS-Adobe-2.0\n",gPSFile);
   fputs("%%Creator: ",gPSFile);
   fputs(creator,gPSFile);
   fputs("\n",gPSFile);
   fputs("%%For: (Andrew Martin Support Library)\n",gPSFile);
   sprintf(sPSBuff,"%%%%Title: (%s)\n",FName);
   fputs(sPSBuff,gPSFile);
   time(&time_value);
   time_struc = localtime(&time_value);
   sprintf(sPSBuff,"%%%%CreationDate: %s",asctime(time_struc));
   fputs(sPSBuff,gPSFile);
   fputs("%%Pages: 1\n",gPSFile);
   fputs("%%DocumentFonts: atend\n",gPSFile);
   
   fputs("%%EndComments\n",gPSFile);
   
   fputs("%%---------------Prologue-----------------\n",gPSFile);
   
   sprintf(sPSBuff,"/xpicsize %7.4f def\n",PSxpicsize);
   fputs(sPSBuff,gPSFile);
   sprintf(sPSBuff,"/ypicsize %7.4f def\n",PSypicsize);
   fputs(sPSBuff,gPSFile);

   sprintf(sPSBuff,"/xoffset %7.4f def\n",PSxoffset);
   fputs(sPSBuff,gPSFile);
   sprintf(sPSBuff,"/yoffset %7.4f def\n",PSyoffset);
   fputs(sPSBuff,gPSFile);

   fputs("/xscale { xpicsize 72 mul mul } def\n",gPSFile);
   fputs("/yscale { ypicsize 72 mul mul } def\n",gPSFile);
   fputs("/xunits { xscale xoffset 72 mul add } def\n",gPSFile);
   fputs("/yunits { yscale yoffset 72 mul add } def\n\n",gPSFile);

   /* A routine to set the font                                         */
   fputs("/font\n",gPSFile);
   fputs("{\n",gPSFile);
   fputs("   findfont exch scalefont setfont\n",gPSFile);
   fputs("}def\n",gPSFile);

   /* The RJust procedure                                               */
   fputs("/rightJustifyText\n",gPSFile);
   fputs("{  /RightColumn exch def\n",gPSFile);
   fputs("   dup\n",gPSFile);
   fputs("   stringwidth pop\n",gPSFile);
   fputs("   RightColumn exch sub\n",gPSFile);
   fputs("   Line moveto\n",gPSFile);
   fputs("   show\n",gPSFile);
   fputs("}  def\n",gPSFile);

   /* The circle procedure                                              */
   fputs("/circle\n",gPSFile);
   fputs("{  0 360 arc\n",gPSFile);
   fputs("}  def\n",gPSFile);
   
   /* The raise procedure                                               */
   fputs("/raise\n",gPSFile);
   fputs("{  sTextHeight 1.5 div FontName font\n",gPSFile);
   fputs("   0 sTextHeight 2 div     rmoveto\n",gPSFile);
   fputs("   show\n",gPSFile);
   fputs("   0 sTextHeight 2 div neg rmoveto\n",gPSFile);
   fputs("   sTextHeight FontName font\n",gPSFile);
   fputs("}  def\n",gPSFile);
   
   /* The lower procedure                                               */
   fputs("/lower\n",gPSFile);
   fputs("{  sTextHeight 1.5 div FontName font\n",gPSFile);
   fputs("   0 sTextHeight 4 div neg rmoveto\n",gPSFile);
   fputs("   show\n",gPSFile);
   fputs("   0 sTextHeight 4 div     rmoveto\n",gPSFile);
   fputs("   sTextHeight FontName font\n",gPSFile);
   fputs("}  def\n",gPSFile);
   
   /* The greek procedure                                               */
   fputs("/greek\n",gPSFile);
   fputs("{  AltFontName findfont sTextHeight scalefont setfont\n",
         gPSFile);
   fputs("   show\n",gPSFile);
   fputs("   sTextHeight FontName font\n",gPSFile);
   fputs("}  def\n",gPSFile);
   
   /* Define the alternate font for Greek                               */
   sprintf(sPSBuff,"/AltFontName  /%s def\n",AltFont);
   fputs(sPSBuff,gPSFile);
   
   /* The max procedure                                                 */
   fputs("/max\n",gPSFile);
   fputs("{  2 copy\n",gPSFile);
   fputs("   lt {exch} if\n",gPSFile);
   fputs("   pop\n",gPSFile);
   fputs("}  def\n",gPSFile);

   fputs("%%EndProlog\n",gPSFile);
   fputs("%%Page 1 1\n",gPSFile);
   fputs("%%---------------Script-----------------\n",gPSFile);
   
   return(TRUE);
}

/************************************************************************/
/*>void blPSThick(REAL thickness)
   ------------------------------
*//**

   \param[in]     thickness     Line thickness

   Set the line thickness

-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSThick(REAL thickness)
{
   sprintf(sPSBuff,"%f setlinewidth\n",thickness);
   fputs(sPSBuff,gPSFile);
}

/************************************************************************/
/*>void blPSMove(REAL X, REAL Y)
   -----------------------------
*//**

   \param[in]     X     X coordinate
   \param[in]     Y     Y coordinate

   Move to X,Y

-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSMove(REAL X,
              REAL Y)
{
   sprintf(sPSBuff,"%7.4f xunits %7.4f yunits moveto\n",X,Y);
   fputs(sPSBuff,gPSFile);
}

/************************************************************************/
/*>void blPSDraw(REAL X, REAL Y)
   -----------------------------
*//**

   \param[in]     X     X coordinate
   \param[in]     Y     Y coordinate

   Draw to X,Y

-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSDraw(REAL X,
              REAL Y)
{
   sprintf(sPSBuff,"%7.4f xunits %7.4f yunits lineto\n",X,Y);
   fputs(sPSBuff,gPSFile);
}

/************************************************************************/
/*>void blPSSetDash(char *linepatt)
   --------------------------------
*//**

   \param[in]     *linepatt    Line pattern (a string of numbers)

   Set a line dash pattern which must be supplied as a string

-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSSetDash(char *linepatt)
{
   sprintf(sPSBuff,"[%s] 0 setdash\n",linepatt);
   fputs(sPSBuff,gPSFile);
}

/************************************************************************/
/*>void blPSClearDash(void)
   ------------------------
*//**

   Clear the dash pattern to a full line

-  30.06.92 Removed redundant sprintf()'s
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSClearDash(void)
{
   fputs("[] 0 setdash\n",gPSFile);
}

/************************************************************************/
/*>void blPSStroke(void)
   ---------------------
*//**

   Actually draw what you've just done onto the paper

-  30.06.92 Removed redundant sprintf()'s
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSStroke(void)
{
   fputs("stroke\n",gPSFile);
}

/************************************************************************/
/*>void blPSFont(char *fontname, REAL size)
   ----------------------------------------
*//**

   \param[in]     *fontname    Font to set
   \param[in]     size         Point size of font

   Set the font and size
-  08.05.92 Changed to support raise and lower
-  23.06.92 Set sTextHeight
-  15.09.92 Changed to support Amiga reencoding
-  06.07.93 Removed Amiga reecncoding for general distribution
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSFont(char  *fontname,
              REAL  size)
{
   char font[80];
   
   blStringToUpper(fontname,font);

   if(!strncmp(font,"TIMES",5))
   {
      if(font[5] != '-')
         strcpy(fontname+5,"-Roman\0");
   }
   
   sprintf(sPSBuff,"/FontName /%s def\n",blPSCorrectCase(fontname));
   fputs(sPSBuff,gPSFile);
   sprintf(sPSBuff,"/sTextHeight %f def\n",size);
   fputs(sPSBuff,gPSFile);
   
   fputs("sTextHeight FontName font\n",gPSFile);

   sTextHeight = size;
}

/************************************************************************/
/*>void blPSLText(REAL X, REAL Y, char *label)
   -------------------------------------------
*//**

   \param[in]     X       X coordinate
   \param[in]     Y       Y coordinate
   \param[in]     *label  Text to be printed

   Left justify text

-  30.06.92 Removed redundant sprintf()'s
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSLText(REAL X,
               REAL Y,
               char *label)
{
   /* Define the current line and position                              */
   sprintf(sPSBuff,"%7.4f xunits\n",X);
   fputs(sPSBuff,gPSFile);
   sprintf(sPSBuff,"/Line %7.4f yunits def\n",Y);
   fputs(sPSBuff,gPSFile);
   fputs("Line moveto\n",gPSFile);
   /* Put the text into the file and display                            */
   blPSShowText(label);
}
   
/************************************************************************/
/*>void blPSCBText(REAL X, REAL Y, REAL Offset, char *label)
   ---------------------------------------------------------
*//**

   \param[in]     X       X coordinate
   \param[in]     Y       Y coordinate
   \param[in]     Offset  Y offset - multiple of font height. Moves up
                          by this quantity
   \param[in]     *label  Text to be printed

   Centers a piece of text with X,Y being the Coords of the BOTTOM centre 
   point

-  30.06.92 Removed redundant sprintf()'s
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSCBText(REAL X,
                REAL Y,
                REAL Offset,
                char *label)
{
   char  *buffer;
   
   buffer = label;
   while(*buffer == ' ') buffer++;
   
   /* Put the string on the stack                                       */
   sprintf(sPSBuff,"(%s)\n",blSimplifyText(buffer));
   fputs(sPSBuff,gPSFile);
   /* Define the current line                                           */
   sprintf(sPSBuff,
           "%7.4f yunits 2 sub /Line exch def\n",Y);
   fputs(sPSBuff,gPSFile);
   /* We are now left with label on the stack.
      Find half its width                                               */
   fputs("dup stringwidth pop 2 div\n",gPSFile);
   /* Put the X position on the stack, subtract width/2                 */
   sprintf(sPSBuff,"%7.4f xunits exch sub\n",X);
   fputs(sPSBuff,gPSFile);
   /* We are left with the label and correct X on the stack
      Put Y on, move there and show the text                            */
   if (Offset == 0.0)
   {
      fputs("Line moveto pop\n",gPSFile);
   }
   else
   {
      sprintf(sPSBuff,"Line sTextHeight %7.4f mul sub moveto pop\n",
              Offset);
      fputs(sPSBuff,gPSFile);
   }
   blPSShowText(buffer);
}

/************************************************************************/
/*>void blPSROffText(REAL X, REAL Y, REAL offset, char *label)
   -----------------------------------------------------------
*//**

   \param[in]     X       X coordinate
   \param[in]     Y       Y coordinate
   \param[in]     offset  X offset in points; text moved to the left by
                          this amount
   \param[in]     *label  Text to be printed

   Right justify text with offset in device coordinates (points).

-  07.05.92 Original
-  30.06.92 Removed use of PSShowText(). This can't be used with the
            current rightJustifyText
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSROffText(REAL X,
                  REAL Y,
                  REAL offset,
                  char *label)
{
   /* Define the current line                                           */
   sprintf(sPSBuff,"/Line %7.4f yunits sTextHeight 2 div sub 2 \
add def\n",Y);
   fputs(sPSBuff,gPSFile);
   /* Put the text into the file                                        */
   sprintf(sPSBuff,"(%s)\n",label);
   fputs(sPSBuff,gPSFile);
   /* Ask to Right Justify                                              */
   sprintf(sPSBuff,"%7.4f xunits 4 sub %f add rightJustifyText\n",
           X,offset);
   fputs(sPSBuff,gPSFile);
}

/************************************************************************/
/*>void blPSLCText(REAL X, REAL Y, char *label)
   --------------------------------------------
*//**

   \param[in]     X       X coordinate
   \param[in]     Y       Y coordinate
   \param[in]     *label  Text to be printed

   Left justify text, centred on Y

-  08.05.92 Original
-  30.06.92 Removed redundant sprintf()'s
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSLCText(REAL X,
                REAL Y,
                char *label)
{
   /* Define the current line and position                              */
   sprintf(sPSBuff,"%7.4f xunits\n",X);
   fputs(sPSBuff,gPSFile);
   sprintf(sPSBuff,"/Line %7.4f yunits sTextHeight 2 div sub 2 \
add def\n",Y);
   fputs(sPSBuff,gPSFile);
   fputs("Line moveto\n",gPSFile);
   /* Put the text into the file and display                            */
   blPSShowText(label);
}

/************************************************************************/
/*>void blPSCTText(REAL X, REAL Y, REAL Offset, char *label)
   ---------------------------------------------------------
*//**

   \param[in]     X       X coordinate
   \param[in]     Y       Y coordinate
   \param[in]     Offset  Y offset in points. Moves down by this quantity
   \param[in]     *label  Text to be printed

   Centers a piece of text with X,Y being the Coords of the TOP centre 
   point

-  26.06.92 Changed strlen() to use SimplifyText()
-  30.06.92 Removed redundant sprintf()'s
-  01.07.92 Changed Offset to be in pts rather than a multiplier of 
            font size.
-  27.07.92 Changed update limits to account better for descenders.
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSCTText(REAL X,
                REAL Y,
                REAL Offset,
                char *label)
{
   char *buffer;
   
   buffer = label;
   while(*buffer == ' ') buffer++;
   
   /* Put the string on the stack                                       */
   sprintf(sPSBuff,"(%s)\n",blSimplifyText(buffer));
   fputs(sPSBuff,gPSFile);
   /* Define the current line                                           */
   sprintf(sPSBuff,
           "sTextHeight %7.4f yunits exch sub /Line exch def\n",Y);
   fputs(sPSBuff,gPSFile);
   /* We are now left with label on the stack.
      Find half its width                                               */
   sprintf(sPSBuff,"dup stringwidth pop 2 div\n");
   fputs(sPSBuff,gPSFile);
   /* Put the X position on the stack, subtract width/2                 */
   sprintf(sPSBuff,"%7.4f xunits exch sub\n",X);
   fputs(sPSBuff,gPSFile);
   /* We are left with the label and correct X on the stack
      Put Y on, move there and show the text                            */
   if (Offset == 0.0)
   {
      fputs("Line moveto pop\n",gPSFile);
   }
   else
   {
      sprintf(sPSBuff,"Line %7.4f add moveto pop\n",Offset);
      fputs(sPSBuff,gPSFile);
   }
   blPSShowText(buffer);
}


/************************************************************************/
/*>void blPSVText(REAL x,       REAL y,        REAL xoff,
                  char *text,   char *font,    REAL size,
                  char *label,  char *lfont,   REAL lsize)
   -------------------------------------------------------
*//**

   \param[in]     x        X coordinate (in data units)
   \param[in]     y        Y coordinate (in data units)
   \param[in]     xoff     X-offset in pts
   \param[in]     *text    Text to be written
   \param[in]     *font    Font in which to write it
   \param[in]     size     Size of font
   \param[in]     *label   Label to be used to calc x offset
   \param[in]     *lfont   Font of this label
   \param[in]     lsize    Size of this label

   Write vertical text centred on x,y offset back along x by the size of
   label and by xoff in pts. Used, for example, to title the y-axis of
   a graph. The `label' specification is used to calculate an amount by
   which to move the text back. Typically this would be the longest data
   label on the graph's Y-axis.

-  08.05.92 Original
-  30.06.92 Removed redundant sprintf()'s
-  15.09.92 Multiply string width by 0.65
-  27.07.93 Floating point precision -> double
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSVText(REAL x,          /* Data coord position (to be offset)     */
               REAL y,          /* Data coord position                    */
               REAL xoff,       /* X-offset in pts                        */
               char *text,      /* Text to be written                     */
               char *font,      /* Font in which to write it              */
               REAL size,       /* Size of font                           */
               char *label,     /* Label to be used to calc x offset      */
               char *lfont,     /* Font of this label                     */
               REAL lsize)      /* Size of this label                     */
{
   /* Put the text on the stack                                         */
   sprintf(sPSBuff,"(%s) dup\n",blSimplifyText(text));
   fputs(sPSBuff,gPSFile);

   /* Find the length of the string/2                                   */
   fputs("stringwidth pop 2 div\n",gPSFile);

   /* Do specified y-pos minus strlen/2                                 */
   sprintf(sPSBuff,"%7.4g yunits exch sub\n",y);
   fputs(sPSBuff,gPSFile);

   /* Set font to the offset label font                                 */
   blPSFont(lfont, lsize);
   
   /* Calculate the x-offset                                            */
   sprintf(sPSBuff,"%7.4g xunits (%s) stringwidth pop sub 5 sub %f \
add exch moveto\n",x,label,xoff);
   fputs(sPSBuff,gPSFile);

   /* Set font back                                                     */
   blPSFont(font, size);
   
   /* Display the actual text                                           */
   fputs("pop 90 rotate ",gPSFile);
   blPSShowText(text);
   fputs(" -90 rotate\n",gPSFile);
}

/************************************************************************/
/*>void blPSShowText(char *text)
   -----------------------------
*//**

   \param[in]     *text    Text to be written

   Displays text, processing it first if any control codes are found. Used
   by the various text positioning routines.

-  08.05.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSShowText(char *text)
{
   char  buffer[MAXBUFF];
   int   i, j;
   
   /* Walk along the string                                             */
   for(i=0, j=0; i<strlen(text) && j<MAXBUFF-1; i++)
   {
      switch(text[i])
      {
      case '\\':           /* Should interpret next character as Greek  */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(sPSBuff,"(%s) show ",buffer);
            fputs(sPSBuff,gPSFile);
            j = 0;
         }
         /* Output the next character as Greek                          */
         sprintf(sPSBuff,"(%c) greek ",text[++i]);
         fputs(sPSBuff,gPSFile);
         break;
      case '^':            /* Should raise next character               */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(sPSBuff,"(%s) show ",buffer);
            fputs(sPSBuff,gPSFile);
            j = 0;
         }
         /* If necessary build string                                   */
         if(text[++i] == '{')
            while(text[++i] != '}' && text[i] != '\0' && j<MAXBUFF-1)
               buffer[j++] = text[i];
         else
            buffer[j++] = text[i];
         /* Output raised string                                        */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(sPSBuff,"(%s) raise ",buffer);
            fputs(sPSBuff,gPSFile);
            j = 0;
         }
         break;
      case '_':            /* Should lower next character               */
         /* Finish the current string                                   */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(sPSBuff,"(%s) show ",buffer);
            fputs(sPSBuff,gPSFile);
            j = 0;
         }
         /* If necessary build string                                   */
         if(text[++i] == '{')
            while(text[++i] != '}' && text[i] != '\0' && j<MAXBUFF-1)
               buffer[j++] = text[i];
         else
            buffer[j++] = text[i];
         /* Output lowered string                                       */
         if(j)
         {
            buffer[j] = '\0';
            sprintf(sPSBuff,"(%s) lower ",buffer);
            fputs(sPSBuff,gPSFile);
            j = 0;
         }
         break;
      case '(':         /* Need to insert a \, before falling through   */
      case ')':
         buffer[j++] = '\\';
      default:          /* An ordinary character                        */
         buffer[j++] = text[i];
         break;
      }
   }
   
   if(j)
   {
      buffer[j] = '\0';
      sprintf(sPSBuff,"(%s) show ",buffer);
      fputs(sPSBuff,gPSFile);
      j = 0;
   }

   if(strlen(text)) fputs("\n",gPSFile);
}

/************************************************************************/
/*>void blPSEnd(void)
   ------------------
*//**

   End of page

-  08.05.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blPSEnd(void)
{
   fputs("showpage\n",gPSFile);
   fputs("%%Trailer\n",gPSFile);
   
   fclose(gPSFile);
}

/************************************************************************/
/*>char *blPSCorrectCase(char *font)
   ---------------------------------
*//**

   \param[in]     *font    Input fontname
   \return                 Case-fixed fontname

   Goes through a fontname and fixes case to match the required standard.

-  08.05.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blPSCorrectCase(char *font)
{
   int i;
   /* Set the first character to UC                                     */
   font[0] = toupper(font[0]);
   /* Set everything else to lower                                      */
   for(i=1;i<strlen(font);i++)
      font[i] = tolower(font[i]);
   /* Now step through and upper the bits which need to be              */
   for(i=1;i<strlen(font);i++)
   {
      /* Anything after a -                                             */
      if(font[i]=='-')
         font[i+1] = toupper(font[i+1]);
      /* Start of the word oblique                                      */
      if(!strncmp(font+i,"oblique",7))
         font[i] = toupper(font[i]);
      /* Start of the word italic                                       */
      if(!strncmp(font+i,"italic",6))
         font[i] = toupper(font[i]);
      /* Start of the word roman                                        */
      if(!strncmp(font+i,"roman",5))
         font[i] = toupper(font[i]);
      /* Start of the word bold                                         */
      if(!strncmp(font+i,"bold",4))
         font[i] = toupper(font[i]);
   }

   return(font);
}


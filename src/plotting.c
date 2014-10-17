/************************************************************************/
/**

   \file       plotting.c
   
   \version    V1.3
   \date       07.07.14
   \brief      Top level HPGL/PS plotting routines
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1992-2014
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

   These functions provide a common interface to either PostScript or
   HPGL output.

   They simplified from a set written for the Amiga which also supports
   IFF-DR2D and Amiga screen output

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  06.04.92 Original    By: ACRM
-  V1.1  01.03.94 First release
-  V1.2  27.02.98 Removed unreachable breaks from switch() statement
-  V1.3  07.07.14 Use bl prefix for functions By: CTP

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Graphics
   #SUBGROUP Plotting
   #ROUTINE  blAMInitPlot()
   Initialise a device ready for plotting.

   #ROUTINE blAMSetPen()
   Change pen

   #ROUTINE blAMMove()
   Move to a position specified in data coordinates.

   #ROUTINE blAMDraw()
   Draw to a position specified in data coordinates.

   #ROUTINE blAMSetLineStyle()
   Set the line style

   #ROUTINE blAMEndLine()
   End a line; required by PostScript actually to draw on the paper.

   #ROUTINE blAMSetFont()
   Sets the current font using PostScript font names. If producing HPGL
   output, a lookup table is used to translate this to an HPGL font
   number

   #ROUTINE blAMText()
   Left/bottom justify text at position in data coordinates

   #ROUTINE blAMCBText()
   Centre-bottom justify text

   #ROUTINE blAMRText()
   Right/centre justify text at position in data coordinates; offset is 
   an x-offset specified in device coordinates (pt)

   #ROUTINE blAMLCText()
   Left/centre height justify text at position in data coordinates

   #ROUTINE blAMCTText()
   Centre/top justify text at position in data coordinates. 

   #ROUTINE blAMEndPlot()
   Close up a device after plotting.

   #ROUTINE blPS2HPGLFont()
   Takes the PostScript font name and works out the best HPGL equivalent
   from a translation table. On the first call, the table is read from 
   disk and space is allocated for it. If the routine is called with a 
   NULL parameter, the space allocated for the table is freed. It is 
   quite safe to call the routine again after this has occurred; the 
   table will simply be re-read from disk.

   #ROUTINE blSimplifyText()
   Removes control codes from a string for screen display. Also used for
   calculating string length. The returned string is stored as static
   within the routine
*/
/************************************************************************/
/* Includes
*/
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "general.h"
#include "macros.h"
#include "plotting.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXTRANS     10
#define TRANS_TABLE  "HPGL.ftr"

#define MAXPEN       6
#define MAXBUFF      160

/************************************************************************/
/* Globals
*/
static struct
{
   REAL     xmin,
            ymin,
            XPScale,
            YPScale;
}  sGraph;

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>BOOL blAMInitPlot(char *filename,  char *title,   int dest, 
                     REAL OutXSize,   REAL OutYSize, 
                     REAL OutXOff,    REAL OutYOff,
                     char *AltFont,   REAL xmargin,  REAL ymargin,
                     REAL DataXMin,   REAL DataYMin, 
                     REAL DataXMax,   REAL DataYMax)
   -------------------------------------------------------------
*//**
   \param[in]     *filename   File to open
   \param[in]     *title      Title for plot
   \param[in]     dest        Destination (DEST_PS or DEST_HPGL)
   \param[in]     OutXSize    Output plot X size (inches)
   \param[in]     OutYSize    Output plot Y size (inches)
   \param[in]     OutXOff     Output plot X offset (inches)
   \param[in]     OutYOff     Output plot Y offset (inches)
   \param[in]     *AltFont    Alternate font name 
   \param[in]     xmargin     Unprintable x margin (inches, HPGL)
                              Sensible default: 0.58
   \param[in]     ymargin     Unprintable y margin (inches, HPGL)
                              Sensible default: 0.1465
   \param[in]     DataXMin    Min data X value
   \param[in]     DataYMin    Min data Y value
   \param[in]     DataXMax    Max data X value
   \param[in]     DataYMax    Max data Y value
   \return                      TRUE: OK, FALSE: Failed

   Initialise a device ready for plotting.
   
-  07.05.92 Original
-  25.06.92 Added HPGL support. Moved setting of PS globals to here.
-  02.07.92 Put ClearWindow() in for screen plotting. Added seek to
            start of file for PS and HPGL file plotting.
-  16.07.92 Added DR2D support
-  17.07.92 Corrected call to InstallDR2DFonts() *before* InitDR2D().
            Changed buffer to [40]. Added bounds calc'n for DR2D.
-  20.07.92 Added alternate font parameter to DR2DInit(). Added check on
            DR2DInit() return.
-  22.07.92 Consider x-axis labelling in finding max x. Corrected to
            consider precision for log axes
-  24.07.92 Added extras parameter to ftostr
-  27.07.92 Removed the specification of EPSF offsets since the
            PSFixBoundingBox() routine takes care of all this. Increased
            size of xmax border for DR2D plots.
-  06.07.93 Changed parameters
-  27.02.98 Removed unreachable breaks from switch() statement
-  07.07.14 Use bl prefix for functions By: CTP
*/
BOOL blAMInitPlot(char   *filename,
                  char   *title,
                  int    dest,
                  REAL   OutXSize, 
                  REAL   OutYSize, 
                  REAL   OutXOff, 
                  REAL   OutYOff,
                  char   *AltFont,
                  REAL   xmargin,
                  REAL   ymargin,
                  REAL   DataXMin,
                  REAL   DataYMin,
                  REAL   DataXMax,
                  REAL   DataYMax)
{
   PSxpicsize     = OutXSize;
   PSypicsize     = OutYSize;
   PSxoffset      = OutXOff;
   PSyoffset      = OutYOff;
   
   sGraph.xmin    = DataXMin;
   sGraph.ymin    = DataYMin;
   
   sGraph.XPScale = 1.0 / (DataXMax - DataXMin);
   sGraph.YPScale = 1.0 / (DataYMax - DataYMin);
   
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      return(blPSInit(filename, title, AltFont));
   case DEST_HPGL:
      return(blHPGLInit(filename, AltFont, xmargin, ymargin));
   default:
      break;
   }

   return(FALSE);
}

/************************************************************************/
/*>void blAMSetPen(int dest, int pen)
   ----------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     pen       Pen number

   Change pen
-  06.04.92 Handles screen
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/   
void blAMSetPen(int   dest,
                int   pen)
{
   static REAL pens[MAXPEN] = {0.5, 0.75, 1.0, 1.25, 1.5, 1.75};

   pen--;
   pen = pen % MAXPEN;

   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      blPSThick(pens[pen]);
      break;
   case DEST_HPGL:
      blHPGLPen(pen+1);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMMove(int dest, REAL x, REAL y)
   ---------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate

   Move to a position specified in data coordinates.
-  06.04.92 Handles screen
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMMove(int  dest,
              REAL x,
              REAL y)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSMove(x,y);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;
      blHPGLMove(x,y);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMDraw(int dest, REAL x, REAL y)
   ---------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate

   Draw to a position specified in data coordinates.

-  06.04.92 Handles screen
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMDraw(int  dest,
              REAL x,
              REAL y)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSDraw(x,y);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLDraw(x,y);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMSetLineStyle(int dest, int style)
   ------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     style     Style number (0--5)

   Set the line style

-  08.04.92 Framework
-  07.05.92 Original (screen & PS)
-  25.06.92 Added HPGL support. Removed static store of style.
-  05.07.92 Corrected line patterns. They don't have commas!
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMSetLineStyle(int   dest,
                      int   style)
{
   static char PSPattern[6][16] = {"",    /* PostScript line patterns   */
                                   "2",
                                   "4 1 2 1",
                                   "4",
                                   "4 3 2 2 2 3",
                                   "4 2 4 2 2 2"};
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      if(style == 0) blPSClearDash();
      else           blPSSetDash(PSPattern[style]);
      break;
   case DEST_HPGL:
      blHPGLSetDash(style);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMEndLine(int dest)
   --------------------------
*//**

   \param[in]     dest      Destination

   End a line; required by PostScript actually to draw on the paper.

-  07.05.92 Original
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMEndLine(int  dest)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      blPSStroke();
      break;
   case DEST_HPGL:
      break;
   default:
      break;
   }
}


/************************************************************************/
/*>void blAMSetFont(int dest, char *PSFontName, REAL FontSize)
   -----------------------------------------------------------
*//**

   \param[in]     dest         Destination
   \param[in]     *PSFontName  PostScript font name
   \param[in]     FontSize     Size (in points) of font

   Sets the current font using PostScript font names. If producing HPGL
   output, a lookup table is used to translate this to an HPGL font
   number

-  07.04.92 Framework
-  05.05.92 Original for Screen
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support.
-  29.06.92 Modified to use new PS2AmigaFont() for HPGL
-  13.07.92 Added DEF_FONT parameter to SetAmigaFont()
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMSetFont(int  dest, 
                 char *PSFontName,
                 REAL FontSize)
{
   int FontNum;
   
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      blPSFont(PSFontName, FontSize);
      break;
   case DEST_HPGL:
      FontNum = blPS2HPGLFont(PSFontName);
      blHPGLFont(FontNum, FontSize);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMText(int dest, REAL x, REAL y, char *text)
   ---------------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate
   \param[in]     *text     Text to write

   Left/bottom justify text at position in data coordinates

-  08.04.92 Handles screen
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMText(int  dest,
              REAL x,
              REAL y,
              char *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSLText(x,y,text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLLText(x,y,text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMCBText(int dest, REAL x, REAL y, char *text)
   -----------------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate
   \param[in]     *text     Text to print

   Centre-bottom justify text

-  07.04.92 Handles screen
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMCBText(int   dest,
                REAL  x,
                REAL  y,
                char  *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSCBText(x, y, 0.0, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLCBText(x, y, 0.0, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMRText(int dest, REAL x, REAL y, REAL offset, char *text)
   -----------------------------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate
   \param[in]     offset    Move left by this amount (in points)
   \param[in]     *text     Text to print

   Right/centre justify text at position in data coordinates; offset is 
   an x-offset specified in device coordinates (pt)
   
-  06.04.92 Handles screen
-  07.04.92 Fix to positioning
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support and offset
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMRText(int    dest,
               REAL   x,
               REAL   y,
               REAL   offset,
               char   *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSROffText(x, y, offset, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLROffText(x, y, offset, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMLCText(int dest, REAL x, REAL y, char *text)
   -----------------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate
   \param[in]     *text     Text to print

   Left/centre height justify text at position in data coordinates

-  08.04.92 Handles screen
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  06.05.92 Fix to height centering
-  07.05.92 Added PS support
-  08.05.92 Corrected Y-pos for PS
-  25.06.92 Added HPGL support
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMLCText(int    dest,
                REAL   x,
                REAL   y,
                char   *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSLCText(x,y,text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;
      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLLCText(x,y,text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text)
   --------------------------------------------------------------------
*//**

   \param[in]     dest      Destination
   \param[in]     x         X coordinate
   \param[in]     y         Y coordinate
   \param[in]     CTOffset  Move down by this amount (points)
   \param[in]     *text     Text to print

   Centre/top justify text at position in data coordinates. 

-  06.04.92 Handles screen
-  07.04.92 Fix to positioning
-  10.04.92 Added log support
-  29.04.92 Added check on log bounds
-  07.05.92 Added PS support
-  25.06.92 Added HPGL support
-  01.07.92 Changed for new versions of PSCTText() and HPGLCTText() which
            take offset in points.
-  16.07.92 Added DR2D support
-  07.06.93 Added CTOffset param
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMCTText(int   dest,
                REAL  x,
                REAL  y,
                REAL  CTOffset,
                char  *text)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
      
   case DEST_PS:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      blPSCTText(x, y, CTOffset, text);
      break;
   case DEST_HPGL:
      x = (x-sGraph.xmin) * sGraph.XPScale;

      y = (y-sGraph.ymin) * sGraph.YPScale;

      blHPGLCTText(x, y, CTOffset, text);
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>void blAMEndPlot(int dest)
   --------------------------
*//**

   \param[in]     dest      Destination

   Close up a device after plotting.
   
-  07.05.92 Original
-  25.06.92 Added HPGL support
-  01.07.92 Added blank WriteMessage() when plotting to screen
-  16.07.92 Added DR2D support
-  07.07.14 Use bl prefix for functions By: CTP
*/
void blAMEndPlot(int  dest)
{
   switch(dest)
   {
   case DEST_SCREEN:
      /* Screen Version                                                 */
      break;
   case DEST_PS:
      blPSEnd();
      break;
   case DEST_HPGL:
      blHPGLEnd();
      break;
   default:
      break;
   }
}

/************************************************************************/
/*>int blPS2HPGLFont(char *font)
   -----------------------------
*//**

   \param[in]     *font    PostScript font name
   \return                    HPGL font number

   Takes the PostScript font name and works out the best HPGL equivalent
   from a translation table. On the first call, the table is read from 
   disk and space is allocated for it. If the routine is called with a 
   NULL parameter, the space allocated for the table is freed. It is 
   quite safe to call the routine again after this has occurred; the 
   table will simply be re-read from disk.
   
   If the requested translation is unsuccessful, 0 will be returned as 
   the font number.
   
-  06.07.93 Original    By: ACRM
-  07.07.14 Use bl prefix for functions By: CTP
*/
int blPS2HPGLFont(char *font)
{
   int         i;
   char        buffer[MAXBUFF];
   static BOOL FirstCall   = TRUE;
   static int  NTrans      = 0;
   static struct                    /* Font translation table           */
   {
      char  *PSFont;
      int   HPGLFont;
   }  FontTable[MAXTRANS];
   
   
   /* If called with a NULL font name, we free the font table and 
      return.
   */
   if(font==NULL)
   {
      for(i=0; i<NTrans; i++)
      {
         if(FontTable[i].PSFont)
         {
            free(FontTable[i].PSFont);
            FontTable[i].PSFont = NULL;
         }
      }
      FirstCall = TRUE;
      NTrans    = 0;
      return(0);
   }
   
   
   /* If it's the first call we read the font translation table allocating
      space for the font names. Should we fail, we just drop out.
   */
   if(FirstCall)
   {
      FILE  *fp   = NULL;  /* Font translation file                     */
      
      FirstCall   = FALSE;
      
      if((fp=fopen(TRANS_TABLE,"r")) != NULL)   /* If found table       */
      {
         char  buffer[MAXBUFF];                 /* Buffer for file      */
         char  FontName[40];                    /* Font name read       */
         int   FontNum;                         /* HPGL number read     */
               
         while(fgets(buffer,MAXBUFF-1,fp))      /* Read file            */
         {
            TERMINATE(buffer);
            sscanf(buffer,"%s %d",FontName,&FontNum);
            
            if(strlen(FontName) && FontName[0] != '!')
            {
               /* We've got a font name and number. Allocate space      */
               FontTable[NTrans].PSFont = 
                  (char *)malloc((strlen(FontName)+1) * sizeof(char));
                  
               /* No room!                                              */
               if(FontTable[NTrans].PSFont == NULL) break;
               
               /* Copy in the info, down casing the font name & removing
                  leading spaces and tabs
               */
               blStringToLower(FontName, buffer);
               strcpy(FontTable[NTrans].PSFont, blKillLeadSpaces(buffer));
               FontTable[NTrans].HPGLFont = FontNum;
               
               /* Increment the translation count                       */
               if(++NTrans > MAXTRANS) break;
            }
         }
         fclose(fp);
      }
   }
   
   /* If we've got some translations, search for the specified font     */
   if(NTrans)
   {
      char  *ptr = NULL;
      
      blStringToLower(font, buffer);
      ptr = blKillLeadSpaces(buffer);

      for(i=0; i<NTrans; i++)
      {
         if(!strcmp(ptr, FontTable[i].PSFont))
         {
            return(FontTable[i].HPGLFont);
         }
      }
   }
   
   /* If no translations, or font not found, just return 0              */
   return(0);
}

/************************************************************************/
/*>char *blSimplifyText(char *string)
   ----------------------------------
*//**

   \param[in]     *string   String containing control codes
   \return                    String with control codes removed

   Removes control codes from a string for screen display. Also used for
   calculating string length. The returned string is stored as static
   within the routine

-  06.05.92 Original
-  07.07.14 Use bl prefix for functions By: CTP
*/
char *blSimplifyText(char *string)
{
   static char retstring[MAXBUFF];
   int         i, j;
   
   /* Just return the string unaltered if its too long                  */
   if(strlen(string) > MAXBUFF-1) return(string);
   
   /* Walk along the string                                             */
   for(i=0, j=0; i<strlen(string) && j<MAXBUFF-1; i++)
   {
      switch(string[i])
      {
      case '\\':           /* Should interpret next character as Greek  */
         retstring[j++] = string[++i];
         break;
      case '^':            /* Should raise next character               */
         if(string[++i] == '{')
            while(string[++i] != '}' && string[i] != '\0' && j<MAXBUFF-1)
               retstring[j++] = string[i];
         else
            retstring[j++] = string[i];
         break;
      case '_':            /* Should lower next character               */
         if(string[++i] == '{')
            while(string[++i] != '}' && string[i] != '\0' && j<MAXBUFF-1)
               retstring[j++] = string[i];
         else
            retstring[j++] = string[i];
         break;
      default:             /* An ordinary character                     */
         retstring[j++] = string[i];
         break;
      }
   }
   
   retstring[j] = '\0';

   return(retstring);
}


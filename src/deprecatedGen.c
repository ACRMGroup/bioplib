/************************************************************************/
/**

   \file       deprecatedGen.c
   
   \version    V1.2
   \date       14.08.14
   \brief      Source code for all deprecated functions.
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2014
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

   Allows code utilising BiopLib to compile without modification to use 
   function names with the 'bl' prefix.(The 'bl' prefix was introduced in 
   July 2014.) 


**************************************************************************

   Usage:
   ======

   Code using deprecated functions can be compiled without modification, 
   however, all deprecated functions use the DEPRECATED macro which prints
   a warning message to stderr when a deprecated function is called.

   The warning message can be suppressed as either an option at the 
   compilation of the library or by setting an environment variable at 
   runtime.

   \par Documentation
   Doxygen is not set to document deprecated functions. To document 
   deprecated functions, add "deprecated" to the ENABLED SECTIONS tag in
   the Doxygen config file:
   bioplib/doc/doxygen/Doxyfile

**************************************************************************

   Revision History:
   =================

-  V1.0  31.07.14 Original By: CTP
-  V1.1  08.08.14 Separated Biop and Gen deprecation  By: ACRM
-  V1.2  14.08.14 Removed unnecessary includes. 
                  Corrected safemem.h function names   By: CTP

*************************************************************************/
/* Includes
*/
#include "deprecated.h"

#include "general.h"
#include "BuffInp.h"
#include "ErrStack.h"
#include "MathUtil.h"
#include "WindIO.h"
#include "angle.h"
#include "array.h"
#include "help.h"
#include "hpgl.h"
#include "matrix.h"
#include "parse.h"
#include "plotting.h"
#include "ps.h"
#include "safemem.h"


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
/** \cond deprecated                                                    */
/************************************************************************/
/* Renamed functions: general.h                                         */


void StringToLower(char *string1, char *string2)
{
   DEPRECATED("StringToLower()","blStringToLower()");
   blStringToLower(string1, string2);
}

void StringToUpper(char *string1, char *string2)
{
   DEPRECATED("StringToUpper()","blStringToUpper()");
   blStringToUpper(string1, string2);
}

char *KillLeadSpaces(char *string)
{
   DEPRECATED("KillLeadSpaces()","blKillLeadSpaces()");
   return(blKillLeadSpaces(string));
}

void KillLine(FILE *fp)
{
   DEPRECATED("KillLine()","blKillLine()");
   blKillLine(fp);
}

void SetExtn(char *File, char *Ext)
{
   DEPRECATED("SetExtn()","blSetExtn()");
   blSetExtn(File, Ext);
}

int chindex(char *string, char ch)
{
   DEPRECATED("chindex()","blChindex()");
   return(blChindex(string, ch));
}

void Word(char *string1, char *string2)
{
   DEPRECATED("Word()","blWord()");
   blWord(string1, string2);
}

void WordN(char *string1, char *string2, int  MaxChar)
{
   DEPRECATED("WordN()","blWordN()");
   blWordN(string1, string2, MaxChar);
}

void padterm(char *string, int length)
{
   DEPRECATED("padterm()","blPadterm()");
   blPadterm(string, length);
}

void padchar(char *string, int length, char ch)
{
   DEPRECATED("padchar()","blPadchar()");
   blPadchar(string, length, ch);
}

BOOL CheckExtn(char *string, char *ext)
{
   DEPRECATED("CheckExtn()","blCheckExtn()");
   return(blCheckExtn(string, ext));
}

char *ftostr(char *str, int maxlen, REAL x, int precision)
{
   DEPRECATED("ftostr()","blFtostr()");
   return(blFtostr(str, maxlen, x, precision));
}

void GetFilestem(char *filename, char *stem)
{
   DEPRECATED("GetFilestem()","blGetFilestem()");
   blGetFilestem(filename, stem);
}

int upstrcmp(char *word1, char *word2)
{
   DEPRECATED("upstrcmp()","blUpstrcmp()");
   return(blUpstrcmp(word1, word2));
}

int upstrncmp(char *word1, char *word2, int ncomp)
{
   DEPRECATED("upstrncmp()","blUpstrncmp()");
   return(blUpstrncmp(word1, word2, ncomp));
}

char *GetWord(char *buffer, char *word, int maxsize)
{
   DEPRECATED("GetWord()","blGetWord()");
   return(blGetWord(buffer, word, maxsize));
}

BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
{
   DEPRECATED("OpenStdFiles()","blOpenStdFiles()");
   return(blOpenStdFiles(infile, outfile, in, out));
}

FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
{
   DEPRECATED("OpenFile()","blOpenFile()");
   return(blOpenFile(filename, envvar, mode, noenv));
}

int countchar(char *string, char ch)
{
   DEPRECATED("countchar()","blCountchar()");
   return(blCountchar(string, ch));
}

char *fgetsany(FILE *fp)
{
   DEPRECATED("fgetsany()","blFgetsany()");
   return(blFgetsany(fp));
}

char *strcatalloc(char *instr, char *catstr)
{
   DEPRECATED("strcatalloc()","blStrcatalloc()");
   return(blStrcatalloc(instr, catstr));
}

STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
{
   DEPRECATED("StoreString()","blStoreString()");
   return(blStoreString(StringList, string));
}

BOOL InStringList(STRINGLIST *StringList, char *string)
{
   DEPRECATED("InStringList()","blInStringList()");
   return(blInStringList(StringList, string));
}

void FreeStringList(STRINGLIST *StringList)
{
   DEPRECATED("FreeStringList()","blFreeStringList()");
   blFreeStringList(StringList);
}

char *QueryStrStr(char *string, char *substring)
{
   DEPRECATED("QueryStrStr()","blQueryStrStr()");
   return(blQueryStrStr(string, substring));
}

void IndexReal(REAL *arrin, int *indx, int n)
{
   DEPRECATED("IndexReal()","blIndexReal()");
   blIndexReal(arrin, indx, n);
}

FILE *OpenOrPipe(char *filename)
{
   DEPRECATED("OpenOrPipe()","blOpenOrPipe()");
   return(blOpenOrPipe(filename));
}

int CloseOrPipe(FILE *fp)
{
   DEPRECATED("CloseOrPipe()","blCloseOrPipe()");
   return(blCloseOrPipe(fp));
}

BOOL WrapString(char *in, char *out, int maxlen)
{
   DEPRECATED("WrapString()","blWrapString()");
   return(blWrapString(in, out, maxlen));
}

BOOL WrapPrint(FILE *out, char *string)
{
   DEPRECATED("WrapPrint()","blWrapPrint()");
   return(blWrapPrint(out, string));
}

void RightJustify(char *string)
{
   DEPRECATED("RightJustify()","blRightJustify()");
   blRightJustify(string);
}

char *GetWordNC(char *buffer, char *word, int maxlen)
{
   DEPRECATED("GetWordNC()","blGetWordNC()");
   return(blGetWordNC(buffer, word, maxlen));
}

void getfield(char *buffer, int start, int width, char *str)
{
   DEPRECATED("getfield()","blGetfield()");
   blGetfield(buffer, start, width, str);
}



/************************************************************************/
/* Renamed functions: BuffInp.h                                         */

INBUFFER *OpenBufferedFile(char *filename, int maxstr)
{
   DEPRECATED("OpenBufferedFile()","blOpenBufferedFile()");
   return(blOpenBufferedFile(filename, maxstr));
}

BOOL ReadBufferedFile(INBUFFER *bfp, char *string, int length)
{
   DEPRECATED("ReadBufferedFile()","blReadBufferedFile()");
   return(blReadBufferedFile(bfp, string, length));
}

BOOL ProbeBufferedFile(INBUFFER *bfp, char *string, int length)
{
   DEPRECATED("ProbeBufferedFile()","blProbeBufferedFile()");
   return(blProbeBufferedFile(bfp, string, length));
}



/************************************************************************/
/* Renamed functions: ErrStack.h                                        */

void StoreError(char *routine, char *error)
{
   DEPRECATED("StoreError()","blStoreError()");
   blStoreError(routine, error);
}

void ShowErrors(void *PrintRoutine(char *), BOOL Trace)
{
   DEPRECATED("ShowErrors()","blShowErrors()");
   blShowErrors(PrintRoutine, Trace);                     /* CHECK THIS */
}



/************************************************************************/
/* Renamed functions: MathUtil.h                                        */

void CalcSD(REAL val, int action, REAL *mean, REAL *SD)
{
   DEPRECATED("CalcSD()","blCalcSD()");
   blCalcSD(val, action, mean, SD);
}

void CalcExtSD(REAL val, int action, REAL *Sx, REAL *SxSq, int *NValues,
               REAL *mean, REAL *SD)
{
   DEPRECATED("CalcExtSD()","blCalcExtSD()");
   blCalcExtSD(val, action, Sx, SxSq, NValues, mean, SD);
}

REAL pearson(REAL *x, REAL *y, int NItem)
{
   DEPRECATED("pearson()","blPearson()");
   return(blPearson(x, y, NItem));
}

REAL pearson1(REAL *x, REAL *y, int NItem)
{
   DEPRECATED("pearson1()","blPearson1()");
   return(blPearson1(x, y, NItem));
}


void CrossProd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("CrossProd3()","blCrossProd3()");
   blCrossProd3(Out, In1, In2);
}

void VecSub3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("VecSub3()","blVecSub3()");
   blVecSub3(Out, In1, In2);
}

void VecAdd3(VEC3F *Out, VEC3F In1, VEC3F In2)
{
   DEPRECATED("VecAdd3()","blVecAdd3()");
   blVecAdd3(Out, In1, In2);
}

REAL VecLen3(VEC3F Vec)
{
   DEPRECATED("VecLen3()","blVecLen3()");
   return(blVecLen3(Vec));
}

REAL DistPtVect(VEC3F Point, VEC3F End1, VEC3F End2)
{
   DEPRECATED("DistPtVect()","blDistPtVect()");
   return(blDistPtVect(Point, End1, End2));
}


REAL PointLineDistance(REAL Px, REAL Py, REAL Pz,
                       REAL P1x, REAL P1y, REAL P1z,
                       REAL P2x, REAL P2y, REAL P2z,
                       REAL *Rx, REAL *Ry, REAL *Rz,
                       REAL *frac)
{
   DEPRECATED("PointLineDistance()","blPointLineDistance()");
   return(blPointLineDistance(Px,  Py,  Pz,
                              P1x, P1y, P1z,
                              P2x, P2y, P2z,
                              Rx,  Ry,  Rz,  frac));
}

ULONG factorial(int n)
{
   DEPRECATED("factorial()","blFactorial()");
   return(blFactorial(n));
}

ULONG factdiv(int n1, int n2)
{
   DEPRECATED("factdiv()","blFactdiv()");
   return(blFactdiv(n1, n2));
}

ULONG NPerm(int n, int r)
{
   DEPRECATED("NPerm()","blNPerm()");
   return(blNPerm(n, r));
}

ULONG NComb(int n, int r)
{
   DEPRECATED("NComb()","blNComb()");
   return(blNComb(n, r));
}



/************************************************************************/
/* Renamed functions: WindIO.h                                          */

void screen(char *string)
{
   DEPRECATED("Screen()","blScreen()");
   blScreen(string);
}

void prompt(char *string)
{
   DEPRECATED("Prompt()","blPrompt()");
   blPrompt(string);
}

void RePrompt(void)
{
   DEPRECATED("RePrompt()","blRePrompt()");
   blRePrompt();
}

void GetKybdString(char *string, int maxlen)
{
   DEPRECATED("GetKybdString()","blGetKybdString()");
   blGetKybdString(string, maxlen);
}

void PagingOn(void)
{
   DEPRECATED("PagingOn()","blPagingOn()");
   blPagingOn();
}

void PagingOff(void)
{
   DEPRECATED("PagingOff()","blPagingOff()");
   blPagingOff();
}

void WindowMode(BOOL mode)
{
   DEPRECATED("WindowMode()","blWindowMode()");
   blWindowMode(mode);
}

void WindowInteractive(BOOL mode)
{
   DEPRECATED("WindowInteractive()","blWindowInteractive()");
   blWindowInteractive(mode);
}

int YorN(char deflt)
{
   DEPRECATED("YorN()","blYorN()");
   return(blYorN(deflt));
}


/************************************************************************/
/* Renamed functions: angle.h                                           */

REAL angle(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj, REAL xk,
           REAL yk, REAL zk)
{
   DEPRECATED("Angle()","blAngle()");
   return(blAngle(xi, yi, zi, xj, yj, zj, xk, yk, zk));
}

REAL phi(REAL xi, REAL yi, REAL zi, REAL xj, REAL yj, REAL zj, REAL xk,
         REAL yk, REAL zk, REAL xl, REAL yl, REAL zl)
{
   DEPRECATED("Phi()","blPhi()");
   return(blPhi(xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl));
}

REAL simpleangle(REAL ang)
{
   DEPRECATED("Simpleangle()","blSimpleangle()");
   return(blSimpleangle(ang));
}

REAL TrueAngle(REAL opp, REAL adj)
{
   DEPRECATED("TrueAngle()","blTrueAngle()");
   return(blTrueAngle(opp, adj));
}

BOOL TorToCoor(VEC3F ant1, VEC3F ant2, VEC3F ant3, REAL bond, REAL theta,
               REAL torsion, VEC3F *coords)
{
   DEPRECATED("TorToCoor()","blTorToCoor()");
   return(blTorToCoor(ant1, ant2, ant3, bond, theta, torsion, coords));
}




/************************************************************************/
/* Renamed functions: array.h                                           */

char **Array2D(int size, int dim1, int dim2)
{
   DEPRECATED("Array2D()","blArray2D()");
   return(blArray2D(size, dim1, dim2));
}

void FreeArray2D(char **array, int dim1, int dim2)
{
   DEPRECATED("FreeArray2D()","blFreeArray2D()");
   blFreeArray2D(array, dim1, dim2);
}


char ***Array3D(int size, int dim1, int dim2, int dim3)
{
   DEPRECATED("Array3D()","blArray3D()");
   return(blArray3D(size, dim1, dim2, dim3));
}

void FreeArray3D(char ***array, int dim1, int dim2, int dim3)
{
   DEPRECATED("FreeArray3D()","blFreeArray3D()");
   blFreeArray3D(array, dim1, dim2, dim3);
}




/************************************************************************/
/* Renamed functions: help.h                                            */

void Help(char *string, char *HelpFile)
{
   DEPRECATED("Help()","blHelp()");
   blHelp(string, HelpFile);
}

void DoHelp(char *string, char *HelpFile)
{
   DEPRECATED("DoHelp()","blDoHelp()");
   blDoHelp(string, HelpFile);
}



/************************************************************************/
/* Renamed functions: hpgl.h                                            */

BOOL HPGLInit(char *filename, char *AltFont, REAL xmargin, REAL ymargin)
{
   DEPRECATED("HPGLInit()","blHPGLInit()");
   return(blHPGLInit(filename, AltFont, xmargin, ymargin));
}

void HPGLPen(int num)
{
   DEPRECATED("HPGLPen()","blHPGLPen()");
   blHPGLPen(num);
}

void HPGLMove(REAL x, REAL y)
{
   DEPRECATED("HPGLMove()","blHPGLMove()");
   blHPGLMove(x, y);
}

void HPGLDraw(REAL x, REAL y)
{
   DEPRECATED("HPGLDraw()","blHPGLDraw()");
   blHPGLDraw(x, y);
}

void HPGLSetDash(int style)
{
   DEPRECATED("HPGLSetDash()","blHPGLSetDash()");
   blHPGLSetDash(style);
}

void HPGLFont(int font, REAL size)
{
   DEPRECATED("HPGLFont()","blHPGLFont()");
   blHPGLFont(font, size);
}

void HPGLLText(REAL x, REAL y, char *string)
{
   DEPRECATED("HPGLLText()","blHPGLLText()");
   blHPGLLText(x, y, string);
}

void HPGLCBText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLCBText()","blHPGLCBText()");
   blHPGLCBText(x, y, offset, text);
}

void HPGLROffText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLROffText()","blHPGLROffText()");
   blHPGLROffText(x, y, offset, text);
}

void HPGLLCText(REAL x, REAL y, char *text)
{
   DEPRECATED("HPGLLCText()","blHPGLLCText()");
   blHPGLLCText(x, y, text);
}

void HPGLCTText(REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("HPGLCTText()","blHPGLCTText()");
   blHPGLCTText(x, y, offset, text);
}

void HPGLVText(REAL x, REAL y, REAL xoff, char *text, int TitleFont, 
               REAL TitleSize, char *label, int LabelFont, REAL LabelSize)
{
   DEPRECATED("HPGLVText()","blHPGLVText()");
   blHPGLVText(x, y, xoff, text, TitleFont, TitleSize, label, LabelFont, 
               LabelSize);
}

void HPGLEnd(void)
{
   DEPRECATED("HPGLEnd()","blHPGLEnd()");
   blHPGLEnd();
}

void HPGLShowText(char *text, BOOL orientation, int XBase, int YBase)
{
   DEPRECATED("HPGLShowText()","blHPGLShowText()");
   blHPGLShowText(text, orientation, XBase, YBase);
}



/************************************************************************/
/* Renamed functions: matrix.h                                          */

void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout)
{
   DEPRECATED("MatMult3_33()","blMatMult3_33()");
   blMatMult3_33(vecin, matin, vecout);
}

void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
{
   DEPRECATED("MatMult33_33()","blMatMult33_33()");
   blMatMult33_33(a, b, out);
}

void invert33(REAL s[3][3], REAL ss[3][3])
{
   DEPRECATED("invert33()","blInvert33()");
   blInvert33(s, ss);
}

void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
{
   DEPRECATED("CreateRotMat()","blCreateRotMat()");
   blCreateRotMat(direction, angle, matrix);
}

REAL VecDist(REAL *a, REAL *b, int len)
{
   DEPRECATED("VecDist()","blVecDist()");
   return(blVecDist(a, b, len));
}



/************************************************************************/
/* Renamed functions: parse.h                                           */

int parse(char *comline, int nkeys, KeyWd *keywords, REAL *REALparam,
          char **strparam)
{
   DEPRECATED("parse()","blParse()");
   return(blParse(comline, nkeys, keywords, REALparam, strparam));
}

int mparse(char *comline, int nkeys, MKeyWd *keywords, REAL *REALparam,
           char **strparam, int *nparams)
{
   DEPRECATED("mparse()","blMparse()");
   return(blMparse(comline, nkeys, keywords, REALparam, strparam, 
                   nparams));
}
           
int match(char *comstring, char *string2, int *nletters)
{
   DEPRECATED("match()","blMatch()");
   return(blMatch(comstring, string2, nletters));
}

int GetString(char *command, char *strparam)
{
   DEPRECATED("GetString()","blGetString()");
   return(blGetString(command, strparam));
}

int GetParam(char  *command, REAL *value, int *nletters)
{
   DEPRECATED("GetParam()","blGetParam()");
   return(blGetParam(command, value, nletters));
}



/************************************************************************/
/* Renamed functions: plotting.h                                        */

BOOL AMInitPlot(char *filename, char *title, int dest, REAL OutXSize, 
                REAL OutYSize, REAL OutXOff, REAL OutYOff,
                char *AltFont, REAL xmargin, REAL ymargin,
                REAL DataXMin, REAL DataYMin, REAL DataXMax,
                REAL DataYMax)
{
   DEPRECATED("AMInitPlot()","blAMInitPlot()");
   return(blAMInitPlot(filename, title, dest, OutXSize, OutYSize, OutXOff,
                       OutYOff, AltFont, xmargin, ymargin, DataXMin, 
                       DataYMin, DataXMax,DataYMax));
}

void AMSetPen(int dest, int pen)
{
   DEPRECATED("AMSetPen()","blAMSetPen()");
   blAMSetPen(dest, pen);
}

void AMMove(int dest, REAL x, REAL y)
{
   DEPRECATED("AMMove()","blAMMove()");
   blAMMove(dest, x, y);
}

void AMDraw(int dest, REAL x, REAL y)
{
   DEPRECATED("AMDraw()","blAMDraw()");
   blAMDraw(dest, x, y);
}

void AMSetLineStyle(int dest, int style)
{
   DEPRECATED("AMSetLineStyle()","blAMSetLineStyle()");
   blAMSetLineStyle(dest, style);
}

void AMEndLine(int dest)
{
   DEPRECATED("AMEndLine()","blAMEndLine()");
   blAMEndLine(dest);
}

void AMSetFont(int dest, char *PSFontName, REAL FontSize)
{
   DEPRECATED("AMSetFont()","blAMSetFont()");
   blAMSetFont(dest, PSFontName, FontSize);
}

void AMText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMText()","blAMText()");
   blAMText(dest, x, y, text);
}

void AMCBText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMCBText()","blAMCBText()");
   blAMCBText(dest, x, y, text);
}

void AMRText(int dest, REAL x, REAL y, REAL offset, char *text)
{
   DEPRECATED("AMRText()","blAMRText()");
   blAMRText(dest, x, y, offset, text);
}

void AMLCText(int dest, REAL x, REAL y, char *text)
{
   DEPRECATED("AMLCText()","blAMLCText()");
   blAMLCText(dest, x, y, text);
}

void AMCTText(int dest, REAL x, REAL y, REAL CTOffset, char *text)
{
   DEPRECATED("AMCTText()","blAMCTText()");
   blAMCTText(dest, x, y, CTOffset, text);
}

void AMEndPlot(int dest)
{
   DEPRECATED("AMEndPlot()","blAMEndPlot()");
   blAMEndPlot(dest);
}

int  PS2HPGLFont(char *font);
char *SimplifyText(char *string)
{
   DEPRECATED("SimplifyText()","blSimplifyText()");
   return(blSimplifyText(string));
}



/************************************************************************/
/* Renamed functions: ps.h                                              */

BOOL PSInit(char *FName, char *creator, char *AltFont)
{
   DEPRECATED("PSInit()","blPSInit()");
   return(blPSInit(FName, creator, AltFont));
}

void PSThick(REAL thickness)
{
   DEPRECATED("PSThick()","blPSThick()");
   blPSThick(thickness);
}

void PSMove(REAL X, REAL Y)
{
   DEPRECATED("PSMove()","blPSMove()");
   blPSMove(X, Y);
}

void PSDraw(REAL X, REAL Y)
{
   DEPRECATED("PSDraw()","blPSDraw()");
   blPSDraw(X, Y);
}

void PSSetDash(char *linepatt)
{
   DEPRECATED("PSSetDash()","blPSSetDash()");
   blPSSetDash(linepatt);
}

void PSClearDash(void)
{
   DEPRECATED("PSClearDash()","blPSClearDash()");
   blPSClearDash();
}

void PSStroke(void)
{
   DEPRECATED("PSStroke()","blPSStroke()");
   blPSStroke();
}

void PSFont(char *fontname, REAL size)
{
   DEPRECATED("PSFont()","blPSFont()");
   blPSFont(fontname, size);
}

void PSLText(REAL X, REAL Y, char *label)
{
   DEPRECATED("PSLText()","blPSLText()");
   blPSLText(X, Y, label);
}

void PSCBText(REAL X, REAL Y, REAL Offset, char *label)
{
   DEPRECATED("PSCBText()","blPSCBText()");
   blPSCBText(X, Y, Offset, label);
}

void PSROffText(REAL X, REAL Y, REAL offset, char *label)
{
   DEPRECATED("PSROffText()","blPSROffText()");
   blPSROffText(X, Y, offset, label);
}

void PSLCText(REAL X, REAL Y, char *label)
{
   DEPRECATED("PSLCText()","blPSLCText()");
   blPSLCText(X, Y, label);
}

void PSCTText(REAL X, REAL Y, REAL Offset, char *label)
{
   DEPRECATED("PSCTText()","blPSCTText()");
   blPSCTText(X, Y, Offset, label);
}

void PSVText(REAL x, REAL y, REAL xoff, char *text, char *font, REAL size,
             char *label, char *lfont, REAL lsize)
{
   DEPRECATED("PSVText()","blPSVText()");
   blPSVText(x, y, xoff, text, font, size, label, lfont, lsize);
}

void PSShowText(char *text)
{
   DEPRECATED("PSShowText()","blPSShowText()");
   blPSShowText(text);
}

void PSEnd(void)
{
   DEPRECATED("PSEnd()","blPSEnd()");
   blPSEnd();
}

char *PSCorrectCase(char *font)
{
   DEPRECATED("PSCorrectCase()","blPSCorrectCase()");
   return(blPSCorrectCase(font));
}



/************************************************************************/
/* Renamed functions: safemem.h                                         */

void *safemalloc(int nbytes)
{
   DEPRECATED("safemalloc()","blSafemalloc()");
   return(blSafemalloc(nbytes));
}

BOOL safefree(void *ptr)
{
   DEPRECATED("safefree()","blSafefree()");
   return(blSafefree(ptr));
}

void safeleaks(void)
{
   DEPRECATED("safeleaks()","blSafeleaks()");
   blSafeleaks();
}

/************************************************************************/
/** \endcond                                                            */
/************************************************************************/

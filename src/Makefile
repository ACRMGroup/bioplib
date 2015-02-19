# makefile for building support libraries under UNIX:

# Destination in which to place the library and the include files
LIBDEST       = ${HOME}/lib
INCDEST       = ${HOME}/include
DATADEST      = ${HOME}/data
SHAREDLIBDEST = ${HOME}/lib

# Doxygen and adddoxy directories
DOCDIR      = ../doc
DOXYDIR     = $(DOCDIR)/doxygen
ADDDOXY     = $(DOCDIR)/adddoxy/adddoxy.pl
ADDDOXYCONF = $(DOCDIR)/adddoxy/adddoxy.conf

# Define the C compiler here
CC = gcc
 
# Define any options for the C compiler here
COPT = -fPIC -ansi -Wall -pedantic -g -DGUNZIP_SUPPORT

# Generate error messages for deprecated functions.
# Default is to display deprecated error message unless environment 
# variable BIOPLIB_DEPRECATED_QUIET is set.
# Uncomment compile options below to ignore the environment and either 
# always display the error message or silence the error message:
#COPT := $(COPT) -D BIOPLIB_DEPRECATED_CHECK
#COPT := $(COPT) -D BIOPLIB_DEPRECATED_QUIET

# Include libxml2
# Note: xml2-config is installed with libxml2
# When you compile code you need to link to libxml2 with -lxml2
# Comment out this line if you do not require PDBML (XML) support
COPT := $(COPT) -D XML_SUPPORT $(shell xml2-config --cflags)

# Use single letter check for filetype
# Only check first character of file when detecting file type (compressed
# file or pdbml).
# For systems were ungetc() fails after pushing back a single character to
# the input stream.
# Note: This option is required for MS Windows.
#COPT := $(COPT) -D SINGLE_CHAR_FILECHECK

# Define the archive/library command here
AR = ar r

# If ranlib required when creating libraries, set here; 
# otherwise set to echo
RANLIB = ranlib
#RANLIB = echo

###########################################################
#          DON'T CHANGE ANYTHING BELOW THIS POINT         #
###########################################################

# Version numbers for shared libraries - strictly the MINOR numbers here
# are the minor and release numbers
# libgen
GMAJOR = 3
GMINOR = 0.3
# libbiop
BMAJOR = 3
BMINOR = 0.3

# Files for libgen.a
OFILESG = BuffInp.o CalcExtSD.o CalcSD.o CheckExtn.o CreateRotMat.o CrossProd3.o \
DistPtVect.o ErrStack.o FreeStringList.o GetFilestem.o GetWord.o \
InStringList.o IndxReal.o KillLeadSpaces.o KillLine.o MatMult33_33.o \
MatMult3_33.o NComb.o NPerm.o OpenFile.o OpenStdFiles.o QueryStrStr.o \
SetExtn.o StoreString.o StringToLower.o StringToUpper.o TorCoor.o TrueAngle.o \
VecAdd3.o VecDist.o VecLen3.o VecSub3.o Word.o WordN.o angle.o \
array2.o chindex.o countchar.o factdiv.o factorial.o fgetsany.o \
fsscanf.o ftostr.o help.o hpgl.o invert33.o openorpipe.o padchar.o \
padterm.o parse.o pearson.o pearson1.o phi.o pldist.o plotting.o \
ps.o safemem.o simpleangle.o strcatalloc.o upstrcmp.o upstrncmp.o \
WindIO.o getfield.o array3.o justify.o wrapprint.o deprecatedGen.o \
eigen.o regression.o


# Files for libbiop.a
OFILESB = AddNTerHs.o ApMatPDB.o AppendPDB.o AtomNameMatch.o \
CalcCellTrans.o CalcChiPDB.o CalcCterCoords.o CalcRMSPDB.o \
CalcTetraHCoords.o CopyPDB.o CopyPDBCoords.o \
DNAtoAA.o DupePDB.o ExtractZonePDB.o FNam2PDB.o FindAtomInRes.o \
FindNextChainPDB.o FindNextResidue.o FindResidue.o FindResidueSpec.o \
FindZonePDB.o FitCaCbPDB.o FitCaPDB.o FitPDB.o FixCterPDB.o \
GetCGPDB.o GetCofGPDBRange.o GetCofGPDBSCRange.o GetCrystPDB.o \
GetPDBByN.o GetPDBChainLabels.o GetPDBCoor.o HAddPDB.o InPDBZone.o \
InPDBZoneSpec.o IndexPDB.o KillPDB.o KillSidechain.o KnownSeqLen.o \
LegalAtomSpec.o MovePDB.o NumericAlign.o OrderPDB.o OriginPDB.o \
PDB2Seq.o ParseRes.o RdSSPDB.o RdSecPDB.o RdSeqPDB.o \
ReadCSSR.o ReadPDB.o ReadPIR.o ReadRawPIR.o ReadSimplePIR.o \
RenumAtomsPDB.o ResolPDB.o RotPDB.o SelAtPDB.o SelectCaPDB.o \
SetChi.o SetResnam.o SplitSeq.o StripHPDB.o TermPDB.o \
TranslatePDB.o TrueSeqLen.o WriteCrystPDB.o WritePDB.o WrtCSSR.o \
align.o fit.o hbond.o rsc.o throne.o FitNCaCPDB.o FindRawAtomInRes.o \
WriteGromosPDB.o WholePDB.o GlyCB.o BuildAtomNeighbourPDBList.o \
FindAtomWildcardInRes.o DupeResiduePDB.o StripWatersPDB.o aalist.o \
StructurePDB.o FindHetatmResidue.o FindHetatmResidueSpec.o access.o \
deprecatedBiop.o BuildConect.o


# Static libraries - the default
all : libgen.a libbiop.a 

libgen.a : $(OFILESG)
	$(AR) libgen.a $? 
	$(RANLIB) libgen.a
libbiop.a : $(OFILESB)
	$(AR) libbiop.a $? 
	$(RANLIB) libbiop.a

# Shared libraries
shared : libgens.so.$(GMAJOR).$(GMINOR) libbiops.so.$(BMAJOR).$(BMINOR)

libgens.so.$(GMAJOR).$(GMINOR) : $(OFILESG)
	$(CC) -shared -fPIC -Wl,-soname,libgens.so.$(GMAJOR) -o libgens.so.$(GMAJOR).$(GMINOR) $? -lc
libbiops.so.$(BMAJOR).$(BMINOR) : $(OFILESB)
	$(CC) -shared -fPIC -Wl,-soname,libbiops.so.$(BMAJOR) -o libbiops.so.$(BMAJOR).$(BMINOR) $? -lc

# C compilation
.c.o : 
	$(CC) $(COPT) -o $@ -c $<

# Cleanup
clean :
	rm -f $(OFILESB) $(OFILESG) libgen.a libbiop.a libgens.so.$(GMAJOR).$(GMINOR) libbiops.so.$(BMAJOR).$(BMINOR)

# Documentation
doxygen :
	$(ADDDOXY) -conf=$(ADDDOXYCONF) *.c *.h
	mkdir -p $(DOXYDIR)/docsrc
	cp $(DOXYDIR)/docsrcinput/*.dox $(DOXYDIR)/docsrc
	(cd $(DOXYDIR); doxygen)
doxyclean :
	rm -rf $(DOCDIR)/html $(DOCDIR)/latex $(DOXYDIR)/docsrc


# Installation
install :
	mkdir -p $(LIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp libgen.a $(LIBDEST)
	cp libbiop.a $(LIBDEST)
	cp *.h $(INCDEST)/bioplib

installdata :
	mkdir -p $(DATADEST)
	(cd ../data; cp -Rcp * $(DATADEST))

installshared : 
	mkdir -p $(SHAREDLIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp *.h $(INCDEST)/bioplib
	cp libgens.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	cp libbiops.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	(cd $(SHAREDLIBDEST); ln -s libgens.so.$(GMAJOR).$(GMINOR) libgens.so.$(GMAJOR); ln -s libgens.so.$(GMAJOR) libgens.so)
	(cd $(SHAREDLIBDEST); ln -s libbiops.so.$(BMAJOR).$(BMINOR) libbiops.so.$(BMAJOR); ln -s libbiops.so.$(BMAJOR) libbiops.so)

# Help
help :
	@echo "make               : Build the static libraries"
	@echo "make install       : Install the static libraries"
	@echo "make installdata   : Install the data files"
	@echo "make shared        : Build shared libraries"
	@echo "make installshared : Install shared libraries"
	@echo "make clean         : Remove object and library files from compilation directory"
	@echo "make doxygen       : Generate doxygen documentation"

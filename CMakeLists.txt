cmake_minimum_required(VERSION 3.1)
project(bioplib)


###############################################################################
#                             PATH DEFINITIONS
###############################################################################

set(DEST $ENV{HOME} CACHE PATH "Installation directory" )
set(LIBDEST ${DEST}/lib CACHE PATH "Installation directory for static libraries")
set(INCDEST ${DEST}/include CACHE PATH "Installation directory for header files")
set(DATADEST ${DEST}/data CACHE PATH "Installation directory for data")
set(SHAREDLIBDEST ${DEST}/lib CACHE PATH "Installation directory for shared libraries")
set(DOCDIR ${PROJECT_SOURCE_DIR}/doc CACHE PATH "Installation directory for documentation")
set(DOXYDIR ${DOCDIR}/doxygen CACHE PATH "Installation directory for doxy documentation")
set(ADDDOXY ${DOCDIR}/adddoxy/adddoxy.pl CACHE PATH "Location of adddoxy.pl script")
set(ADDDOXYCONF ${DOCDIR}/adddoxy/adddoxy.conf CACHE PATH "Location of adddoxy.conf configuration file")

message("-- Installation directory                       : ${DEST}")
message("-- Installation directory for static libraries  : ${LIBDEST}")
message("-- Installation directory for header files      : ${INCDEST}")
message("-- Installation directory for data              : ${DATADEST}")
message("-- Installation directory for shared libraries  : ${SHAREDLIBDEST}")
message("-- Installation directory for documentation     : ${DOCDIR}")
message("-- Installation directory for doxy documentation: ${DOXYDIR}")
message("-- Location of adddoxy.pl script                : ${ADDDOXY}")
message("-- Location of adddoxy.conf configuration file  : ${ADDDOXYCONF}")

###############################################################################
#                               FIND PACKAGES
###############################################################################
find_package(LibXml2 REQUIRED)

include_directories(${LIBXML2_INCLUDE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CheckCSourceCompiles)
include(CheckCSourceRuns)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckSymbolExists)
include(CheckTypeSize)

###############################################################################
#                               CHECK HEADERS
###############################################################################

set(INCLUDES "")
macro(ck_check_include_file header var)
    check_include_files("${INCLUDES};${header}" ${var})
    if(${var})
        set(INCLUDES ${INCLUDES} ${header})
    endif(${var})
endmacro(ck_check_include_file)

ck_check_include_file("stdlib.h" HAVE_STDLIB_H)

###############################################################################
#                           CHECK FOR INTEGER TYPES
###############################################################################

check_type_size(intmax_t INTMAX_T)
check_type_size(uintmax_t UINTMAX_T)

check_type_size(pid_t PID_T)
if(NOT HAVE_PID_T)
    if(WIN32)
        set(pid_t "int")
    else(WIN32)
        MESSAGE(FATAL_ERROR "pid_t doesn't exist on this platform?")
    endif(WIN32)
endif(NOT HAVE_PID_T)


###############################################################################

# Generate "config.h" from "cmake/config.h.cmake"
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in
        ${CMAKE_CURRENT_BINARY_DIR}/config.h)
include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-DHAVE_CONFIG_H)
set(CONFIG_HEADER ${CMAKE_CURRENT_BINARY_DIR}/config.h)


###############################################################################


###############################################################################
#                            VERSION DEFINITIONS
###############################################################################
# libgen
set(GMAJOR 3)
set(GMINOR 5.0)
# libbiop
set(BMAJOR 3)
set(BMINOR 5.0)

###############################################################################
#                           COMPILER AND LINKER FLAGS
###############################################################################

option(XML_SUPPORT ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -ansi -Wall -pedantic -g -DGUNZIP_SUPPORT -finline-functions")

if (${LIBXML2_FOUND})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D XML_SUPPORT")
else()
    message("-- found xml2: false")
endif ()

###############################################################################
#                         DEFINE SOURCE FILES
###############################################################################

set(GEN_SOURCES src/BuffInp.c
        src/CalcExtSD.c
        src/CalcSD.c
        src/CheckExtn.c
        src/CreateRotMat.c
        src/CrossProd3.c
        src/DistPtLine.c
        src/ErrStack.c
        src/FreeStringList.c
        src/GetFilestem.c
        src/GetWord.c
        src/InStringList.c
        src/IndxReal.c
        src/KillLeadSpaces.c
        src/KillLine.c
        src/MatMult33_33.c
        src/MatMult3_33.c
        src/NComb.c
        src/NPerm.c
        src/OpenFile.c
        src/OpenStdFiles.c
        src/QueryStrStr.c
        src/SetExtn.c
        src/StoreString.c
        src/StringToLower.c
        src/StringToUpper.c
        src/TorCoor.c
        src/TrueAngle.c
        src/VecAdd3.c
        src/VecDist.c
        src/VecLen3.c
        src/VecSub3.c
        src/Word.c
        src/WordN.c
        src/angle.c
        src/array2.c
        src/chindex.c
        src/countchar.c
        src/factdiv.c
        src/factorial.c
        src/fgetsany.c
        src/fsscanf.c
        src/ftostr.c
        src/help.c
        src/hpgl.c
        src/invert33.c
        src/openorpipe.c
        src/padchar.c
        src/padterm.c
        src/parse.c
        src/pearson.c
        src/pearson1.c
        src/phi.c
        src/pldist.c
        src/plotting.c
        src/ps.c
        src/safemem.c
        src/simpleangle.c
        src/strcatalloc.c
        src/upstrcmp.c
        src/upstrncmp.c
        src/WindIO.c
        src/getfield.c
        src/array3.c
        src/justify.c
        src/wrapprint.c
        src/deprecatedGen.c
        src/eigen.c
        src/regression.c
        src/filename.c
        src/stringcat.c
        src/stringutil.c
        src/hash.c
        src/prime.c
        )

set(BIOP_SOURCES
        src/00bioplib.c
        src/aalist.c
        src/access.c
        src/AddNTerHs.c
        src/align.c
        src/align.c
        src/ApMatPDB.c
        src/AppendPDB.c
        src/AtomNameMatch.c
        src/atomtype.c
        src/BuildAtomNeighbourPDBList.c
        src/BuildConect.c
        src/CalcCellTrans.c
        src/CalcChiPDB.c
        src/CalcCterCoords.c
        src/CalcRMSPDB.c
        src/CalcTetraHCoords.c
        src/CopyPDB.c
        src/CopyPDBCoords.c
        src/deprecatedBiop.c
        src/DNAtoAA.c
        src/DupePDB.c
        src/DupeResiduePDB.c
        src/ExtractZonePDB.c
        src/FindAtomInRes.c
        src/FindHetatmResidue.c
        src/FindAtomWildcardInRes.c
        src/FindHetatmResidueSpec.c
        src/FindNextChainPDB.c
        src/FindNextResidue.c
        src/FindRawAtomInRes.c
        src/FindResidue.c
        src/FindResidueSpec.c
        src/FindZonePDB.c
        src/fit.c
        src/FitCaCbPDB.c
        src/FitCaPDB.c
        src/FitNCaCPDB.c
        src/FitPDB.c
        src/FixCterPDB.c
        src/FNam2PDB.c
        src/GetCGPDB.c
        src/GetCofGPDBRange.c
        src/GetCofGPDBSCRange.c
        src/GetCrystPDB.c
        src/GetPDBByN.c
        src/GetPDBChainAsCopy.c
        src/GetPDBChainLabels.c
        src/GetPDBCoor.c
        src/GlyCB.c
        src/HAddPDB.c
        src/hbond.c
        src/IndexPDB.c
        src/InPDBZone.c
        src/InPDBZoneSpec.c
        src/KillPDB.c
        src/KillSidechain.c
        src/KnownSeqLen.c
        src/LegalAtomSpec.c
        src/MovePDB.c
        src/NumericAlign.c
        src/OrderPDB.c
        src/ParseRes.c
        src/PDB2Seq.c
        src/PDBHeaderInfo.c
        src/RdSecPDB.c
        src/RdSeqPDB.c
        src/RdSSPDB.c
        src/ReadCSSR.c
        src/ReadPDB.c
        src/ReadPIR.c
        src/ReadRawPIR.c
        src/ReadSimplePIR.c
        src/RenumAtomsPDB.c
        src/ResolPDB.c
        src/RotPDB.c
        src/rsc.c
        src/secstr.c
        src/SelAtPDB.c
        src/SelectCaPDB.c
        src/SetChi.c
        src/SetResnam.c
        src/SplitSeq.c
        src/StripHPDB.c
        src/StripWatersPDB.c
        src/StructurePDB.c
        src/TermPDB.c
        src/test.c
        src/throne.c
        src/TranslatePDB.c
        src/TrueSeqLen.c
        src/WriteCrystPDB.c
        src/WritePDB.c
        src/WritePIR.c
        src/WrtCSSR.c
        )

###############################################################################
#                          DEFINE HEADER FILES
###############################################################################
set(HEADER_FILES src/array.h
        src/cssr.h
        src/deprecated.h
        src/general.h
        src/macros.h
        src/MathType.h
        src/MathUtil.h
        src/matrix.h
        src/OriginPDB.c
        src/pdb.h
        src/port.h
        src/seq.h
        src/SysDefs.h)

###############################################################################
#                   DEFINE STATIC AND DYNAMIC LIBRARIES
###############################################################################

# build libraries libgen and libbiop
add_library(gen OBJECT ${GEN_SOURCES} ${HEADER_FILES})
add_library(biop OBJECT ${BIOP_SOURCES} ${HEADER_FILES})

# build static libraries
add_library(gen-static STATIC  $<TARGET_OBJECTS:gen>)
add_library(biop-static STATIC $<TARGET_OBJECTS:biop>)
target_link_libraries(biop-static ${LIBXML2_LIBRARIES})

# build dynamics libraries
add_library(gen-shared SHARED $<TARGET_OBJECTS:gen>)
add_library(biop-shared SHARED $<TARGET_OBJECTS:biop>)
target_link_libraries(biop-shared ${LIBXML2_LIBRARIES})

set_target_properties(gen-static PROPERTIES OUTPUT_NAME "libgen")
set_target_properties(biop-static PROPERTIES OUTPUT_NAME "libbiop")

set_target_properties(gen-shared
        PROPERTIES
            OUTPUT_NAME "gens"
            VERSION ${GMAJOR}.${GMINOR})
set_target_properties(biop-shared
        PROPERTIES
            OUTPUT_NAME "biops"
            VERSION ${BMAJOR}.${BMINOR})

###############################################################################
#                           MAKE RULES
###############################################################################

# install rules
install(TARGETS gen-static LIBRARY ARCHIVE DESTINATION ${LIBDEST})
install(TARGETS biop-static LIBRARY ARCHIVE DESTINATION ${LIBDEST})

# shared rules
add_custom_target(shared
        COMMAND cp ${PROJECT_BINARY_DIR}/libgen.so.${GMAJOR}.${GMINOR} ${LIBDEST}
        COMMAND cp ${PROJECT_BINARY_DIR}/libbiop.so.${BMAJOR}.${BMINOR} ${LIBDEST})


# installdata rules
add_custom_target(installdata
        COMMAND mkdir -p ${DATADEST}
        COMMAND cp -Rp data/* ${DATADEST}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

# doxygen rules
add_custom_target(doxygen
        COMMAND ${ADDDOXY} -conf=${ADDDOXYCONF} *.c *.h
        COMMAND mkdir -p $(DOXYDIR)/docsrc
        COMMAND cp ${DOXYDIR}/docsrcinput/*.dox ${DOXYDIR}/docsrc
        COMMAND cd ${DOXYDIR}
        COMMAND doxygen)

# doxyclean rules
add_custom_target(doxyclean
        COMMAND rm -rf ${DOCDIR}/html ${DOCDIR}/latex ${DOXYDIR}/docsrc)

###############################################################################
#                                  TESTS
###############################################################################

add_subdirectory(test)

enable_testing()
add_test(NAME BioplibTESTS
        COMMAND bioplib_tests
        WORKING_DIRECTORY ${bioplibTests_SOURCE_DIR})

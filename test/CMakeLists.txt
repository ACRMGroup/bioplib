cmake_minimum_required(VERSION 3.1)
project(bioplibTests)

###############################################################################
#                               FIND PACKAGES
###############################################################################

find_package(LibXml2 REQUIRED)
find_package(Check REQUIRED)
find_package(Threads REQUIRED)

include_directories(${LIBXML2_INCLUDE_DIR}
        ${CHECK_INCLUDE_DIRS}
        ${bioplibTests_SOURCE_DIR}/src
        ${bioplib_SOURCE_DIR}/src)

###############################################################################
#                         DEFINE SOURCE FILES
###############################################################################

set(TestSuite_SOURCES
        src/conect_suite.c
        src/findzone_suite.c
        src/getpdbchainlabels_suite.c
        src/header_suite.c
        src/inpdbzone_suite.c
        src/parseresspec_suite.c
        src/readpdbml_suite.c
        src/wholepdb_suite.c
        src/writepdbml_suite.c)

set(TestSuite_HEADERS
        src/connect_suite.h
        src/findzone_suite.h
        src/getpdbchainlabels_suite.h
        src/header_suite.h
        src/inpdbzone_suite.h
        src/parseresspec_suite.h
        src/readpdbml_suite.h
        src/wholepdb_suite.h
        src/writepdbml_suite.h)

set(TARGET_NAME bioplib_tests)

###############################################################################
#                                BUILD TESTS
###############################################################################


add_executable(bioplib_tests src/main.c ${TestSuite_SOURCES})

add_custom_command(TARGET bioplib_tests PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${bioplibTests_SOURCE_DIR}/data  ${bioplibTests_BINARY_DIR}/data)

target_link_libraries(bioplib_tests PRIVATE
        biop-shared
        gen-shared
        ${CHECK_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})